#!/usr/bin/env python3
"""Local ReplicaLab server: static GUI, project storage, DCD audit, and RMSD jobs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import statistics
import struct
import subprocess
import threading
import time
import traceback
import uuid
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PROJECTS_FILE = DATA / "projects.json"
JOBS_FILE = DATA / "jobs.json"
RESULTS = ROOT / "local-results"
LOCK = threading.RLock()
VMD = shutil.which("vmd")


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path, default):
    with LOCK:
        if not path.exists():
            return default
        try:
            return json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            return default


def write_json(path: Path, value) -> None:
    with LOCK:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(value, indent=2) + "\n")
        temporary.replace(path)


def safe_resolve(raw: str, expected: str | None = None) -> Path:
    if not raw or "\x00" in raw:
        raise ValueError("A non-empty local path is required.")
    path = Path(raw).expanduser().resolve()
    if expected == "file" and not path.is_file():
        raise ValueError(f"File not found: {path}")
    if expected == "directory" and not path.is_dir():
        raise ValueError(f"Directory not found: {path}")
    return path


def natural_key(path: Path):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", path.name)]


def dcd_header(path: Path) -> dict:
    """Read standard CHARMM/NAMD DCD header fields without external libraries."""
    result = {"frames": None, "start_step": None, "save_frequency": None}
    try:
        with path.open("rb") as handle:
            raw = handle.read(96)
        if len(raw) < 24:
            return result
        marker_le = struct.unpack("<i", raw[:4])[0]
        endian = "<" if marker_le in (84, 164) else ">"
        marker = struct.unpack(endian + "i", raw[:4])[0]
        if marker not in (84, 164) or raw[4:8] not in (b"CORD", b"VELD"):
            return result
        result["frames"] = struct.unpack(endian + "i", raw[8:12])[0]
        result["start_step"] = struct.unpack(endian + "i", raw[12:16])[0]
        result["save_frequency"] = struct.unpack(endian + "i", raw[16:20])[0]
    except (OSError, struct.error):
        pass
    return result


def validate_replica(replica: dict) -> dict:
    directory = safe_resolve(replica.get("trajectory_dir", ""), "directory")
    pattern = replica.get("pattern", "*.dcd").strip() or "*.dcd"
    if "/" in pattern or "\\" in pattern or pattern.startswith("."):
        raise ValueError("Trajectory pattern must be a filename glob such as md-*.dcd.")
    files = sorted((p for p in directory.glob(pattern) if p.is_file()), key=natural_key)
    if not files:
        raise ValueError(f"No trajectories matching {pattern!r} were found in {directory}.")
    try:
        interval_ns = float(replica.get("frame_interval_ns", 0))
    except (TypeError, ValueError):
        interval_ns = 0
    if interval_ns <= 0:
        raise ValueError("Frame interval must be greater than zero.")
    sizes = [path.stat().st_size for path in files]
    median_size = statistics.median(sizes)
    segments = []
    total_frames = 0
    readable_frames = True
    for path, size in zip(files, sizes):
        header = dcd_header(path)
        frames = header["frames"]
        if frames is None or frames < 0:
            readable_frames = False
        else:
            total_frames += frames
        ratio = size / median_size if median_size else 1
        segments.append({
            "path": str(path),
            "name": path.name,
            "bytes": size,
            "frames": frames,
            "start_step": header["start_step"],
            "save_frequency": header["save_frequency"],
            "size_ratio": round(ratio, 3),
            "warning": ratio < 0.65,
        })
    return {
        "id": replica.get("id") or f"rep_{uuid.uuid4().hex[:10]}",
        "name": replica.get("name", "Replica").strip() or "Replica",
        "trajectory_dir": str(directory),
        "pattern": pattern,
        "frame_interval_ns": interval_ns,
        "start_time_ns": float(replica.get("start_time_ns", 0) or 0),
        "segments": segments,
        "segment_count": len(segments),
        "total_frames": total_frames if readable_frames else None,
        "estimated_duration_ns": (
            round(total_frames * interval_ns, 6) if readable_frames else None
        ),
        "warnings": sum(item["warning"] for item in segments),
        "validated_at": utcnow(),
    }


def validate_system(payload: dict) -> dict:
    psf = safe_resolve(payload.get("psf", ""), "file")
    pdb = safe_resolve(payload.get("pdb", ""), "file")
    name = str(payload.get("name", "")).strip()
    if not name:
        raise ValueError("System name is required.")
    replicas = payload.get("replicas") or []
    if not replicas:
        raise ValueError("At least one replica is required.")
    validated = [validate_replica(replica) for replica in replicas]
    return {
        "id": payload.get("id") or f"sys_{uuid.uuid4().hex[:10]}",
        "name": name,
        "description": str(payload.get("description", "")).strip(),
        "psf": str(psf),
        "pdb": str(pdb),
        "ligand_selection": str(payload.get("ligand_selection", "resname LIG")).strip(),
        "replicas": validated,
        "validated_at": utcnow(),
    }


def update_job(job_id: str, **changes) -> dict:
    jobs = read_json(JOBS_FILE, [])
    for job in jobs:
        if job["id"] == job_id:
            job.update(changes)
            write_json(JOBS_FILE, jobs)
            return job
    raise KeyError(job_id)


def append_job(job: dict) -> None:
    jobs = read_json(JOBS_FILE, [])
    jobs.insert(0, job)
    write_json(JOBS_FILE, jobs[:200])


def make_manifest(replica: dict, destination: Path, stride: int) -> None:
    start = float(replica.get("start_time_ns", 0))
    interval = float(replica["frame_interval_ns"])
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "start_ns", "interval_ns", "stride"])
        cumulative = start
        for segment in replica["segments"]:
            writer.writerow([segment["path"], cumulative, interval, stride])
            frames = segment.get("frames")
            if frames:
                cumulative += frames * interval


def summarize_rmsd(csv_path: Path, output_path: Path) -> dict:
    columns: dict[str, list[float]] = {}
    frames = 0
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            frames += 1
            for key, value in row.items():
                if key in {"time_ns", "segment", "frame"}:
                    continue
                try:
                    columns.setdefault(key, []).append(float(value))
                except (TypeError, ValueError):
                    continue
    summary = {"frames": frames, "quantities": {}}
    for name, values in columns.items():
        if values:
            summary["quantities"][name] = {
                "mean": statistics.fmean(values),
                "sd": statistics.stdev(values) if len(values) > 1 else 0,
                "min": min(values),
                "max": max(values),
            }
    write_json(output_path, summary)
    return summary


def run_rmsd_job(job_id: str, system: dict, replicas: list[dict], config: dict) -> None:
    try:
        if not VMD:
            raise RuntimeError("VMD is not installed or is not available on PATH.")
        output_root = RESULTS / job_id
        output_root.mkdir(parents=True, exist_ok=True)
        selections = config.get("selections") or [
            {"name": "protein_backbone", "selection": "protein and backbone"}
        ]
        selection_spec = ";".join(
            f"{item['name']}|{item['selection']}" for item in selections
        )
        completed = 0
        outputs = []
        for replica in replicas:
            update_job(
                job_id,
                status="running",
                message=f"Analyzing {system['name']} · {replica['name']}",
                progress=round(100 * completed / len(replicas)),
            )
            replica_dir = output_root / replica["id"]
            replica_dir.mkdir(parents=True, exist_ok=True)
            manifest = replica_dir / "trajectory_manifest.csv"
            make_manifest(replica, manifest, max(1, int(config.get("stride", 1))))
            output_csv = replica_dir / "rmsd.csv"
            log_path = replica_dir / "vmd.log"
            environment = os.environ.copy()
            environment.update({
                "RL_PSF": system["psf"],
                "RL_PDB": system["pdb"],
                "RL_MANIFEST": str(manifest),
                "RL_OUTPUT": str(output_csv),
                "RL_ALIGNMENT": config.get(
                    "alignment_selection", "protein and name CA"
                ),
                "RL_SELECTIONS": selection_spec,
            })
            with log_path.open("w") as log:
                process = subprocess.run(
                    [VMD, "-dispdev", "text", "-e", str(ROOT / "rmsd_runner.tcl")],
                    cwd=replica_dir,
                    env=environment,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if process.returncode != 0 or not output_csv.exists():
                raise RuntimeError(
                    f"VMD failed for {replica['name']}; see {log_path}"
                )
            summary_path = replica_dir / "summary.json"
            summary = summarize_rmsd(output_csv, summary_path)
            outputs.append({
                "replica_id": replica["id"],
                "replica_name": replica["name"],
                "csv": str(output_csv),
                "summary": str(summary_path),
                "frames": summary["frames"],
                "log": str(log_path),
            })
            completed += 1
            update_job(job_id, progress=round(100 * completed / len(replicas)))
        update_job(
            job_id,
            status="completed",
            progress=100,
            message="RMSD analysis completed",
            outputs=outputs,
            completed_at=utcnow(),
        )
    except Exception as error:
        update_job(
            job_id,
            status="failed",
            message=str(error),
            error=traceback.format_exc(),
            completed_at=utcnow(),
        )


class ReplicaLabHandler(SimpleHTTPRequestHandler):
    server_version = "ReplicaLab/0.2"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def log_message(self, fmt, *args):
        print(f"[{self.log_date_time_string()}] {fmt % args}")

    def json_response(self, value, status=HTTPStatus.OK):
        body = json.dumps(value).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def request_json(self):
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            length = 0
        if length > 2_000_000:
            raise ValueError("Request is too large.")
        raw = self.rfile.read(length) if length else b"{}"
        return json.loads(raw)

    def do_GET(self):
        route = urlparse(self.path).path
        if route == "/api/status":
            return self.json_response({
                "ok": True,
                "version": "0.2.0",
                "vmd_available": bool(VMD),
                "vmd_path": VMD,
                "systems": len(read_json(PROJECTS_FILE, [])),
                "jobs": len(read_json(JOBS_FILE, [])),
            })
        if route == "/api/systems":
            return self.json_response({"systems": read_json(PROJECTS_FILE, [])})
        if route == "/api/jobs":
            return self.json_response({"jobs": read_json(JOBS_FILE, [])})
        return super().do_GET()

    def do_POST(self):
        route = urlparse(self.path).path
        try:
            payload = self.request_json()
            if route == "/api/systems/validate":
                return self.json_response({"system": validate_system(payload)})
            if route == "/api/systems":
                system = validate_system(payload)
                systems = read_json(PROJECTS_FILE, [])
                systems = [item for item in systems if item["id"] != system["id"]]
                systems.append(system)
                write_json(PROJECTS_FILE, systems)
                return self.json_response({"system": system}, HTTPStatus.CREATED)
            if route == "/api/rmsd/run":
                systems = read_json(PROJECTS_FILE, [])
                system = next(
                    (item for item in systems if item["id"] == payload.get("system_id")),
                    None,
                )
                if not system:
                    raise ValueError("Choose a saved local system.")
                requested = set(payload.get("replica_ids") or [])
                replicas = [
                    item for item in system["replicas"]
                    if not requested or item["id"] in requested
                ]
                if not replicas:
                    raise ValueError("Choose at least one replica.")
                job = {
                    "id": f"job_{uuid.uuid4().hex[:12]}",
                    "type": "rmsd",
                    "system_id": system["id"],
                    "system_name": system["name"],
                    "replica_ids": [item["id"] for item in replicas],
                    "status": "queued",
                    "progress": 0,
                    "message": "Waiting for local runner",
                    "config": payload.get("config") or {},
                    "created_at": utcnow(),
                }
                append_job(job)
                thread = threading.Thread(
                    target=run_rmsd_job,
                    args=(job["id"], system, replicas, job["config"]),
                    daemon=True,
                )
                thread.start()
                return self.json_response({"job": job}, HTTPStatus.ACCEPTED)
            return self.json_response(
                {"error": "Unknown API route"}, HTTPStatus.NOT_FOUND
            )
        except (ValueError, json.JSONDecodeError) as error:
            return self.json_response({"error": str(error)}, HTTPStatus.BAD_REQUEST)
        except Exception as error:
            return self.json_response(
                {"error": str(error)}, HTTPStatus.INTERNAL_SERVER_ERROR
            )


def main():
    parser = argparse.ArgumentParser(description="Run ReplicaLab locally.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8766)
    args = parser.parse_args()
    DATA.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    server = ThreadingHTTPServer((args.host, args.port), ReplicaLabHandler)
    print(f"ReplicaLab is running at http://{args.host}:{args.port}")
    print(f"VMD adapter: {'ready at ' + VMD if VMD else 'not available'}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
