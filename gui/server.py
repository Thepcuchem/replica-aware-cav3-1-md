#!/usr/bin/env python3
"""ReplicaLab local web server and safe analysis-job runner."""

from __future__ import annotations

import argparse, json, mimetypes, os, shutil, subprocess, sys, threading, traceback, uuid
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parent
REPOSITORY, SOURCE = ROOT.parent, ROOT.parent / "src"
DATA, RESULTS = ROOT / "data", ROOT / "local-results"
JOBS_FILE, LOCK = DATA / "jobs.json", threading.RLock()

def field(key, label, kind="directory", default="", **extra):
    return {"key": key, "label": label, "type": kind, "default": str(default), **extra}

FEATURE_DIR = REPOSITORY / "processed_data" / "common_ca_distances"
AUXILIARY = REPOSITORY / "processed_data" / "auxiliary_inputs" / "analysis_replicas"
RESULTS_DIR = REPOSITORY / "results"
ANALYSES = {
    "distance-validation": {
        "name": "Distance-map ML validation", "group": "Machine learning",
        "summary": "PCA and pooled-frame versus held-out-replica classification.",
        "script": "analyze_nine_replica_distances.py",
        "fields": [field("feature-dir", "Feature directory", default=FEATURE_DIR), field("components", "Maximum PCA components", "integer", 50, min=2, max=200), field("seed", "Random seed", "integer", 2026, min=0)],
        "outputs": "PCA coordinates, loadings, holdout scores, predictions, confusion matrix, figures"},
    "state-discovery": {
        "name": "Hierarchical state discovery", "group": "Machine learning",
        "summary": "Replica-level PCA, clustering, state matching, recurrence, and dwell statistics.",
        "script": "hierarchical_state_discovery.py",
        "fields": [field("feature-dir", "Feature directory", default=FEATURE_DIR), field("clusters", "Clusters per replica", "integer", 2, min=2, max=12), field("components", "PCA components", "integer", 30, min=2, max=200), field("seed", "Random seed", "integer", 2026, min=0)],
        "outputs": "State assignments, populations, dwell statistics, matching table, figures"},
    "distance-determinants": {
        "name": "Reproducible distance determinants", "group": "Machine learning",
        "summary": "Ranks residue-pair distance effects that recur across independent replicas.",
        "script": "reproducible_distance_determinants.py",
        "fields": [field("feature-dir", "Feature directory", default=FEATURE_DIR), field("top", "Top pairs per comparison", "integer", 25, min=1, max=500)],
        "outputs": "Replica-matched effects, ranked determinant table, comparison figure"},
    "uncertainty": {
        "name": "Determinant uncertainty", "group": "Validation",
        "summary": "Block-bootstrap confidence intervals and replica-level effect uncertainty.",
        "script": "validate_determinant_uncertainty.py",
        "fields": [field("feature-dir", "Feature directory", default=FEATURE_DIR), field("determinants", "Determinant CSV", "file", RESULTS_DIR / "reproducible_distance_determinants" / "top_reproducible_distance_determinants.csv"), field("top-per-comparison", "Pairs per comparison", "integer", 25, min=1, max=500), field("block-ns", "Block length", "number", 10.0, min=0.01, unit="ns"), field("frame-interval-ns", "Frame interval", "number", 0.2, min=0.0001, unit="ns"), field("bootstrap-iterations", "Bootstrap iterations", "integer", 2000, min=100, max=100000), field("seed", "Random seed", "integer", 2026, min=0)],
        "outputs": "Confidence intervals, uncertainty forest plot, analysis summary"},
    "rmsd-validation": {
        "name": "Replica-aware RMSD validation", "group": "Validation",
        "summary": "PCA and held-out-replica classification from matched RMSD feature tables.",
        "script": "replica_rmsd_validation.py",
        "fields": [field("analysis-root", "RMSD input root", default=AUXILIARY), field("interval-ns", "Frame interval", "number", 0.2, min=0.0001, unit="ns"), field("duration-ns", "Matched duration", "number", 300.0, min=0.1, unit="ns"), field("seed", "Random seed", "integer", 2026, min=0)],
        "outputs": "Matched features, PCA coordinates/loadings, holdout scores, confusion matrix, figures"},
    "ligand-coupling": {
        "name": "Ligand-geometry coupling", "group": "Mechanistic analysis",
        "summary": "Frame-matched distance coupling to ligand COM, torsions, and hydration.",
        "script": "frame_matched_ligand_coupling.py",
        "fields": [field("analysis-root", "Auxiliary analysis root", default=AUXILIARY), field("feature-dir", "Feature directory", default=FEATURE_DIR), field("top-pairs-per-comparison", "Pairs per comparison", "integer", 25, min=1, max=500)],
        "outputs": "COM, torsion, and water correlations with replica-consensus summaries"},
    "prolif-coupling": {
        "name": "ProLIF contact coupling", "group": "Mechanistic analysis",
        "summary": "Couples structural determinants to frame-matched interaction fingerprints.",
        "script": "frame_matched_prolif_coupling.py",
        "fields": [field("analysis-root", "Auxiliary analysis root", default=AUXILIARY), field("feature-dir", "Feature directory", default=FEATURE_DIR), field("top-pairs-per-comparison", "Pairs per comparison", "integer", 25, min=1, max=500), field("minimum-occupancy", "Minimum occupancy", "number", 0.05, min=0, max=1), field("maximum-occupancy", "Maximum occupancy", "number", 0.95, min=0, max=1)],
        "outputs": "Distance-contact correlations, occupancy table, consensus summary, figures"},
    "evidence-integration": {
        "name": "Mechanistic evidence integration", "group": "Interpretation",
        "summary": "Combines determinant, interaction, hydration, ligand, and energetic evidence.",
        "script": "integrate_mechanistic_evidence.py",
        "fields": [field("analysis-root", "Auxiliary analysis root", default=AUXILIARY), field("distance-results", "Distance determinant results", default=RESULTS_DIR / "reproducible_distance_determinants"), field("pairs-per-comparison", "Pairs per comparison", "integer", 25, min=1, max=500)],
        "outputs": "Integrated pair/residue evidence, interpretation report, summary figure"},
    "structural-mapping": {
        "name": "Structural determinant mapping", "group": "Interpretation",
        "summary": "Maps ranked residue pairs to PDB, PyMOL, VMD, and network outputs.",
        "script": "map_determinants_to_structure.py",
        "fields": [field("determinants", "Determinant CSV", "file", RESULTS_DIR / "reproducible_distance_determinants" / "top_reproducible_distance_determinants.csv"), field("top-pairs", "Top pairs", "integer", 25, min=1, max=500)],
        "outputs": "Annotated PDB files, PyMOL/VMD scripts, mapping tables, network figures"},
    "common-landscape": {
        "name": "Common-landscape baseline", "group": "Feature extraction",
        "summary": "Extracts selected C-alpha distances, then performs PCA and clustering.",
        "script": "common_landscape_baseline.py",
        "fields": [field("project", "Input project root", default=REPOSITORY), field("stride", "Frame stride", "integer", 1, min=1), field("frame-interval-ns", "Frame interval", "number", 0.1, min=0.0001, unit="ns"), field("max-components", "Maximum PCA components", "integer", 50, min=2, max=200), field("clusters", "Clusters (blank = automatic)", "integer", "", optional=True, min=2, max=12), field("seed", "Random seed", "integer", 2026, min=0)],
        "outputs": "Distance matrices, PCA coordinates/loadings, cluster diagnostics, figures"},
}

def utcnow(): return datetime.now(timezone.utc).isoformat()

def read_jobs():
    with LOCK:
        try: return json.loads(JOBS_FILE.read_text()) if JOBS_FILE.exists() else []
        except (OSError, json.JSONDecodeError): return []

def write_jobs(jobs):
    with LOCK:
        DATA.mkdir(parents=True, exist_ok=True)
        temporary = JOBS_FILE.with_suffix(".tmp")
        temporary.write_text(json.dumps(jobs, indent=2) + "\n"); temporary.replace(JOBS_FILE)

def update_job(job_id, **changes):
    jobs = read_jobs()
    for job in jobs:
        if job["id"] == job_id:
            job.update(changes); write_jobs(jobs); return job
    raise KeyError(job_id)

def validate_value(spec, value):
    if value in (None, ""):
        if spec.get("optional"): return None
        raise ValueError(f"{spec['label']} is required.")
    kind = spec["type"]
    if kind in {"directory", "file"}:
        path = Path(str(value)).expanduser().resolve()
        if kind == "directory" and not path.is_dir(): raise ValueError(f"Directory not found: {path}")
        if kind == "file" and not path.is_file(): raise ValueError(f"File not found: {path}")
        return str(path)
    parsed = int(value) if kind == "integer" else float(value) if kind == "number" else str(value)
    if "min" in spec and parsed < spec["min"]: raise ValueError(f"{spec['label']} must be at least {spec['min']}.")
    if "max" in spec and parsed > spec["max"]: raise ValueError(f"{spec['label']} must not exceed {spec['max']}.")
    return parsed

def build_command(analysis_id, values, output_dir):
    analysis = ANALYSES.get(analysis_id)
    if not analysis: raise ValueError("Unknown analysis type.")
    command, normalized = [sys.executable, str(SOURCE / analysis["script"])], {}
    for spec in analysis["fields"]:
        value = validate_value(spec, values.get(spec["key"], spec.get("default")))
        if value is not None:
            normalized[spec["key"]] = value; command.extend([f"--{spec['key']}", str(value)])
    command.extend(["--output-dir", str(output_dir)])
    return analysis, command, normalized

def collect_outputs(directory):
    return [{"path": str(p.relative_to(directory)), "bytes": p.stat().st_size}
            for p in sorted(directory.rglob("*")) if p.is_file() and p.name != "run.log"][:500]

def run_job(job_id, command, output_dir):
    log_path = output_dir / "run.log"
    try:
        update_job(job_id, status="running", started_at=utcnow(), message="Analysis running")
        environment = os.environ.copy(); environment.setdefault("MPLBACKEND", "Agg")
        with log_path.open("w", encoding="utf-8") as log:
            process = subprocess.run(command, cwd=REPOSITORY, stdout=log, stderr=subprocess.STDOUT, env=environment, text=True)
        outputs, tail = collect_outputs(output_dir), log_path.read_text(errors="replace")[-3000:]
        update_job(job_id, status="complete" if process.returncode == 0 else "failed", finished_at=utcnow(), message=f"Created {len(outputs)} output files" if process.returncode == 0 else f"Exited with code {process.returncode}", log_tail=tail, outputs=outputs)
    except Exception as error:
        update_job(job_id, status="failed", finished_at=utcnow(), message=str(error), log_tail=traceback.format_exc()[-3000:], outputs=collect_outputs(output_dir))

class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs): super().__init__(*args, directory=str(ROOT), **kwargs)
    def json_response(self, payload, status=HTTPStatus.OK):
        encoded = json.dumps(payload).encode(); self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8"); self.send_header("Content-Length", str(len(encoded))); self.end_headers(); self.wfile.write(encoded)
    def read_json(self): return json.loads(self.rfile.read(int(self.headers.get("Content-Length", "0"))) or b"{}")
    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/api/status":
            deps = {name: bool(shutil.which(name)) for name in ("vmd", "namd2")}; deps["python"] = sys.version.split()[0]
            try: import numpy, sklearn, matplotlib; deps["ml_stack"] = True
            except ImportError: deps["ml_stack"] = False
            self.json_response({"online": True, "repository": str(REPOSITORY), "results": str(RESULTS), "dependencies": deps}); return
        if path == "/api/analyses": self.json_response({"analyses": [{"id": key, **value} for key, value in ANALYSES.items()]}); return
        if path == "/api/jobs": self.json_response({"jobs": read_jobs()}); return
        if path.startswith("/api/jobs/"):
            job = next((x for x in read_jobs() if x["id"] == path.rsplit("/", 1)[-1]), None)
            self.json_response(job or {"error": "Job not found"}, HTTPStatus.OK if job else HTTPStatus.NOT_FOUND); return
        if path.startswith("/download/"):
            parts = unquote(path).split("/", 3)
            if len(parts) != 4: self.send_error(HTTPStatus.NOT_FOUND); return
            base, relative = (RESULTS / parts[2]).resolve(), parts[3]; target = (base / relative).resolve()
            if base not in target.parents or not target.is_file(): self.send_error(HTTPStatus.NOT_FOUND); return
            self.send_response(HTTPStatus.OK); self.send_header("Content-Type", mimetypes.guess_type(target.name)[0] or "application/octet-stream"); self.send_header("Content-Length", str(target.stat().st_size)); self.send_header("Content-Disposition", f'attachment; filename="{target.name}"'); self.end_headers()
            with target.open("rb") as handle: shutil.copyfileobj(handle, self.wfile)
            return
        super().do_GET()
    def do_POST(self):
        if urlparse(self.path).path != "/api/run": self.json_response({"error": "Not found"}, HTTPStatus.NOT_FOUND); return
        try:
            payload = self.read_json(); job_id = datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
            output_dir = RESULTS / job_id; output_dir.mkdir(parents=True)
            analysis, command, values = build_command(payload.get("analysis"), payload.get("values") or {}, output_dir)
            job = {"id": job_id, "analysis": payload["analysis"], "name": analysis["name"], "status": "queued", "message": "Queued", "created_at": utcnow(), "output_dir": str(output_dir), "values": values, "command": command, "outputs": [], "log_tail": ""}
            jobs = read_jobs(); jobs.insert(0, job); write_jobs(jobs[:200]); threading.Thread(target=run_job, args=(job_id, command, output_dir), daemon=True).start(); self.json_response(job, HTTPStatus.ACCEPTED)
        except (ValueError, TypeError, json.JSONDecodeError) as error: self.json_response({"error": str(error)}, HTTPStatus.BAD_REQUEST)
        except Exception as error: self.json_response({"error": str(error)}, HTTPStatus.INTERNAL_SERVER_ERROR)
    def log_message(self, fmt, *args):
        if not args or not str(args[0]).startswith("GET /api/jobs"): super().log_message(fmt, *args)

def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--host", default="127.0.0.1"); parser.add_argument("--port", type=int, default=8765); args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True); DATA.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer((args.host, args.port), Handler); print(f"ReplicaLab running at http://{args.host}:{args.port}")
    try: server.serve_forever()
    except KeyboardInterrupt: pass

if __name__ == "__main__": main()
