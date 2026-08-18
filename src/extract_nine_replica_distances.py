#!/usr/bin/env python3
"""Checkpointed extraction of common C-alpha distances from nine MD replicas."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / ".deps"))

import MDAnalysis as mda
import numpy as np
from scipy.spatial.distance import pdist

ANALYSIS_ROOT = PROJECT / "raw_data"
DEFAULT_OUTPUT = PROJECT / "processed_data" / "common_ca_distances"
POCKET_RESIDS = (
    384, 387, 388, 391, 868, 872, 875, 876, 916, 917, 918, 920, 921,
    922, 948, 950, 951, 952, 953, 955, 956, 957, 959, 960, 1462, 1495,
    1498, 1499, 1502, 1505, 1506, 1509, 1510, 1816, 1819, 1820, 1823,
    1824,
)
FILTER_RESIDS = tuple(range(351, 358)) + tuple(range(919, 927)) + tuple(
    range(1459, 1467)
) + tuple(range(1776, 1783))
ANALYSIS_RESIDS = tuple(sorted(set(POCKET_RESIDS + FILTER_RESIDS)))

CONFIGS = {
    (system, replica): {
        "rmsd": f"maps/{system.lower()}_run{replica}_rmsd.dat",
        "psf": str(ANALYSIS_ROOT / system.lower() / f"run{replica}" / "system.psf"),
        "trajectory": str(ANALYSIS_ROOT / system.lower() / f"run{replica}"),
    }
    for system in ("Apo", "Z944", "mZ944")
    for replica in (1, 2, 3)
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--system", choices=("Apo", "Z944", "mZ944"))
    parser.add_argument("--replica", type=int, choices=(1, 2, 3))
    parser.add_argument("--duration-ns", type=float, default=300.0)
    parser.add_argument("--interval-ns", type=float, default=0.2)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_frame_map(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().lstrip("#").split()
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            values = line.split()
            parsed = dict(zip(header, values))
            rows.append(
                {
                    "time_ns": float(parsed["Time_ns"]),
                    "dcd": int(parsed["DCD"]),
                    "frame": int(parsed["Frame"]),
                }
            )
    return rows


def select_frames(
    rows: list[dict[str, float | int]], duration_ns: float, interval_ns: float
) -> list[dict[str, float | int]]:
    """Select unique source frames nearest a regular grid within the final window."""
    end = float(rows[-1]["time_ns"])
    start = end - duration_ns
    eligible = [row for row in rows if float(row["time_ns"]) >= start - interval_ns]
    times = np.asarray([float(row["time_ns"]) for row in eligible])
    targets = np.arange(start, end + interval_ns * 0.25, interval_ns)
    selected: list[dict[str, float | int]] = []
    previous_key: tuple[int, int] | None = None
    for target in targets:
        index = int(np.argmin(np.abs(times - target)))
        row = eligible[index]
        key = (int(row["dcd"]), int(row["frame"]))
        if key == previous_key:
            continue
        selected.append(row)
        previous_key = key
    if len(selected) < int(duration_ns / interval_ns * 0.95):
        raise ValueError(
            f"Only {len(selected)} unique frames selected; expected about "
            f"{duration_ns / interval_ns:.0f}"
        )
    return selected


def dcd_path(directory: Path, number: int) -> Path:
    return directory / f"md-{number:02d}.dcd"


def validate_checkpoint(path: Path, expected_system: str, expected_replica: int) -> bool:
    try:
        with np.load(path, allow_pickle=False) as data:
            return (
                str(data["system"].item()) == expected_system
                and int(data["replica"].item()) == expected_replica
                and data["distances"].ndim == 2
                and data["distances"].shape[1] == 2016
                and len(data["time_ns"]) == len(data["distances"])
            )
    except Exception:
        return False


def extract_one(
    system: str,
    replica: int,
    output_dir: Path,
    duration_ns: float,
    interval_ns: float,
    dry_run: bool,
    force: bool,
) -> dict[str, object]:
    config = CONFIGS[(system, replica)]
    psf = Path(str(config["psf"]))
    trajectory_dir = Path(str(config["trajectory"]))
    map_path = ANALYSIS_ROOT / str(config["rmsd"])
    output_path = output_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
    if output_path.exists() and not force:
        if validate_checkpoint(output_path, system, replica):
            print(f"SKIP verified checkpoint: {output_path}", flush=True)
            return {"system": system, "replica": replica, "status": "skipped"}
        raise ValueError(f"Invalid existing checkpoint (use --force): {output_path}")

    rows = read_frame_map(map_path)
    selected = select_frames(rows, duration_ns, interval_ns)
    dcd_numbers = sorted({int(row["dcd"]) for row in selected})
    missing = [
        str(dcd_path(trajectory_dir, number))
        for number in dcd_numbers
        if not dcd_path(trajectory_dir, number).is_file()
    ]
    if not psf.is_file() or missing:
        raise FileNotFoundError(
            f"{system} run{replica}: psf_exists={psf.is_file()}, missing={missing}"
        )
    plan = {
        "system": system,
        "replica": replica,
        "selected_frames": len(selected),
        "source_dcds": len(dcd_numbers),
        "start_ns": float(selected[0]["time_ns"]),
        "end_ns": float(selected[-1]["time_ns"]),
        "output": str(output_path),
    }
    print(json.dumps(plan), flush=True)
    if dry_run:
        return {**plan, "status": "dry-run"}

    output_dir.mkdir(parents=True, exist_ok=True)
    matrix = np.empty((len(selected), 2016), dtype=np.float32)
    times = np.empty(len(selected), dtype=np.float64)
    source_dcd = np.empty(len(selected), dtype=np.int16)
    source_frame = np.empty(len(selected), dtype=np.int32)
    row_lookup: dict[int, list[tuple[int, dict[str, float | int]]]] = {}
    for output_index, row in enumerate(selected):
        row_lookup.setdefault(int(row["dcd"]), []).append((output_index, row))

    universe: mda.Universe | None = None
    resids: np.ndarray | None = None
    started = time.monotonic()
    for dcd_index, number in enumerate(dcd_numbers, start=1):
        path = dcd_path(trajectory_dir, number)
        if universe is None:
            universe = mda.Universe(str(psf), str(path))
            atoms = universe.select_atoms(
                "protein and name CA and resid " + " ".join(map(str, ANALYSIS_RESIDS))
            )
            resids = atoms.resids.astype(np.int32)
            if len(resids) != len(ANALYSIS_RESIDS) or len(set(resids.tolist())) != len(resids):
                raise ValueError(
                    f"{system} run{replica}: expected {len(ANALYSIS_RESIDS)} unique "
                    f"C-alpha atoms; found {len(resids)}"
                )
        else:
            universe.load_new(str(path))
            atoms = universe.select_atoms(
                "protein and name CA and resid " + " ".join(map(str, ANALYSIS_RESIDS))
            )
        for output_index, row in row_lookup[number]:
            frame = int(row["frame"])
            if frame >= len(universe.trajectory):
                raise IndexError(f"{path}: requested frame {frame}, has {len(universe.trajectory)}")
            universe.trajectory[frame]
            matrix[output_index] = pdist(atoms.positions).astype(np.float32)
            times[output_index] = float(row["time_ns"])
            source_dcd[output_index] = number
            source_frame[output_index] = frame
        elapsed = time.monotonic() - started
        print(
            f"{system} run{replica}: DCD {dcd_index}/{len(dcd_numbers)} "
            f"({number:02d}), elapsed {elapsed / 60:.1f} min",
            flush=True,
        )

    assert resids is not None
    feature_names = np.asarray(
        [
            f"ca_dist_{resids[i]}_{resids[j]}_A"
            for i in range(len(resids))
            for j in range(i + 1, len(resids))
        ]
    )
    temporary = output_path.with_suffix(".npz.partial")
    with temporary.open("wb") as handle:
        np.savez_compressed(
            handle,
            distances=matrix,
            time_ns=times,
            source_dcd=source_dcd,
            source_frame=source_frame,
            resids=resids,
            feature_names=feature_names,
            system=np.asarray(system),
            replica=np.asarray(replica),
            duration_ns=np.asarray(duration_ns),
            target_interval_ns=np.asarray(interval_ns),
        )
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(output_path)
    if not validate_checkpoint(output_path, system, replica):
        raise RuntimeError(f"Checkpoint failed post-write validation: {output_path}")
    print(f"COMPLETE {output_path}", flush=True)
    return {**plan, "status": "complete"}


def main() -> int:
    args = parse_args()
    if (args.system is None) != (args.replica is None):
        raise ValueError("--system and --replica must be provided together")
    keys = [(args.system, args.replica)] if args.system else list(CONFIGS)
    results = [
        extract_one(
            str(system),
            int(replica),
            args.output_dir,
            args.duration_ns,
            args.interval_ns,
            args.dry_run,
            args.force,
        )
        for system, replica in keys
    ]
    if not args.dry_run:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "extraction_status.json").write_text(
            json.dumps(results, indent=2) + "\n", encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
