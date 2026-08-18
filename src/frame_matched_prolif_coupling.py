#!/usr/bin/env python3
"""Couple reproducible protein distances to residue-level ProLIF contacts."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / ".deps"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

ANALYSIS_ROOT = PROJECT / "processed_data" / "auxiliary_inputs" / "analysis_replicas"
FEATURE_DIR = PROJECT / "processed_data" / "common_ca_distances"
SYSTEMS = ("Z944", "mZ944")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, default=ANALYSIS_ROOT)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "frame_matched_prolif_coupling",
    )
    parser.add_argument("--top-pairs-per-comparison", type=int, default=25)
    parser.add_argument("--minimum-occupancy", type=float, default=0.05)
    parser.add_argument("--maximum-occupancy", type=float, default=0.95)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_fingerprint(path: Path) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 5:
        raise ValueError(f"Invalid fingerprint file: {path}")
    proteins = rows[1][1:]
    residue_columns: dict[int, list[int]] = {}
    for column, label in enumerate(proteins):
        match = re.search(r"([A-Z]{3})(\d+)", label)
        if not match:
            continue
        residue_columns.setdefault(int(match.group(2)), []).append(column)
    frames: list[int] = []
    contacts: dict[int, list[bool]] = {residue: [] for residue in residue_columns}
    for row in rows[4:]:
        if not row or not row[0].strip():
            continue
        frames.append(int(float(row[0])))
        values = [value.strip().lower() == "true" for value in row[1:]]
        for residue, columns in residue_columns.items():
            contacts[residue].append(any(values[column] for column in columns))
    return np.asarray(frames), {
        residue: np.asarray(values, dtype=bool) for residue, values in contacts.items()
    }


def fingerprint_path(root: Path, system: str, replica: int) -> Path:
    return root / f"prolif/{system}/run{replica}/interaction_fingerprint.csv"


def fingerprint_times(
    root: Path, system: str, replica: int, frames: np.ndarray
) -> np.ndarray:
    sampled = root / f"prolif/{system}/run{replica}/sampled_frames.csv"
    if sampled.is_file():
        with sampled.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        times = np.asarray([float(row["time_ns"]) for row in rows])
        if len(times) != len(frames):
            raise ValueError(f"Fingerprint/time length mismatch for {system} run{replica}")
        return times
    if system == "Z944" and replica == 1:
        # Native frame spacing is 0.01 ns. Align the final 300 ns of the
        # 318-ns fingerprint window to the final-300-ns distance checkpoint.
        raw = frames.astype(float) * 0.01
        return raw - (raw[-1] - 300.0)
    raise FileNotFoundError(f"No time mapping for {system} run{replica}")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ranked_path = (
        PROJECT
        / "results/reproducible_distance_determinants/"
        "top_reproducible_distance_determinants.csv"
    )
    with ranked_path.open(encoding="utf-8") as handle:
        ranked = list(csv.DictReader(handle))
    selected = [
        row for row in ranked if int(row["rank"]) <= args.top_pairs_per_comparison
    ]
    features = sorted({row["feature"] for row in selected})
    endpoints: dict[str, tuple[int, int]] = {}
    for feature in features:
        match = re.fullmatch(r"ca_dist_(\d+)_(\d+)_A", feature)
        if not match:
            raise ValueError(feature)
        endpoints[feature] = tuple(map(int, match.groups()))

    raw_rows: list[dict[str, object]] = []
    occupancy_rows: list[dict[str, object]] = []
    for system in SYSTEMS:
        for replica in (1, 2, 3):
            fp_frames, contacts = parse_fingerprint(
                fingerprint_path(args.analysis_root, system, replica)
            )
            fp_time = fingerprint_times(
                args.analysis_root, system, replica, fp_frames
            )
            checkpoint = (
                args.feature_dir
                / f"{system.lower()}_run{replica}_common_ca_distances.npz"
            )
            with np.load(checkpoint, allow_pickle=False) as data:
                names = data["feature_names"].astype(str)
                name_index = {name: index for index, name in enumerate(names)}
                checkpoint_time = data["time_ns"].astype(float)
                if replica == 1:
                    checkpoint_time = checkpoint_time - checkpoint_time[0]
                matrix = data["distances"]
                valid = (fp_time >= checkpoint_time[0]) & (
                    fp_time <= checkpoint_time[-1]
                )
                matched_time = fp_time[valid]
                for feature in features:
                    first, second = endpoints[feature]
                    distance = np.interp(
                        matched_time,
                        checkpoint_time,
                        matrix[:, name_index[feature]].astype(float),
                    )
                    for residue in (first, second):
                        if residue not in contacts:
                            continue
                        observed = contacts[residue][valid].astype(float)
                        occupancy = float(observed.mean())
                        occupancy_rows.append(
                            {
                                "system": system,
                                "replica": replica,
                                "residue": residue,
                                "frames": len(observed),
                                "occupancy": occupancy,
                            }
                        )
                        if not (
                            args.minimum_occupancy
                            <= occupancy
                            <= args.maximum_occupancy
                        ):
                            continue
                        rho = spearmanr(distance, observed).statistic
                        raw_rows.append(
                            {
                                "system": system,
                                "replica": replica,
                                "feature": feature,
                                "contact_residue": residue,
                                "frames": len(observed),
                                "contact_occupancy": occupancy,
                                "spearman_rho": float(rho),
                            }
                        )
    # Remove duplicated occupancy records caused by residues appearing in many pairs.
    occupancy_unique = {
        (row["system"], row["replica"], row["residue"]): row
        for row in occupancy_rows
    }
    write_csv(
        args.output_dir / "ProLIF_residue_occupancy_in_matched_window.csv",
        list(occupancy_unique.values()),
    )
    write_csv(args.output_dir / "distance_ProLIF_correlations_by_replica.csv", raw_rows)

    reproducible: list[dict[str, object]] = []
    keys = sorted(
        {
            (str(row["system"]), str(row["feature"]), int(row["contact_residue"]))
            for row in raw_rows
        }
    )
    for system, feature, residue in keys:
        rows = [
            row
            for row in raw_rows
            if row["system"] == system
            and row["feature"] == feature
            and int(row["contact_residue"]) == residue
        ]
        by_replica = {int(row["replica"]): row for row in rows}
        if set(by_replica) != {1, 2, 3}:
            continue
        values = np.asarray(
            [float(by_replica[replica]["spearman_rho"]) for replica in (1, 2, 3)]
        )
        reproducible.append(
            {
                "system": system,
                "feature": feature,
                "contact_residue": residue,
                "run1_rho": values[0],
                "run2_rho": values[1],
                "run3_rho": values[2],
                "sign_consistent": bool(np.all(values > 0) or np.all(values < 0)),
                "mean_rho": float(values.mean()),
                "minimum_absolute_rho": float(np.min(np.abs(values))),
            }
        )
    reproducible.sort(
        key=lambda row: (
            bool(row["sign_consistent"]),
            float(row["minimum_absolute_rho"]),
        ),
        reverse=True,
    )
    write_csv(
        args.output_dir / "reproducible_distance_ProLIF_coupling.csv",
        reproducible,
    )

    consistent = [row for row in reproducible if row["sign_consistent"]][:20]
    if consistent:
        values = np.asarray(
            [[row["run1_rho"], row["run2_rho"], row["run3_rho"]] for row in consistent],
            dtype=float,
        )
        labels = [
            f"{row['feature'].replace('ca_dist_', '').replace('_A', '')}"
            f" | contact {row['contact_residue']}"
            for row in consistent
        ]
        figure, axis = plt.subplots(
            figsize=(9, max(5, len(consistent) * 0.35)), constrained_layout=True
        )
        image = axis.imshow(values, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
        axis.set_yticks(range(len(labels)), labels, fontsize=8)
        axis.set_xticks((0, 1, 2), ("run1", "run2", "run3"))
        axis.set_title("Replica-consistent distance/ProLIF contact coupling")
        figure.colorbar(image, ax=axis, label="Point-biserial Spearman correlation")
        figure.savefig(
            args.output_dir / "reproducible_distance_ProLIF_coupling.png", dpi=300
        )
        figure.savefig(
            args.output_dir / "reproducible_distance_ProLIF_coupling.pdf"
        )
        plt.close(figure)

    summary = {
        "protein_distances_screened": len(features),
        "contact_definition": "Any ProLIF interaction type for the endpoint residue",
        "occupancy_filter": [
            args.minimum_occupancy,
            args.maximum_occupancy,
        ],
        "three_replica_testable_relationships": len(reproducible),
        "sign_consistent_relationships": sum(
            bool(row["sign_consistent"]) for row in reproducible
        ),
        "strong_relationships_min_abs_rho_0.2": sum(
            bool(row["sign_consistent"])
            and float(row["minimum_absolute_rho"]) >= 0.2
            for row in reproducible
        ),
        "caveat": (
            "ProLIF sampling is sparse (approximately 3-6 ns). Correlations are "
            "descriptive and limited to contacts with 5-95% occupancy in every run."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
