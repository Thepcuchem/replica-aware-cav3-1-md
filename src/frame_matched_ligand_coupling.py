#!/usr/bin/env python3
"""Frame-match reproducible protein distances to ligand and water descriptors."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / ".deps"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr

ANALYSIS_ROOT = PROJECT / "processed_data" / "auxiliary_inputs" / "analysis_replicas"
FEATURE_DIR = PROJECT / "processed_data" / "common_ca_distances"
SYSTEMS = ("Z944", "mZ944")
CHI = tuple(f"chi{i}_deg" for i in range(1, 8))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, default=ANALYSIS_ROOT)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "frame_matched_ligand_coupling",
    )
    parser.add_argument("--top-pairs-per-comparison", type=int, default=25)
    return parser.parse_args()


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_distance_series(
    feature_dir: Path, system: str, replica: int, selected_features: set[str]
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    path = feature_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
    with np.load(path, allow_pickle=False) as data:
        names = data["feature_names"].astype(str)
        index = {name: position for position, name in enumerate(names)}
        missing = selected_features - set(index)
        if missing:
            raise ValueError(f"{path}: missing features {sorted(missing)}")
        relative_time = data["time_ns"].astype(float) - float(data["time_ns"][0])
        series = {
            feature: data["distances"][:, index[feature]].astype(float)
            for feature in selected_features
        }
    return relative_time, series


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ranked = read_csv(
        PROJECT
        / "results/reproducible_distance_determinants/"
        "top_reproducible_distance_determinants.csv"
    )
    chosen = [
        row
        for row in ranked
        if int(row["rank"]) <= args.top_pairs_per_comparison
    ]
    features = sorted({row["feature"] for row in chosen})

    com_rows = read_csv(
        args.analysis_root
        / "com_distance/Z944_vs_mZ944/com_distances_all_replicas_wide.csv"
    )
    dihedral_rows = read_csv(
        args.analysis_root
        / "ligand_dihedrals/Z944_vs_mZ944/dihedral_angles_all_replicas.csv"
    )

    distance_cache: dict[tuple[str, int], tuple[np.ndarray, dict[str, np.ndarray]]] = {}
    com_correlations: list[dict[str, object]] = []
    circular_correlations: list[dict[str, object]] = []
    for system in SYSTEMS:
        for replica in (1, 2, 3):
            time, series = load_distance_series(
                args.feature_dir, system, replica, set(features)
            )
            distance_cache[(system, replica)] = (time, series)
            com_selected = [
                row
                for row in com_rows
                if row["system"] == system and int(row["run"]) == replica
            ]
            dihedral_selected = [
                row
                for row in dihedral_rows
                if row["system"] == system and int(row["run"]) == replica
            ]
            com_time = np.asarray([float(row["time_ns"]) for row in com_selected])
            dih_time = np.asarray([float(row["time_ns"]) for row in dihedral_selected])
            for feature, values in series.items():
                at_com = np.interp(com_time, time, values)
                for metric in ("pocket_distance_A", "filter_distance_A"):
                    descriptor = np.asarray(
                        [float(row[metric]) for row in com_selected]
                    )
                    rho = spearmanr(at_com, descriptor).statistic
                    com_correlations.append(
                        {
                            "system": system,
                            "replica": replica,
                            "feature": feature,
                            "descriptor": metric,
                            "frames": len(descriptor),
                            "spearman_rho": float(rho),
                        }
                    )
                at_dih = np.interp(dih_time, time, values)
                for chi in CHI:
                    radians = np.deg2rad(
                        [float(row[chi]) for row in dihedral_selected]
                    )
                    r_sin = pearsonr(at_dih, np.sin(radians)).statistic
                    r_cos = pearsonr(at_dih, np.cos(radians)).statistic
                    circular_correlations.append(
                        {
                            "system": system,
                            "replica": replica,
                            "feature": feature,
                            "dihedral": chi,
                            "frames": len(radians),
                            "r_sin": float(r_sin),
                            "r_cos": float(r_cos),
                            "circular_linear_strength": float(
                                np.sqrt(r_sin * r_sin + r_cos * r_cos)
                            ),
                            "phase_rad": float(np.arctan2(r_sin, r_cos)),
                        }
                    )
    write_csv(args.output_dir / "distance_COM_correlations_by_replica.csv", com_correlations)
    write_csv(
        args.output_dir / "distance_dihedral_correlations_by_replica.csv",
        circular_correlations,
    )

    com_reproducible: list[dict[str, object]] = []
    for system in SYSTEMS:
        for feature in features:
            for descriptor in ("pocket_distance_A", "filter_distance_A"):
                rows = [
                    row
                    for row in com_correlations
                    if row["system"] == system
                    and row["feature"] == feature
                    and row["descriptor"] == descriptor
                ]
                values = np.asarray([float(row["spearman_rho"]) for row in rows])
                consistent = bool(np.all(values > 0) or np.all(values < 0))
                com_reproducible.append(
                    {
                        "system": system,
                        "feature": feature,
                        "descriptor": descriptor,
                        "run1_rho": values[0],
                        "run2_rho": values[1],
                        "run3_rho": values[2],
                        "sign_consistent": consistent,
                        "mean_rho": float(values.mean()),
                        "minimum_absolute_rho": float(np.min(np.abs(values))),
                    }
                )
    com_reproducible.sort(
        key=lambda row: (
            bool(row["sign_consistent"]),
            float(row["minimum_absolute_rho"]),
        ),
        reverse=True,
    )
    write_csv(args.output_dir / "reproducible_distance_COM_coupling.csv", com_reproducible)

    dihedral_reproducible: list[dict[str, object]] = []
    for system in SYSTEMS:
        for feature in features:
            for chi in CHI:
                rows = [
                    row
                    for row in circular_correlations
                    if row["system"] == system
                    and row["feature"] == feature
                    and row["dihedral"] == chi
                ]
                strengths = np.asarray(
                    [float(row["circular_linear_strength"]) for row in rows]
                )
                phases = np.asarray([float(row["phase_rad"]) for row in rows])
                phase_consistency = float(abs(np.mean(np.exp(1j * phases))))
                dihedral_reproducible.append(
                    {
                        "system": system,
                        "feature": feature,
                        "dihedral": chi,
                        "run1_strength": strengths[0],
                        "run2_strength": strengths[1],
                        "run3_strength": strengths[2],
                        "minimum_strength": float(strengths.min()),
                        "mean_strength": float(strengths.mean()),
                        "phase_consistency": phase_consistency,
                    }
                )
    dihedral_reproducible.sort(
        key=lambda row: (
            float(row["phase_consistency"]),
            float(row["minimum_strength"]),
        ),
        reverse=True,
    )
    write_csv(
        args.output_dir / "reproducible_distance_dihedral_coupling.csv",
        dihedral_reproducible,
    )

    water_correlations: list[dict[str, object]] = []
    for replica in (2, 3):
        water_path = (
            args.analysis_root
            / f"water_bridges/Z944/run{replica}/bridge_frame_timeseries.tsv"
        )
        water_data = read_csv(water_path, delimiter="\t")
        absolute_time = np.asarray([float(row["time_ns"]) for row in water_data])
        time, series = distance_cache[("Z944", replica)]
        # Checkpoint time is absolute for these runs; reconstruct from its first time.
        checkpoint_path = (
            args.feature_dir / f"z944_run{replica}_common_ca_distances.npz"
        )
        with np.load(checkpoint_path, allow_pickle=False) as data:
            checkpoint_absolute_time = data["time_ns"].astype(float)
        for feature, values in series.items():
            matched = np.interp(absolute_time, checkpoint_absolute_time, values)
            for descriptor in ("has_bridge", "event_count"):
                observed = np.asarray([float(row[descriptor]) for row in water_data])
                rho = spearmanr(matched, observed).statistic
                water_correlations.append(
                    {
                        "replica": replica,
                        "feature": feature,
                        "descriptor": descriptor,
                        "frames": len(observed),
                        "spearman_rho": float(rho),
                    }
                )
    write_csv(
        args.output_dir / "Z944_distance_water_correlations_by_replica.csv",
        water_correlations,
    )
    water_reproducible: list[dict[str, object]] = []
    for feature in features:
        for descriptor in ("has_bridge", "event_count"):
            rows = [
                row
                for row in water_correlations
                if row["feature"] == feature and row["descriptor"] == descriptor
            ]
            values = np.asarray([float(row["spearman_rho"]) for row in rows])
            water_reproducible.append(
                {
                    "feature": feature,
                    "descriptor": descriptor,
                    "run2_rho": values[0],
                    "run3_rho": values[1],
                    "sign_consistent": bool(
                        np.all(values > 0) or np.all(values < 0)
                    ),
                    "mean_rho": float(values.mean()),
                    "minimum_absolute_rho": float(np.min(np.abs(values))),
                }
            )
    water_reproducible.sort(
        key=lambda row: (
            bool(row["sign_consistent"]),
            float(row["minimum_absolute_rho"]),
        ),
        reverse=True,
    )
    write_csv(
        args.output_dir / "reproducible_Z944_distance_water_coupling.csv",
        water_reproducible,
    )

    figure, axes = plt.subplots(1, 2, figsize=(13, 7), constrained_layout=True)
    for axis, system in zip(axes, SYSTEMS):
        rows = [
            row
            for row in com_reproducible
            if row["system"] == system and row["sign_consistent"]
        ][:10]
        labels = [
            f"{row['feature'].replace('ca_dist_', '').replace('_A', '')}\n"
            f"{row['descriptor'].replace('_distance_A', '')}"
            for row in rows
        ]
        values = np.asarray(
            [[row["run1_rho"], row["run2_rho"], row["run3_rho"]] for row in rows],
            dtype=float,
        )
        image = axis.imshow(values, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
        axis.set_yticks(range(len(rows)), labels, fontsize=8)
        axis.set_xticks((0, 1, 2), ("run1", "run2", "run3"))
        axis.set_title(system)
    figure.colorbar(image, ax=axes, label="Spearman correlation")
    figure.suptitle("Replica-consistent protein-distance/ligand-COM coupling")
    figure.savefig(args.output_dir / "reproducible_distance_COM_coupling.png", dpi=300)
    figure.savefig(args.output_dir / "reproducible_distance_COM_coupling.pdf")
    plt.close(figure)

    summary = {
        "unique_protein_distances": len(features),
        "COM_frames_per_replica": 600,
        "dihedral_frames_per_replica": 600,
        "COM_sign_consistent_relationships": sum(
            bool(row["sign_consistent"]) for row in com_reproducible
        ),
        "strong_COM_relationships_min_abs_rho_0.2": sum(
            bool(row["sign_consistent"])
            and float(row["minimum_absolute_rho"]) >= 0.2
            for row in com_reproducible
        ),
        "water_bridge_runs": [2, 3],
        "water_sign_consistent_relationships": sum(
            bool(row["sign_consistent"]) for row in water_reproducible
        ),
        "caveat": (
            "Correlations are descriptive and temporally autocorrelated. COM and "
            "dihedral coupling require direction/phase consistency across all three "
            "replicas; water coupling currently covers Z944 runs 2 and 3 only."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
