#!/usr/bin/env python3
"""Quantify temporal and across-replica uncertainty for top distance determinants."""

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
from scipy.stats import t

FEATURE_DIR = PROJECT / "processed_data" / "common_ca_distances"
COMPARISONS = (("Z944", "Apo"), ("mZ944", "Apo"), ("mZ944", "Z944"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument(
        "--determinants",
        type=Path,
        default=PROJECT
        / "results/reproducible_distance_determinants/"
        "top_reproducible_distance_determinants.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "determinant_uncertainty_validation",
    )
    parser.add_argument("--top-per-comparison", type=int, default=25)
    parser.add_argument("--block-ns", type=float, default=10.0)
    parser.add_argument("--frame-interval-ns", type=float, default=0.2)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def block_means(values: np.ndarray, block_frames: int) -> np.ndarray:
    complete = len(values) // block_frames
    if complete < 2:
        raise ValueError("Fewer than two complete blocks")
    trimmed = values[: complete * block_frames]
    return trimmed.reshape(complete, block_frames, values.shape[1]).mean(axis=1)


def bootstrap_difference(
    first_blocks: np.ndarray,
    second_blocks: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    first_indices = rng.integers(
        0, len(first_blocks), size=(iterations, len(first_blocks))
    )
    second_indices = rng.integers(
        0, len(second_blocks), size=(iterations, len(second_blocks))
    )
    differences = (
        first_blocks[first_indices].mean(axis=1)
        - second_blocks[second_indices].mean(axis=1)
    )
    return np.quantile(differences, 0.025, axis=0), np.quantile(
        differences, 0.975, axis=0
    )


def excludes_zero(low: float, high: float) -> bool:
    return low > 0 or high < 0


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.determinants.open(encoding="utf-8") as handle:
        ranked = list(csv.DictReader(handle))
    chosen = [
        row for row in ranked if int(row["rank"]) <= args.top_per_comparison
    ]
    features = sorted({row["feature"] for row in chosen})
    matrices: dict[tuple[str, int], np.ndarray] = {}
    names = None
    for system in ("Apo", "Z944", "mZ944"):
        for replica in (1, 2, 3):
            path = args.feature_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
            with np.load(path, allow_pickle=False) as data:
                if names is None:
                    names = data["feature_names"].astype(str)
                matrices[(system, replica)] = data["distances"].astype(float)
    assert names is not None
    index = {name: position for position, name in enumerate(names)}
    selected_indices = np.asarray([index[feature] for feature in features])
    block_frames = int(round(args.block_ns / args.frame_interval_ns))
    rng = np.random.default_rng(args.seed)

    run_intervals: dict[tuple[str, int], tuple[np.ndarray, np.ndarray]] = {}
    run_differences: dict[str, np.ndarray] = {}
    for numerator, denominator in COMPARISONS:
        label = f"{numerator}_vs_{denominator}"
        differences = []
        for replica in (1, 2, 3):
            first = matrices[(numerator, replica)][:, selected_indices]
            second = matrices[(denominator, replica)][:, selected_indices]
            differences.append(first.mean(axis=0) - second.mean(axis=0))
            low, high = bootstrap_difference(
                block_means(first, block_frames),
                block_means(second, block_frames),
                args.bootstrap_iterations,
                rng,
            )
            run_intervals[(label, replica)] = (low, high)
        run_differences[label] = np.vstack(differences)

    output_rows: list[dict[str, object]] = []
    t_critical = float(t.ppf(0.975, df=2))
    for row in chosen:
        label = row["comparison"]
        feature_position = features.index(row["feature"])
        effects = run_differences[label][:, feature_position]
        mean = float(effects.mean())
        sd = float(effects.std(ddof=1))
        half_width = t_critical * sd / np.sqrt(3)
        across_low, across_high = mean - half_width, mean + half_width
        intervals = [
            (
                float(run_intervals[(label, replica)][0][feature_position]),
                float(run_intervals[(label, replica)][1][feature_position]),
            )
            for replica in (1, 2, 3)
        ]
        within_all = all(excludes_zero(low, high) for low, high in intervals)
        same_interval_direction = all(
            (low > 0 if mean > 0 else high < 0) for low, high in intervals
        )
        across_excludes = excludes_zero(across_low, across_high)
        output_rows.append(
            {
                "comparison": label,
                "original_rank": int(row["rank"]),
                "feature": row["feature"],
                "run1_mean_difference_A": effects[0],
                "run1_block_CI_low_A": intervals[0][0],
                "run1_block_CI_high_A": intervals[0][1],
                "run2_mean_difference_A": effects[1],
                "run2_block_CI_low_A": intervals[1][0],
                "run2_block_CI_high_A": intervals[1][1],
                "run3_mean_difference_A": effects[2],
                "run3_block_CI_low_A": intervals[2][0],
                "run3_block_CI_high_A": intervals[2][1],
                "replica_mean_difference_A": mean,
                "replica_SD_A": sd,
                "replica_t_CI_low_A": across_low,
                "replica_t_CI_high_A": across_high,
                "all_run_block_CIs_exclude_zero": within_all,
                "all_run_block_CIs_same_direction": same_interval_direction,
                "replica_t_CI_excludes_zero": across_excludes,
                "dual_uncertainty_supported": (
                    same_interval_direction and across_excludes
                ),
            }
        )
    output_rows.sort(
        key=lambda row: (
            bool(row["dual_uncertainty_supported"]),
            abs(float(row["replica_mean_difference_A"])),
        ),
        reverse=True,
    )
    write_csv(
        args.output_dir / "determinant_uncertainty_validation.csv", output_rows
    )

    figure, axes = plt.subplots(1, 3, figsize=(16, 8), constrained_layout=True)
    for axis, (numerator, denominator) in zip(axes, COMPARISONS):
        label = f"{numerator}_vs_{denominator}"
        rows = [
            row for row in output_rows if row["comparison"] == label
        ][:15]
        rows = list(reversed(rows))
        y = np.arange(len(rows))
        means = np.asarray([float(row["replica_mean_difference_A"]) for row in rows])
        lows = np.asarray([float(row["replica_t_CI_low_A"]) for row in rows])
        highs = np.asarray([float(row["replica_t_CI_high_A"]) for row in rows])
        colors = [
            "#2CA02C" if row["dual_uncertainty_supported"] else "#7F7F7F"
            for row in rows
        ]
        for position, mean_value, low, high, color in zip(
            y, means, lows, highs, colors
        ):
            axis.errorbar(
                mean_value,
                position,
                xerr=([mean_value - low], [high - mean_value]),
                fmt="none",
                ecolor=color,
                elinewidth=2,
                capsize=3,
            )
        axis.scatter(means, y, c=colors, s=35, zorder=3)
        axis.axvline(0, color="black", linewidth=0.8)
        axis.set_yticks(
            y,
            [
                row["feature"].replace("ca_dist_", "").replace("_A", "")
                for row in rows
            ],
            fontsize=8,
        )
        axis.set_xlabel("Mean distance difference (Å), 95% replica t interval")
        axis.set_title(label.replace("_", " "))
    figure.savefig(
        args.output_dir / "determinant_replica_uncertainty_forest.png", dpi=300
    )
    figure.savefig(
        args.output_dir / "determinant_replica_uncertainty_forest.pdf"
    )
    plt.close(figure)

    summary = {
        "top_features_per_comparison": args.top_per_comparison,
        "block_length_ns": args.block_ns,
        "complete_blocks_per_trajectory": 1501 // block_frames,
        "bootstrap_iterations": args.bootstrap_iterations,
        "independent_replicas": 3,
        "replica_t_degrees_of_freedom": 2,
        "dual_uncertainty_supported_counts": {
            f"{numerator}_vs_{denominator}": sum(
                row["comparison"] == f"{numerator}_vs_{denominator}"
                and bool(row["dual_uncertainty_supported"])
                for row in output_rows
            )
            for numerator, denominator in COMPARISONS
        },
        "interpretation": (
            "Dual support requires every run's 10-ns block-bootstrap interval to "
            "exclude zero in the common direction and the across-replica t interval "
            "to exclude zero."
        ),
        "caveat": (
            "Only three independent replicas are available; df=2 t intervals are "
            "wide and distribution-sensitive. Block bootstrap quantifies temporal "
            "uncertainty within runs, not biological-replica uncertainty."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
