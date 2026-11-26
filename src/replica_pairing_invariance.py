#!/usr/bin/env python3
"""Test determinant conclusions across every cross-system replica pairing."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import t


PROJECT = Path(__file__).resolve().parents[1]
PONE = PROJECT
DEFAULT_FEATURE_DIR = PROJECT / "processed_data/common_ca_distances"
COMPARISONS = (("Z944", "Apo"), ("mZ944", "Apo"), ("mZ944", "Z944"))
REPLICAS = (1, 2, 3)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, default=DEFAULT_FEATURE_DIR)
    parser.add_argument(
        "--original-results",
        type=Path,
        default=PROJECT
        / "results/reproducible_distance_determinants/"
        "top_reproducible_distance_determinants.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PONE / "revision_analysis/results/replica_pairing_invariance",
    )
    parser.add_argument("--top", type=int, default=25)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def welch_interval(first: np.ndarray, second: np.ndarray) -> tuple[float, float, float]:
    """Return a 95% Welch interval and degrees of freedom for two sets of run means."""
    n1, n2 = len(first), len(second)
    v1, v2 = float(first.var(ddof=1)), float(second.var(ddof=1))
    standard_error = np.sqrt(v1 / n1 + v2 / n2)
    numerator = (v1 / n1 + v2 / n2) ** 2
    denominator = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
    degrees_freedom = numerator / denominator if denominator > 0 else float("inf")
    difference = float(first.mean() - second.mean())
    half_width = float(t.ppf(0.975, degrees_freedom) * standard_error)
    return difference - half_width, difference + half_width, degrees_freedom


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    means: dict[tuple[str, int], np.ndarray] = {}
    sds: dict[tuple[str, int], np.ndarray] = {}
    feature_names: np.ndarray | None = None
    for system in ("Apo", "Z944", "mZ944"):
        for replica in REPLICAS:
            path = args.feature_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
            with np.load(path, allow_pickle=False) as data:
                names = data["feature_names"].astype(str)
                values = data["distances"].astype(float)
            if feature_names is None:
                feature_names = names
            elif not np.array_equal(feature_names, names):
                raise ValueError(f"Feature mismatch in {path}")
            means[(system, replica)] = values.mean(axis=0)
            sds[(system, replica)] = values.std(axis=0, ddof=1)
    assert feature_names is not None

    original: dict[tuple[str, str], dict[str, str]] = {}
    if args.original_results.exists():
        with args.original_results.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                original[(row["comparison"], row["feature"])] = row

    all_rows: list[dict[str, object]] = []
    assignment_rows: list[dict[str, object]] = []
    summaries: dict[str, dict[str, object]] = {}

    for numerator, denominator in COMPARISONS:
        label = f"{numerator}_vs_{denominator}"
        pair_differences = np.stack(
            [
                means[(numerator, i)] - means[(denominator, j)]
                for i in REPLICAS
                for j in REPLICAS
            ]
        ).reshape(3, 3, -1)
        pair_effects = np.stack(
            [
                np.divide(
                    means[(numerator, i)] - means[(denominator, j)],
                    np.sqrt((sds[(numerator, i)] ** 2 + sds[(denominator, j)] ** 2) / 2),
                    out=np.zeros_like(means[(numerator, i)]),
                    where=np.sqrt(
                        (sds[(numerator, i)] ** 2 + sds[(denominator, j)] ** 2) / 2
                    )
                    > 1e-12,
                )
                for i in REPLICAS
                for j in REPLICAS
            ]
        ).reshape(3, 3, -1)

        all_positive = np.all(pair_differences > 0, axis=(0, 1))
        all_negative = np.all(pair_differences < 0, axis=(0, 1))
        invariant = all_positive | all_negative
        minimum_effect = np.min(np.abs(pair_effects), axis=(0, 1))
        numerator_means = np.stack([means[(numerator, i)] for i in REPLICAS])
        denominator_means = np.stack([means[(denominator, j)] for j in REPLICAS])
        system_mean_difference = numerator_means.mean(axis=0) - denominator_means.mean(axis=0)

        comparison_rows: list[dict[str, object]] = []
        for k, feature in enumerate(feature_names):
            original_row = original.get((label, str(feature)))
            ci_low, ci_high, ci_df = welch_interval(
                numerator_means[:, k], denominator_means[:, k]
            )
            row: dict[str, object] = {
                "comparison": label,
                "feature": str(feature),
                "system_mean_difference_A": float(system_mean_difference[k]),
                "minimum_pair_difference_A": float(np.min(pair_differences[:, :, k])),
                "maximum_pair_difference_A": float(np.max(pair_differences[:, :, k])),
                "all_9_pairs_same_direction": bool(invariant[k]),
                "minimum_absolute_standardized_effect_all_9_pairs": float(minimum_effect[k]),
                "unpaired_replica_95_CI_low_A": ci_low,
                "unpaired_replica_95_CI_high_A": ci_high,
                "welch_degrees_of_freedom": ci_df,
                "unpaired_replica_CI_excludes_zero": bool(ci_low > 0 or ci_high < 0),
                "original_top25": original_row is not None,
                "original_rank": int(original_row["rank"]) if original_row else "",
            }
            for i in REPLICAS:
                for j in REPLICAS:
                    row[f"run{i}_vs_run{j}_difference_A"] = float(
                        pair_differences[i - 1, j - 1, k]
                    )
            comparison_rows.append(row)
            all_rows.append(row)

        robust = [row for row in comparison_rows if row["all_9_pairs_same_direction"]]
        robust.sort(
            key=lambda row: (
                float(row["minimum_absolute_standardized_effect_all_9_pairs"]),
                abs(float(row["system_mean_difference_A"])),
            ),
            reverse=True,
        )
        robust_top = {str(row["feature"]) for row in robust[: args.top]}
        write_csv(
            args.output_dir / f"{label}_top_pairing_invariant_features.csv",
            robust[: args.top],
        )

        original_top = [
            row for row in comparison_rows if bool(row["original_top25"])
        ]
        for permutation in itertools.permutations(REPLICAS):
            selected = np.stack(
                [pair_differences[i - 1, permutation[i - 1] - 1] for i in REPLICAS]
            )
            selected_consistent = np.all(selected > 0, axis=0) | np.all(
                selected < 0, axis=0
            )
            assignment_rows.append(
                {
                    "comparison": label,
                    "denominator_assignment_for_numerator_runs_1_2_3": "-".join(
                        map(str, permutation)
                    ),
                    "sign_consistent_features_all_2016": int(selected_consistent.sum()),
                    "original_top25_retained": sum(
                        bool(selected_consistent[np.where(feature_names == row["feature"])[0][0]])
                        for row in original_top
                    ),
                }
            )

        summaries[label] = {
            "features_tested": int(len(feature_names)),
            "features_invariant_across_all_9_pairs": int(invariant.sum()),
            "original_top25_invariant_count": sum(
                bool(row["all_9_pairs_same_direction"]) for row in original_top
            ),
            "original_top25_invariant_features": [
                str(row["feature"])
                for row in original_top
                if bool(row["all_9_pairs_same_direction"])
            ],
            "top25_under_pairing_invariant_ranking": sorted(robust_top),
            "pairing_invariant_features_with_unpaired_CI_excluding_zero": sum(
                bool(row["unpaired_replica_CI_excludes_zero"]) for row in robust
            ),
        }

    write_csv(args.output_dir / "all_pairing_invariance_results.csv", all_rows)
    write_csv(args.output_dir / "one_to_one_assignment_sensitivity.csv", assignment_rows)

    colors = {
        "Z944_vs_Apo": "#2F5DA8",
        "mZ944_vs_Apo": "#D1495B",
        "mZ944_vs_Z944": "#5B8E3E",
    }
    labels = {
        "Z944_vs_Apo": "Z944 vs apo",
        "mZ944_vs_Apo": "mZ944 vs apo",
        "mZ944_vs_Z944": "mZ944 vs Z944",
    }
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    assignments = [
        row["denominator_assignment_for_numerator_runs_1_2_3"]
        for row in assignment_rows
        if row["comparison"] == "Z944_vs_Apo"
    ]
    x = np.arange(len(assignments))
    for comparison in colors:
        selected = [row for row in assignment_rows if row["comparison"] == comparison]
        axes[0].plot(
            x,
            [row["sign_consistent_features_all_2016"] for row in selected],
            marker="o",
            linewidth=2,
            color=colors[comparison],
            label=labels[comparison],
        )
        axes[1].plot(
            x,
            [row["original_top25_retained"] for row in selected],
            marker="o",
            linewidth=2,
            color=colors[comparison],
            label=labels[comparison],
        )
    for axis in axes:
        axis.set_xticks(x, assignments, rotation=45, ha="right")
        axis.set_xlabel("Denominator-run assignment to numerator runs 1–3")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
    axes[0].set_ylabel("Sign-consistent distances (of 2,016)")
    axes[0].set_title("A  Sensitivity of all distances")
    axes[1].set_ylabel("Original top-25 distances retained")
    axes[1].set_ylim(0, 26)
    axes[1].set_title("B  Sensitivity of reported top features")
    axes[1].legend(frameon=False, fontsize=9)
    figure.savefig(args.output_dir / "replica_assignment_sensitivity.png", dpi=300)
    figure.savefig(args.output_dir / "replica_assignment_sensitivity.pdf")
    plt.close(figure)

    figure, axes = plt.subplots(1, 3, figsize=(15, 10), constrained_layout=True)
    for axis, comparison in zip(axes, colors):
        rows = sorted(
            [
                row
                for row in all_rows
                if row["comparison"] == comparison and row["original_top25"]
            ],
            key=lambda row: int(row["original_rank"]),
            reverse=True,
        )
        y = np.arange(len(rows))
        lows = np.asarray([row["minimum_pair_difference_A"] for row in rows], dtype=float)
        highs = np.asarray([row["maximum_pair_difference_A"] for row in rows], dtype=float)
        means_plot = np.asarray([row["system_mean_difference_A"] for row in rows], dtype=float)
        point_colors = [
            "#20854E" if row["all_9_pairs_same_direction"] else "#8A8A8A"
            for row in rows
        ]
        axis.hlines(y, lows, highs, color=point_colors, linewidth=2)
        axis.scatter(lows, y, color=point_colors, marker="|", s=45)
        axis.scatter(highs, y, color=point_colors, marker="|", s=45)
        axis.scatter(means_plot, y, color=point_colors, s=22, zorder=3)
        axis.axvline(0, color="black", linewidth=0.9)
        axis.set_yticks(
            y,
            [
                str(row["feature"]).replace("ca_dist_", "").replace("_A", "").replace("_", "–")
                for row in rows
            ],
            fontsize=7.5,
        )
        axis.set_xlabel("Distance difference across replica pairs (Å)")
        axis.set_title(labels[comparison])
        axis.grid(axis="x", color="#E5E5E5", linewidth=0.6)
    figure.suptitle(
        "All-pairs ranges for the original top-25 determinants\n"
        "Green: invariant direction across all nine pairs; gray: range crosses zero",
        fontsize=12,
    )
    figure.savefig(args.output_dir / "all_pair_effect_ranges_original_top25.png", dpi=300)
    figure.savefig(args.output_dir / "all_pair_effect_ranges_original_top25.pdf")
    plt.close(figure)

    summary = {
        "criterion": (
            "A feature is pairing-invariant only when all 3 x 3 = 9 cross-system "
            "replica mean differences have the same nonzero direction. This also "
            "guarantees directional consistency under all 3! = 6 one-to-one labelings."
        ),
        "ranking": (
            "Pairing-invariant features are ranked by the minimum absolute standardized "
            "effect across all nine replica pairs, then absolute system-mean difference."
        ),
        "comparisons": summaries,
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
