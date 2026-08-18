#!/usr/bin/env python3
"""Rank residue-pair distance effects that reproduce across three matched runs."""

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

SYSTEMS = ("Apo", "Z944", "mZ944")
COMPARISONS = (("Z944", "Apo"), ("mZ944", "Apo"), ("mZ944", "Z944"))
FEATURE_DIR = PROJECT / "processed_data" / "common_ca_distances"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "reproducible_distance_determinants",
    )
    parser.add_argument("--top", type=int, default=25)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrices: dict[tuple[str, int], np.ndarray] = {}
    feature_names = None
    for system in SYSTEMS:
        for replica in (1, 2, 3):
            path = args.feature_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
            with np.load(path, allow_pickle=False) as data:
                matrices[(system, replica)] = data["distances"].astype(float)
                if feature_names is None:
                    feature_names = data["feature_names"]
                elif not np.array_equal(feature_names, data["feature_names"]):
                    raise ValueError(f"Feature mismatch in {path}")
    assert feature_names is not None

    all_rows: list[dict[str, object]] = []
    top_by_comparison: dict[str, list[dict[str, object]]] = {}
    for numerator, denominator in COMPARISONS:
        differences = []
        effect_sizes = []
        for replica in (1, 2, 3):
            first = matrices[(numerator, replica)]
            second = matrices[(denominator, replica)]
            difference = first.mean(axis=0) - second.mean(axis=0)
            pooled_sd = np.sqrt((first.var(axis=0, ddof=1) + second.var(axis=0, ddof=1)) / 2)
            effect = np.divide(
                difference,
                pooled_sd,
                out=np.zeros_like(difference),
                where=pooled_sd > 1e-12,
            )
            differences.append(difference)
            effect_sizes.append(effect)
        differences_array = np.vstack(differences)
        effects_array = np.vstack(effect_sizes)
        sign_consistent = np.all(differences_array > 0, axis=0) | np.all(
            differences_array < 0, axis=0
        )
        label = f"{numerator}_vs_{denominator}"
        rows: list[dict[str, object]] = []
        for index, feature in enumerate(feature_names):
            mean_effect = float(effects_array[:, index].mean())
            row = {
                "comparison": label,
                "feature": str(feature),
                "run1_difference_A": float(differences_array[0, index]),
                "run2_difference_A": float(differences_array[1, index]),
                "run3_difference_A": float(differences_array[2, index]),
                "mean_difference_A": float(differences_array[:, index].mean()),
                "sd_difference_A": float(differences_array[:, index].std(ddof=1)),
                "sign_consistent": bool(sign_consistent[index]),
                "run1_standardized_effect": float(effects_array[0, index]),
                "run2_standardized_effect": float(effects_array[1, index]),
                "run3_standardized_effect": float(effects_array[2, index]),
                "mean_standardized_effect": mean_effect,
                "minimum_absolute_standardized_effect": float(
                    np.min(np.abs(effects_array[:, index]))
                ),
            }
            rows.append(row)
            all_rows.append(row)
        reproducible = [row for row in rows if row["sign_consistent"]]
        reproducible.sort(
            key=lambda row: (
                float(row["minimum_absolute_standardized_effect"]),
                abs(float(row["mean_standardized_effect"])),
            ),
            reverse=True,
        )
        top_by_comparison[label] = reproducible[: args.top]

    write_csv(args.output_dir / "all_replica_matched_distance_effects.csv", all_rows)
    top_rows = [
        {"rank": rank, **row}
        for label in top_by_comparison
        for rank, row in enumerate(top_by_comparison[label], start=1)
    ]
    write_csv(args.output_dir / "top_reproducible_distance_determinants.csv", top_rows)

    figure, axes = plt.subplots(1, 3, figsize=(16, 8), constrained_layout=True)
    for axis, (label, rows) in zip(axes, top_by_comparison.items()):
        shown = list(reversed(rows[:15]))
        values = np.asarray(
            [
                [
                    row["run1_standardized_effect"],
                    row["run2_standardized_effect"],
                    row["run3_standardized_effect"],
                ]
                for row in shown
            ],
            dtype=float,
        )
        limit = max(1.0, float(np.max(np.abs(values))))
        image = axis.imshow(
            values, aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit
        )
        axis.set_yticks(range(len(shown)))
        axis.set_yticklabels([str(row["feature"]).replace("ca_dist_", "").replace("_A", "")
                              for row in shown], fontsize=8)
        axis.set_xticks((0, 1, 2), ("run1", "run2", "run3"))
        axis.set_title(label.replace("_", " "))
        figure.colorbar(image, ax=axis, shrink=0.65, label="Standardized effect")
    figure.savefig(args.output_dir / "top_reproducible_distance_effects.png", dpi=300)
    figure.savefig(args.output_dir / "top_reproducible_distance_effects.pdf")
    plt.close(figure)

    counts = {
        label: sum(
            row["comparison"] == label and bool(row["sign_consistent"])
            for row in all_rows
        )
        for label in top_by_comparison
    }
    summary = {
        "replica_matching": "run1 vs run1, run2 vs run2, run3 vs run3",
        "features_tested_per_comparison": len(feature_names),
        "sign_consistent_feature_counts": counts,
        "ranking": (
            "Descending minimum absolute standardized effect across the three "
            "replica-matched comparisons, then absolute mean effect."
        ),
        "interpretation_limit": (
            "These are descriptive reproducibility rankings, not independent-frame "
            "significance tests. Temporal autocorrelation precludes treating frames "
            "as independent replicates."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
