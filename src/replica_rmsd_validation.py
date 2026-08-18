#!/usr/bin/env python3
"""Replica-held-out validation using common RMSD descriptors from nine MD runs."""

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
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

FEATURES = ("PD_RMSD_A", "S6_RMSD_A", "SF_RMSD_A", "Pocket_RMSD_A", "Protein_RMSD_A")
SYSTEMS = ("Apo", "Z944", "mZ944")
COLORS = {"Apo": "#555555", "Z944": "#2537D8", "mZ944": "#E63737"}
FILES = {
    ("Z944", 1): "rmsd_run1_50_650ns.dat",
    ("Z944", 2): "rmsd_run2_0_500ns.dat",
    ("Z944", 3): "rmsd_run3_0_500ns.dat",
    ("mZ944", 1): "rmsd_mZ944_run1_0_600ns.dat",
    ("mZ944", 2): "rmsd_mZ944_run2_0_500ns.dat",
    ("mZ944", 3): "rmsd_mZ944_run3_0_500ns.dat",
    ("Apo", 1): "rmsd_apo_run1_0_600ns.dat",
    ("Apo", 2): "rmsd_apo_run2_0_500ns.dat",
    ("Apo", 3): "rmsd_apo_run3_0_500ns.dat",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=PROJECT / "processed_data" / "auxiliary_inputs" / "analysis_replicas",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--interval-ns", type=float, default=0.2)
    parser.add_argument("--duration-ns", type=float, default=300.0)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def read_rmsd(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().lstrip("#").split()
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            values = line.split()
            row = {name: float(value) for name, value in zip(header, values)}
            rows.append(row)
    return rows


def matched_window(
    rows: list[dict[str, float]], duration_ns: float, interval_ns: float
) -> list[dict[str, float]]:
    end = rows[-1]["Time_ns"]
    start = end - duration_ns
    count = int(round(duration_ns / interval_ns)) + 1
    targets = start + np.arange(count) * interval_ns
    times = np.asarray([row["Time_ns"] for row in rows])
    if np.any(np.diff(times) <= 0):
        raise ValueError("RMSD time values must be strictly increasing")
    selected: list[dict[str, float]] = []
    for index, target in enumerate(targets):
        selected.append(
            {
                "Time_ns": float(target),
                **{
                    feature: float(
                        np.interp(target, times, [row[feature] for row in rows])
                    )
                    for feature in FEATURES
                },
            }
        )
    return selected


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_pca(rows: list[dict[str, object]], output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for system in SYSTEMS:
        selected = [row for row in rows if row["system"] == system]
        axes[0].scatter(
            [row["PC1"] for row in selected],
            [row["PC2"] for row in selected],
            s=5, alpha=0.2, linewidths=0, color=COLORS[system], label=system,
        )
    axes[0].set_title("Matched RMSD frames")
    axes[0].legend(frameon=False, markerscale=3)

    markers = {1: "o", 2: "s", 3: "^"}
    for system in SYSTEMS:
        for run in (1, 2, 3):
            selected = [
                row for row in rows if row["system"] == system and row["replica"] == run
            ]
            axes[1].scatter(
                np.mean([row["PC1"] for row in selected]),
                np.mean([row["PC2"] for row in selected]),
                s=110, marker=markers[run], color=COLORS[system],
                edgecolor="white", linewidth=0.8,
                label=f"{system} run{run}",
            )
    axes[1].set_title("Replica centroids")
    axes[1].legend(frameon=False, fontsize=8, ncol=2)
    for axis in axes:
        axis.set_xlabel("PC1")
        axis.set_ylabel("PC2")
    figure.savefig(output, dpi=300)
    figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)


def main() -> int:
    args = parse_args()
    output = args.output_dir or PROJECT / "results" / "replica_rmsd_validation"
    output.mkdir(parents=True, exist_ok=True)

    raw = {key: read_rmsd(args.analysis_root / filename) for key, filename in FILES.items()}
    duration = min(
        args.duration_ns,
        min(rows[-1]["Time_ns"] - rows[0]["Time_ns"] for rows in raw.values()),
    )
    duration = np.floor(duration / args.interval_ns) * args.interval_ns

    records: list[dict[str, object]] = []
    for (system, replica), rows in raw.items():
        selected = matched_window(rows, duration, args.interval_ns)
        start = selected[0]["Time_ns"]
        for index, row in enumerate(selected):
            records.append(
                {
                    "system": system,
                    "replica": replica,
                    "frame_index": index,
                    "source_time_ns": row["Time_ns"],
                    "relative_time_ns": row["Time_ns"] - start,
                    **{feature: row[feature] for feature in FEATURES},
                }
            )
    write_csv(output / "matched_common_rmsd_features.csv", records)

    matrix = np.asarray([[float(row[name]) for name in FEATURES] for row in records])
    labels = np.asarray([str(row["system"]) for row in records])
    groups = np.asarray([int(row["replica"]) for row in records])
    scaled = StandardScaler().fit_transform(matrix)
    pca = PCA(n_components=5).fit(scaled)
    coordinates = pca.transform(scaled)
    latent_rows: list[dict[str, object]] = []
    for index, row in enumerate(records):
        latent_rows.append(
            {
                "system": row["system"],
                "replica": row["replica"],
                "source_time_ns": row["source_time_ns"],
                "relative_time_ns": row["relative_time_ns"],
                "PC1": float(coordinates[index, 0]),
                "PC2": float(coordinates[index, 1]),
            }
        )
    write_csv(output / "rmsd_pca_coordinates.csv", latent_rows)

    fold_rows: list[dict[str, object]] = []
    pooled_truth: list[str] = []
    pooled_prediction: list[str] = []
    for held_out in (1, 2, 3):
        train = groups != held_out
        test = groups == held_out
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                max_iter=3000, class_weight="balanced", random_state=args.seed
            ),
        )
        model.fit(matrix[train], labels[train])
        predicted = model.predict(matrix[test])
        pooled_truth.extend(labels[test])
        pooled_prediction.extend(predicted)
        fold_rows.append(
            {
                "held_out_replica": held_out,
                "test_frames": int(test.sum()),
                "accuracy": float(accuracy_score(labels[test], predicted)),
                "balanced_accuracy": float(
                    balanced_accuracy_score(labels[test], predicted)
                ),
            }
        )
    write_csv(output / "replica_holdout_scores.csv", fold_rows)
    confusion = confusion_matrix(pooled_truth, pooled_prediction, labels=list(SYSTEMS))
    confusion_rows = [
        {
            "true_system": true_system,
            **{f"predicted_{system}": int(confusion[i, j])
               for j, system in enumerate(SYSTEMS)},
        }
        for i, true_system in enumerate(SYSTEMS)
    ]
    write_csv(output / "pooled_confusion_matrix.csv", confusion_rows)

    loading_rows = []
    for component_index in range(2):
        for feature_index, feature in enumerate(FEATURES):
            loading_rows.append(
                {
                    "component": f"PC{component_index + 1}",
                    "feature": feature,
                    "loading": float(pca.components_[component_index, feature_index]),
                }
            )
    write_csv(output / "rmsd_pca_loadings.csv", loading_rows)
    plot_pca(latent_rows, output / "replica_rmsd_pca.png")

    summary = {
        "matched_duration_ns": float(duration),
        "sampling_interval_ns": args.interval_ns,
        "frames_per_replica": len(records) // 9,
        "total_frames": len(records),
        "features": list(FEATURES),
        "pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "mean_holdout_balanced_accuracy": float(
            np.mean([row["balanced_accuracy"] for row in fold_rows])
        ),
        "interpretation_limit": (
            "RMSDs are relative to run-specific references. This validates "
            "replica-level dynamic signatures, not absolute conformational states."
        ),
    }
    (output / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
