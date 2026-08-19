#!/usr/bin/env python3
"""Replica-held-out analysis of nine common C-alpha distance checkpoints."""

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
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SYSTEMS = ("Apo", "Z944", "mZ944")
COLORS = {"Apo": "#555555", "Z944": "#2537D8", "mZ944": "#E63737"}
DEFAULT_FEATURE_DIR = PROJECT / "processed_data" / "common_ca_distances"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, default=DEFAULT_FEATURE_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "nine_replica_distance_validation",
    )
    parser.add_argument("--components", type=int, default=50)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_checkpoints(
    feature_dir: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    matrices: list[np.ndarray] = []
    systems: list[str] = []
    replicas: list[int] = []
    times: list[float] = []
    reference_names: np.ndarray | None = None
    for system in SYSTEMS:
        for replica in (1, 2, 3):
            path = feature_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
            with np.load(path, allow_pickle=False) as data:
                matrix = data["distances"].astype(np.float32, copy=False)
                names = data["feature_names"]
                if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 1:
                    raise ValueError(f"{path}: expected a non-empty 2D distance matrix")
                if len(names) != matrix.shape[1]:
                    raise ValueError(
                        f"{path}: {matrix.shape[1]} columns but {len(names)} feature names"
                    )
                if reference_names is None:
                    reference_names = names
                elif not np.array_equal(reference_names, names):
                    raise ValueError(f"{path}: feature ordering differs")
                matrices.append(matrix)
                systems.extend([system] * len(matrix))
                replicas.extend([replica] * len(matrix))
                times.extend(data["time_ns"].astype(float).tolist())
    assert reference_names is not None
    return (
        np.vstack(matrices),
        np.asarray(systems),
        np.asarray(replicas),
        np.asarray(times),
        reference_names,
    )


def plot_landscape(rows: list[dict[str, object]], output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for system in SYSTEMS:
        selected = [row for row in rows if row["system"] == system]
        axes[0].scatter(
            [row["PC1"] for row in selected],
            [row["PC2"] for row in selected],
            s=5,
            alpha=0.2,
            linewidths=0,
            color=COLORS[system],
            label=system,
        )
    axes[0].set_title("All matched frames")
    axes[0].legend(frameon=False, markerscale=3)

    markers = {1: "o", 2: "s", 3: "^"}
    for system in SYSTEMS:
        for replica in (1, 2, 3):
            selected = [
                row
                for row in rows
                if row["system"] == system and row["replica"] == replica
            ]
            axes[1].scatter(
                np.mean([row["PC1"] for row in selected]),
                np.mean([row["PC2"] for row in selected]),
                s=110,
                marker=markers[replica],
                color=COLORS[system],
                edgecolor="white",
                linewidth=0.8,
                label=f"{system} run{replica}",
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
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix, labels, replicas, times, feature_names = load_checkpoints(args.feature_dir)
    print(f"Loaded {matrix.shape[0]} frames x {matrix.shape[1]} features", flush=True)

    scaler = StandardScaler()
    scaled = scaler.fit_transform(matrix)
    n_components = min(args.components, matrix.shape[0] - 1, matrix.shape[1])
    pca = PCA(
        n_components=n_components, svd_solver="randomized", random_state=args.seed
    )
    coordinates = pca.fit_transform(scaled)
    latent_rows = [
        {
            "system": str(labels[index]),
            "replica": int(replicas[index]),
            "time_ns": float(times[index]),
            "PC1": float(coordinates[index, 0]),
            "PC2": float(coordinates[index, 1]),
        }
        for index in range(len(matrix))
    ]
    write_csv(args.output_dir / "distance_pca_coordinates.csv", latent_rows)
    plot_landscape(latent_rows, args.output_dir / "nine_replica_distance_pca.png")

    variance_rows = [
        {
            "component": index + 1,
            "explained_variance_ratio": float(value),
            "cumulative_explained_variance": float(
                np.sum(pca.explained_variance_ratio_[: index + 1])
            ),
        }
        for index, value in enumerate(pca.explained_variance_ratio_)
    ]
    write_csv(args.output_dir / "pca_explained_variance.csv", variance_rows)
    loading_rows: list[dict[str, object]] = []
    for component_index in range(2):
        ranked = np.argsort(np.abs(pca.components_[component_index]))[::-1][:50]
        for rank, feature_index in enumerate(ranked, start=1):
            loading_rows.append(
                {
                    "component": f"PC{component_index + 1}",
                    "rank": rank,
                    "feature": str(feature_names[feature_index]),
                    "loading": float(pca.components_[component_index, feature_index]),
                }
            )
    write_csv(args.output_dir / "top_pca_distance_loadings.csv", loading_rows)

    fold_rows: list[dict[str, object]] = []
    trajectory_rows: list[dict[str, object]] = []
    pooled_truth: list[str] = []
    pooled_prediction: list[str] = []
    for held_out in (1, 2, 3):
        train = replicas != held_out
        test = replicas == held_out
        model = make_pipeline(
            StandardScaler(),
            PCA(
                n_components=n_components,
                svd_solver="randomized",
                random_state=args.seed,
            ),
            LogisticRegression(
                max_iter=3000, class_weight="balanced", random_state=args.seed
            ),
        )
        model.fit(matrix[train], labels[train])
        predicted = model.predict(matrix[test])
        probabilities = model.predict_proba(matrix[test])
        classes = model.classes_
        truth = labels[test]
        pooled_truth.extend(truth.tolist())
        pooled_prediction.extend(predicted.tolist())
        fold_rows.append(
            {
                "held_out_replica": held_out,
                "test_frames": int(test.sum()),
                "frame_accuracy": float(accuracy_score(truth, predicted)),
                "frame_balanced_accuracy": float(
                    balanced_accuracy_score(truth, predicted)
                ),
            }
        )
        test_labels = labels[test]
        for system in SYSTEMS:
            system_mask = test_labels == system
            mean_probability = probabilities[system_mask].mean(axis=0)
            trajectory_prediction = str(classes[int(np.argmax(mean_probability))])
            trajectory_rows.append(
                {
                    "replica": held_out,
                    "true_system": system,
                    "predicted_system": trajectory_prediction,
                    "correct": trajectory_prediction == system,
                    **{
                        f"mean_probability_{label}": float(
                            mean_probability[np.where(classes == label)[0][0]]
                        )
                        for label in SYSTEMS
                    },
                }
            )
        print(f"Completed held-out replica {held_out}", flush=True)

    write_csv(args.output_dir / "replica_holdout_scores.csv", fold_rows)
    write_csv(args.output_dir / "held_out_trajectory_predictions.csv", trajectory_rows)
    confusion = confusion_matrix(pooled_truth, pooled_prediction, labels=list(SYSTEMS))
    confusion_rows = [
        {
            "true_system": true_system,
            **{
                f"predicted_{system}": int(confusion[i, j])
                for j, system in enumerate(SYSTEMS)
            },
        }
        for i, true_system in enumerate(SYSTEMS)
    ]
    write_csv(args.output_dir / "pooled_frame_confusion_matrix.csv", confusion_rows)

    contact_rows: list[dict[str, object]] = []
    contact_trajectory_correct: list[bool] = []
    for held_out in (1, 2, 3):
        train = replicas != held_out
        test = replicas == held_out
        training_contacts = matrix[train] < 8.0
        occupancy = training_contacts.mean(axis=0)
        keep = (occupancy > 0.01) & (occupancy < 0.99)
        contact_model = SGDClassifier(
            loss="log_loss",
            penalty="l2",
            alpha=0.0001,
            max_iter=2000,
            tol=1e-4,
            class_weight="balanced",
            random_state=args.seed,
            average=True,
        )
        contact_model.fit(training_contacts[:, keep], labels[train])
        test_contacts = (matrix[test] < 8.0)[:, keep]
        contact_prediction = contact_model.predict(test_contacts)
        contact_probability = contact_model.predict_proba(test_contacts)
        test_labels = labels[test]
        contact_rows.append(
            {
                "held_out_replica": held_out,
                "selected_switching_contacts": int(keep.sum()),
                "frame_balanced_accuracy": float(
                    balanced_accuracy_score(test_labels, contact_prediction)
                ),
            }
        )
        for system in SYSTEMS:
            mean_probability = contact_probability[test_labels == system].mean(axis=0)
            guess = str(contact_model.classes_[int(np.argmax(mean_probability))])
            contact_trajectory_correct.append(guess == system)
    write_csv(args.output_dir / "contact_map_holdout_scores.csv", contact_rows)

    summary = {
        "feature_directory": str(args.feature_dir),
        "frames_per_replica": 1501,
        "total_frames": len(matrix),
        "distance_features": matrix.shape[1],
        "pca_components": n_components,
        "pca_explained_variance_ratio_first_five": (
            pca.explained_variance_ratio_[:5].tolist()
        ),
        "mean_frame_holdout_balanced_accuracy": float(
            np.mean([row["frame_balanced_accuracy"] for row in fold_rows])
        ),
        "trajectory_holdout_accuracy": float(
            np.mean([bool(row["correct"]) for row in trajectory_rows])
        ),
        "contact_map_cutoff_A": 8.0,
        "mean_contact_map_frame_holdout_balanced_accuracy": float(
            np.mean([row["frame_balanced_accuracy"] for row in contact_rows])
        ),
        "contact_map_trajectory_holdout_accuracy": float(
            np.mean(contact_trajectory_correct)
        ),
        "validation_design": (
            "Each fold holds out run N for all three systems and trains on the "
            "other six trajectories."
        ),
        "interpretation_limit": (
            "Frame observations are temporally correlated; trajectory-level "
            "predictions are the primary generalization check."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
