#!/usr/bin/env python3
"""Discover candidate states per replica and quantify within-system recurrence."""

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
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

SYSTEMS = ("Apo", "Z944", "mZ944")
COLORS = {"Apo": "#555555", "Z944": "#2537D8", "mZ944": "#E63737"}
FEATURE_DIR = PROJECT / "processed_data" / "common_ca_distances"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "hierarchical_state_discovery",
    )
    parser.add_argument("--clusters", type=int, default=2)
    parser.add_argument("--components", type=int, default=30)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_data(
    feature_dir: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    matrices, systems, replicas, times = [], [], [], []
    names = None
    for system in SYSTEMS:
        for replica in (1, 2, 3):
            path = feature_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
            with np.load(path, allow_pickle=False) as data:
                matrix = data["distances"].astype(np.float32, copy=False)
                if names is None:
                    names = data["feature_names"]
                elif not np.array_equal(names, data["feature_names"]):
                    raise ValueError(f"Feature mismatch in {path}")
                matrices.append(matrix)
                systems.extend([system] * len(matrix))
                replicas.extend([replica] * len(matrix))
                times.extend(data["time_ns"].tolist())
    assert names is not None
    return (
        np.vstack(matrices),
        np.asarray(systems),
        np.asarray(replicas),
        np.asarray(times, dtype=float),
        names,
    )


def dwell_statistics(labels: np.ndarray, interval_ns: float = 0.2) -> tuple[int, float, float]:
    lengths: list[int] = []
    current = 1
    for previous, value in zip(labels[:-1], labels[1:]):
        if value == previous:
            current += 1
        else:
            lengths.append(current)
            current = 1
    lengths.append(current)
    dwell = np.asarray(lengths) * interval_ns
    return len(lengths), float(np.median(dwell)), float(np.max(dwell))


def plot_system(
    system: str,
    coordinates: np.ndarray,
    system_mask: np.ndarray,
    replicas: np.ndarray,
    aligned_labels: np.ndarray,
    output: Path,
) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    markers = ("o", "s", "^", "D", "P", "X")
    for axis, replica in zip(axes, (1, 2, 3)):
        mask = system_mask & (replicas == replica)
        for state in sorted(set(aligned_labels[mask].tolist())):
            state_mask = mask & (aligned_labels == state)
            axis.scatter(
                coordinates[state_mask, 0],
                coordinates[state_mask, 1],
                s=7,
                alpha=0.38,
                linewidths=0,
                marker=markers[state % len(markers)],
                label=f"Matched state {state}",
            )
        axis.set_title(f"{system} run {replica}")
        axis.set_xlabel("Global PC1")
        axis.set_ylabel("Global PC2")
        axis.legend(frameon=False, fontsize=8)
    figure.savefig(output, dpi=300)
    figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)


def main() -> int:
    args = parse_args()
    if args.clusters < 2 or args.clusters > 6:
        raise ValueError("--clusters must be between 2 and 6")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix, systems, replicas, times, feature_names = load_data(args.feature_dir)
    scaled = StandardScaler().fit_transform(matrix)
    pca = PCA(
        n_components=args.components, svd_solver="randomized", random_state=args.seed
    )
    coordinates = pca.fit_transform(scaled)
    cluster_dimensions = min(
        int(np.searchsorted(np.cumsum(pca.explained_variance_ratio_), 0.80) + 1),
        args.components,
    )
    model_space = coordinates[:, :cluster_dimensions]

    raw_labels = np.full(len(matrix), -1, dtype=int)
    aligned_labels = np.full(len(matrix), -1, dtype=int)
    population_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    recurrence_rows: list[dict[str, object]] = []
    dwell_rows: list[dict[str, object]] = []

    for system in SYSTEMS:
        system_mask = systems == system
        centroids: dict[int, np.ndarray] = {}
        models: dict[int, KMeans] = {}
        for replica in (1, 2, 3):
            mask = system_mask & (replicas == replica)
            replica_space = model_space[mask]
            for candidate_k in range(2, 7):
                candidate = KMeans(
                    n_clusters=candidate_k, n_init=30, random_state=args.seed
                ).fit(replica_space)
                selection_rows.append(
                    {
                        "system": system,
                        "replica": replica,
                        "candidate_clusters": candidate_k,
                        "silhouette_score": float(
                            silhouette_score(
                                replica_space,
                                candidate.labels_,
                                sample_size=min(1500, len(replica_space)),
                                random_state=args.seed,
                            )
                        ),
                    }
                )
            model = KMeans(
                n_clusters=args.clusters, n_init=50, random_state=args.seed
            ).fit(replica_space)
            models[replica] = model
            centroids[replica] = model.cluster_centers_
            raw_labels[mask] = model.labels_

        reference = centroids[1]
        mappings: dict[int, dict[int, int]] = {
            1: {state: state for state in range(args.clusters)}
        }
        reference_pair_distances = np.linalg.norm(
            reference[:, None, :] - reference[None, :, :], axis=2
        )
        reference_separation = float(
            np.median(
                reference_pair_distances[
                    np.triu_indices(args.clusters, k=1)
                ]
            )
        )
        for replica in (2, 3):
            distances = np.linalg.norm(
                centroids[replica][:, None, :] - reference[None, :, :], axis=2
            )
            source, target = linear_sum_assignment(distances)
            mappings[replica] = {int(s): int(t) for s, t in zip(source, target)}
            for source_state, target_state in zip(source, target):
                recurrence_rows.append(
                    {
                        "system": system,
                        "reference_replica": 1,
                        "comparison_replica": replica,
                        "reference_state": int(target_state),
                        "comparison_state": int(source_state),
                        "centroid_distance": float(distances[source_state, target_state]),
                        "reference_within_replica_state_separation": reference_separation,
                        "distance_to_separation_ratio": float(
                            distances[source_state, target_state] / reference_separation
                        ),
                    }
                )

        for replica in (1, 2, 3):
            mask = system_mask & (replicas == replica)
            local = raw_labels[mask]
            aligned = np.asarray([mappings[replica][int(label)] for label in local])
            aligned_labels[mask] = aligned
            for state in range(args.clusters):
                state_mask = aligned == state
                # State-specific dwell runs, not dwell runs of the full label sequence.
                binary = state_mask.astype(int)
                state_lengths: list[int] = []
                length = 0
                for value in binary:
                    if value:
                        length += 1
                    elif length:
                        state_lengths.append(length)
                        length = 0
                if length:
                    state_lengths.append(length)
                dwell = np.asarray(state_lengths, dtype=float) * 0.2
                population_rows.append(
                    {
                        "system": system,
                        "replica": replica,
                        "matched_state": state,
                        "frames": int(state_mask.sum()),
                        "population": float(state_mask.mean()),
                    }
                )
                dwell_rows.append(
                    {
                        "system": system,
                        "replica": replica,
                        "matched_state": state,
                        "episodes": len(state_lengths),
                        "median_dwell_ns": float(np.median(dwell)) if len(dwell) else 0.0,
                        "max_dwell_ns": float(np.max(dwell)) if len(dwell) else 0.0,
                    }
                )
        plot_system(
            system,
            coordinates,
            system_mask,
            replicas,
            aligned_labels,
            args.output_dir / f"{system.lower()}_replica_candidate_states.png",
        )

    assignment_rows = [
        {
            "system": str(systems[index]),
            "replica": int(replicas[index]),
            "time_ns": float(times[index]),
            "PC1": float(coordinates[index, 0]),
            "PC2": float(coordinates[index, 1]),
            "local_state": int(raw_labels[index]),
            "matched_state": int(aligned_labels[index]),
        }
        for index in range(len(matrix))
    ]
    write_csv(args.output_dir / "candidate_state_assignments.csv", assignment_rows)
    write_csv(args.output_dir / "candidate_state_populations.csv", population_rows)
    write_csv(args.output_dir / "candidate_state_dwell_statistics.csv", dwell_rows)
    write_csv(args.output_dir / "cluster_number_diagnostics.csv", selection_rows)
    write_csv(args.output_dir / "state_centroid_matching.csv", recurrence_rows)

    summary = {
        "systems": list(SYSTEMS),
        "replicas_per_system": 3,
        "frames_per_replica": 1501,
        "candidate_states_per_replica": args.clusters,
        "global_pca_components": args.components,
        "components_used_for_clustering": cluster_dimensions,
        "variance_captured_for_clustering": float(
            np.sum(pca.explained_variance_ratio_[:cluster_dimensions])
        ),
        "matching_method": (
            "Within each system, run2/run3 state centroids were assigned one-to-one "
            "to run1 centroids by minimum Euclidean distance in global PCA space."
        ),
        "caveat": (
            "Matched labels indicate nearest structural centroids, not validated "
            "metastable equivalence. Kinetics require lag-time convergence tests."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
