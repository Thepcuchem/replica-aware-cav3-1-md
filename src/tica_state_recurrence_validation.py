#!/usr/bin/env python3
"""Validate the PCA cluster-recurrence conclusion with time-lagged ICA."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
PONE = PROJECT
FEATURE_DIR = PROJECT / "processed_data/common_ca_distances"
SYSTEMS = ("Apo", "Z944", "mZ944")
REPLICAS = (1, 2, 3)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PONE / "revision_analysis/results/tica_state_recurrence_validation",
    )
    parser.add_argument("--lag-ns", type=float, nargs="+", default=(1.0, 5.0, 10.0))
    parser.add_argument("--frame-interval-ns", type=float, default=0.2)
    parser.add_argument("--pca-components", type=int, default=100)
    parser.add_argument("--tica-components", type=int, default=5)
    parser.add_argument("--clusters", type=int, default=2)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_trajectories(feature_dir: Path) -> tuple[list[np.ndarray], list[tuple[str, int]]]:
    trajectories: list[np.ndarray] = []
    identities: list[tuple[str, int]] = []
    names: np.ndarray | None = None
    for system in SYSTEMS:
        for replica in REPLICAS:
            path = feature_dir / f"{system.lower()}_run{replica}_common_ca_distances.npz"
            with np.load(path, allow_pickle=False) as data:
                current_names = data["feature_names"].astype(str)
                values = data["distances"].astype(np.float32)
            if names is None:
                names = current_names
            elif not np.array_equal(names, current_names):
                raise ValueError(f"Feature mismatch in {path}")
            trajectories.append(values)
            identities.append((system, replica))
    return trajectories, identities


def tica_coordinates(
    trajectories: list[np.ndarray], lag_frames: int, components: int
) -> tuple[list[np.ndarray], np.ndarray]:
    starts = np.vstack([trajectory[:-lag_frames] for trajectory in trajectories])
    ends = np.vstack([trajectory[lag_frames:] for trajectory in trajectories])
    kinetic_mean = (starts.sum(axis=0) + ends.sum(axis=0)) / (2 * len(starts))
    starts = starts - kinetic_mean
    ends = ends - kinetic_mean
    c00 = (starts.T @ starts + ends.T @ ends) / (2 * len(starts))
    c0t = (starts.T @ ends + ends.T @ starts) / (2 * len(starts))
    regularization = max(float(np.trace(c00)) / len(c00) * 1e-8, 1e-10)
    values, vectors = eigh(c0t, c00 + regularization * np.eye(len(c00)))
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    keep = min(components, len(values))
    return [
        (trajectory - kinetic_mean) @ vectors[:, :keep] for trajectory in trajectories
    ], values[:keep]


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw, identities = load_trajectories(args.feature_dir)

    lengths = [len(trajectory) for trajectory in raw]
    combined = np.vstack(raw)
    scaled = StandardScaler().fit_transform(combined)
    pca = PCA(
        n_components=args.pca_components,
        svd_solver="randomized",
        random_state=args.seed,
    )
    reduced = pca.fit_transform(scaled)
    split_at = np.cumsum(lengths)[:-1]
    reduced_trajectories = list(np.split(reduced, split_at))

    recurrence_rows: list[dict[str, object]] = []
    timescale_rows: list[dict[str, object]] = []
    assignment_rows: list[dict[str, object]] = []

    for lag_ns in args.lag_ns:
        lag_frames = int(round(lag_ns / args.frame_interval_ns))
        if lag_frames < 1 or lag_frames >= min(lengths):
            raise ValueError(f"Invalid lag {lag_ns} ns ({lag_frames} frames)")
        coordinates, eigenvalues = tica_coordinates(
            reduced_trajectories, lag_frames, args.tica_components
        )
        for component, eigenvalue in enumerate(eigenvalues, start=1):
            implied = (
                float(-lag_ns / np.log(abs(eigenvalue)))
                if 0 < abs(eigenvalue) < 1
                else float("nan")
            )
            timescale_rows.append(
                {
                    "lag_ns": lag_ns,
                    "component": component,
                    "tica_eigenvalue": float(eigenvalue),
                    "implied_timescale_ns": implied,
                }
            )

        for system in SYSTEMS:
            system_indices = [
                index for index, identity in enumerate(identities) if identity[0] == system
            ]
            centroids: dict[int, np.ndarray] = {}
            labels: dict[int, np.ndarray] = {}
            for index in system_indices:
                replica = identities[index][1]
                model = KMeans(
                    n_clusters=args.clusters, n_init=50, random_state=args.seed
                ).fit(coordinates[index])
                centroids[replica] = model.cluster_centers_
                labels[replica] = model.labels_
                for frame, (coordinate, label) in enumerate(
                    zip(coordinates[index], model.labels_)
                ):
                    assignment_rows.append(
                        {
                            "lag_ns": lag_ns,
                            "system": system,
                            "replica": replica,
                            "frame": frame,
                            "time_ns": frame * args.frame_interval_ns,
                            "local_state": int(label),
                            **{
                                f"tIC{component + 1}": float(value)
                                for component, value in enumerate(coordinate)
                            },
                        }
                    )

            reference = centroids[1]
            within = np.linalg.norm(reference[:, None] - reference[None, :], axis=2)
            reference_separation = float(
                np.median(within[np.triu_indices(args.clusters, k=1)])
            )
            for replica in (2, 3):
                distances = np.linalg.norm(
                    centroids[replica][:, None] - reference[None, :], axis=2
                )
                source, target = linear_sum_assignment(distances)
                for source_state, target_state in zip(source, target):
                    recurrence_rows.append(
                        {
                            "lag_ns": lag_ns,
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

    write_csv(args.output_dir / "tica_implied_timescales.csv", timescale_rows)
    write_csv(args.output_dir / "tica_state_centroid_recurrence.csv", recurrence_rows)
    write_csv(args.output_dir / "tica_state_assignments.csv", assignment_rows)

    colors = {"Apo": "#555555", "Z944": "#2F5DA8", "mZ944": "#D1495B"}
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    offsets = {"Apo": -0.12, "Z944": 0.0, "mZ944": 0.12}
    lag_positions = {lag: index for index, lag in enumerate(args.lag_ns)}
    for system in SYSTEMS:
        selected = [row for row in recurrence_rows if row["system"] == system]
        axes[0].scatter(
            [lag_positions[row["lag_ns"]] + offsets[system] for row in selected],
            [row["distance_to_separation_ratio"] for row in selected],
            color=colors[system],
            s=34,
            alpha=0.82,
            label=system,
        )
    axes[0].axhline(1, color="black", linestyle="--", linewidth=1.2)
    axes[0].set_xticks(range(len(args.lag_ns)), [str(lag) for lag in args.lag_ns])
    axes[0].set_xlabel("tICA lag time (ns)")
    axes[0].set_ylabel("Cross-replica centroid distance /\nwithin-run state separation")
    axes[0].set_title("A  Candidate-state recurrence")
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(axis="y", color="#DDDDDD", linewidth=0.7)

    component_colors = plt.cm.viridis(np.linspace(0.08, 0.9, args.tica_components))
    for component, color in zip(range(1, args.tica_components + 1), component_colors):
        selected = [
            row for row in timescale_rows if row["component"] == component
        ]
        axes[1].plot(
            [row["lag_ns"] for row in selected],
            [row["implied_timescale_ns"] for row in selected],
            marker="o",
            linewidth=1.8,
            color=color,
            label=f"tIC{component}",
        )
    axes[1].set_xlabel("tICA lag time (ns)")
    axes[1].set_ylabel("Implied timescale (ns)")
    axes[1].set_title("B  Lag-time sensitivity")
    axes[1].legend(frameon=False, fontsize=8, ncol=2)
    axes[1].grid(color="#DDDDDD", linewidth=0.7)
    figure.savefig(args.output_dir / "tica_recurrence_and_timescales.png", dpi=300)
    figure.savefig(args.output_dir / "tica_recurrence_and_timescales.pdf")
    plt.close(figure)

    ratios = np.asarray([row["distance_to_separation_ratio"] for row in recurrence_rows])
    by_lag = {
        str(lag): {
            "minimum_recurrence_ratio": float(
                min(
                    row["distance_to_separation_ratio"]
                    for row in recurrence_rows
                    if row["lag_ns"] == lag
                )
            ),
            "maximum_recurrence_ratio": float(
                max(
                    row["distance_to_separation_ratio"]
                    for row in recurrence_rows
                    if row["lag_ns"] == lag
                )
            ),
            "matched_centroids_with_ratio_below_1": sum(
                row["lag_ns"] == lag and row["distance_to_separation_ratio"] < 1
                for row in recurrence_rows
            ),
        }
        for lag in args.lag_ns
    }
    summary = {
        "method": (
            "Reversible tICA fitted to time-lagged pairs formed within each trajectory; "
            "trajectory boundaries were never connected. PCA to 100 dimensions was used "
            "only as covariance regularization before the generalized eigenproblem."
        ),
        "input_features": 2016,
        "pca_regularization_components": args.pca_components,
        "pca_regularization_variance_fraction": float(
            pca.explained_variance_ratio_.sum()
        ),
        "tica_components_clustered": args.tica_components,
        "candidate_states_per_replica": args.clusters,
        "lags_ns": list(args.lag_ns),
        "recurrence_by_lag": by_lag,
        "overall_minimum_recurrence_ratio": float(ratios.min()),
        "overall_maximum_recurrence_ratio": float(ratios.max()),
        "interpretation": (
            "A ratio below one is necessary for a matched cross-replica centroid to be "
            "closer than the two candidate states in reference run 1. Candidate clusters "
            "are not called metastable states without recurrence and lag-time convergence."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
