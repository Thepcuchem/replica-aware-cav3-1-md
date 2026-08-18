#!/usr/bin/env python3
"""Common-protein conformational-landscape baseline for apo, Z944, and mZ944.

The script intentionally uses only protein features shared by all three systems.
Ligand-contact, water-bridge, and ligand-dihedral descriptors belong in a
separate bound-system analysis and can later be used to explain these states.
"""

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
import MDAnalysis as mda
import numpy as np
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

SYSTEMS = ("apo", "z944", "mz944")
DISPLAY_NAMES = {"apo": "Apo", "z944": "Z944", "mz944": "mZ944"}
COLORS = {"apo": "#555555", "z944": "#2537D8", "mz944": "#E63737"}
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract common C-alpha distance maps and run PCA/clustering"
    )
    parser.add_argument("--project", type=Path, default=PROJECT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--frame-interval-ns", type=float, default=0.1)
    parser.add_argument("--max-components", type=int, default=50)
    parser.add_argument("--clusters", type=int, default=None)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def input_paths(project: Path, system: str) -> tuple[Path, Path]:
    folder = project / system
    return (
        folder / f"{system}_dl.psf",
        folder / f"{system}_final300ns_100ps.dcd",
    )


def extract_distances(
    project: Path, stride: int, frame_interval_ns: float
) -> tuple[np.ndarray, list[dict[str, object]], list[str]]:
    matrices: list[np.ndarray] = []
    metadata: list[dict[str, object]] = []
    reference_resids: list[int] | None = None

    for system in SYSTEMS:
        psf, dcd = input_paths(project, system)
        if not psf.is_file() or not dcd.is_file():
            raise FileNotFoundError(f"Missing prepared input for {system}: {psf} / {dcd}")
        universe = mda.Universe(str(psf), str(dcd))
        atoms = universe.select_atoms(
            "protein and name CA and resid " + " ".join(map(str, ANALYSIS_RESIDS))
        )
        resids = atoms.resids.astype(int).tolist()
        if len(resids) != len(ANALYSIS_RESIDS) or len(set(resids)) != len(resids):
            raise ValueError(
                f"{system}: expected {len(ANALYSIS_RESIDS)} unique C-alpha atoms, "
                f"found {len(resids)}"
            )
        if reference_resids is None:
            reference_resids = resids
        elif resids != reference_resids:
            raise ValueError(f"{system}: C-alpha residue order differs from other systems")

        frame_indices = range(0, len(universe.trajectory), stride)
        system_matrix = np.empty(
            (len(frame_indices), len(resids) * (len(resids) - 1) // 2),
            dtype=np.float32,
        )
        for row_index, frame_index in enumerate(frame_indices):
            universe.trajectory[frame_index]
            system_matrix[row_index] = pdist(atoms.positions).astype(np.float32)
            metadata.append(
                {
                    "system": DISPLAY_NAMES[system],
                    "frame": frame_index,
                    "time_ns": frame_index * frame_interval_ns,
                }
            )
        matrices.append(system_matrix)
        print(f"{DISPLAY_NAMES[system]}: extracted {len(system_matrix)} frames")

    assert reference_resids is not None
    feature_names = [
        f"ca_dist_{reference_resids[i]}_{reference_resids[j]}_A"
        for i in range(len(reference_resids))
        for j in range(i + 1, len(reference_resids))
    ]
    return np.vstack(matrices), metadata, feature_names


def select_cluster_count(
    coordinates: np.ndarray, seed: int, requested: int | None
) -> tuple[int, list[dict[str, float]]]:
    if requested is not None:
        if requested < 2:
            raise ValueError("--clusters must be at least 2")
        return requested, []
    scores: list[dict[str, float]] = []
    sample_size = min(5000, len(coordinates))
    for n_clusters in range(2, 9):
        model = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed)
        labels = model.fit_predict(coordinates)
        score = silhouette_score(
            coordinates, labels, sample_size=sample_size, random_state=seed
        )
        scores.append({"n_clusters": n_clusters, "silhouette_score": float(score)})
    best = max(scores, key=lambda row: row["silhouette_score"])
    return int(best["n_clusters"]), scores


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_latent(rows: list[dict[str, object]], output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for system in ("Apo", "Z944", "mZ944"):
        selected = [row for row in rows if row["system"] == system]
        key = system.lower()
        axes[0].scatter(
            [row["PC1"] for row in selected],
            [row["PC2"] for row in selected],
            s=7,
            alpha=0.42,
            linewidths=0,
            label=system,
            color=COLORS[key],
        )
    axes[0].set_title("Common protein landscape by system")
    axes[0].legend(frameon=False, markerscale=2)

    clusters = sorted({int(row["state"]) for row in rows})
    palette = plt.get_cmap("tab10")
    for cluster in clusters:
        selected = [row for row in rows if int(row["state"]) == cluster]
        axes[1].scatter(
            [row["PC1"] for row in selected],
            [row["PC2"] for row in selected],
            s=7,
            alpha=0.42,
            linewidths=0,
            label=f"State {cluster}",
            color=palette(cluster % 10),
        )
    axes[1].set_title("Unsupervised state assignments")
    axes[1].legend(frameon=False, markerscale=2)
    for axis in axes:
        axis.set_xlabel("PC1")
        axis.set_ylabel("PC2")
    figure.savefig(output, dpi=300)
    figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)


def main() -> int:
    args = parse_args()
    if args.stride < 1:
        raise ValueError("--stride must be at least 1")
    output = args.output_dir or args.project / "results" / "common_protein_baseline"
    output.mkdir(parents=True, exist_ok=True)

    distances, metadata, feature_names = extract_distances(
        args.project, args.stride, args.frame_interval_ns
    )
    np.savez_compressed(
        output / "common_ca_distance_features.npz",
        distances=distances,
        feature_names=np.asarray(feature_names),
    )

    scaler = StandardScaler()
    scaled = scaler.fit_transform(distances)
    n_components = min(args.max_components, scaled.shape[0] - 1, scaled.shape[1])
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=args.seed)
    coordinates = pca.fit_transform(scaled)
    cumulative = np.cumsum(pca.explained_variance_ratio_)
    cluster_dimensions = min(int(np.searchsorted(cumulative, 0.90) + 1), n_components)
    n_clusters, silhouette_rows = select_cluster_count(
        coordinates[:, :cluster_dimensions], args.seed, args.clusters
    )
    clusterer = KMeans(n_clusters=n_clusters, n_init=50, random_state=args.seed)
    labels = clusterer.fit_predict(coordinates[:, :cluster_dimensions])

    latent_rows: list[dict[str, object]] = []
    for index, meta in enumerate(metadata):
        latent_rows.append(
            {
                **meta,
                "PC1": float(coordinates[index, 0]),
                "PC2": float(coordinates[index, 1]),
                "state": int(labels[index]),
            }
        )
    write_csv(output / "latent_coordinates.csv", latent_rows)
    write_csv(output / "cluster_selection.csv", silhouette_rows)

    variance_rows = [
        {
            "component": index + 1,
            "explained_variance_ratio": float(value),
            "cumulative_explained_variance": float(cumulative[index]),
        }
        for index, value in enumerate(pca.explained_variance_ratio_)
    ]
    write_csv(output / "pca_explained_variance.csv", variance_rows)
    loading_rows: list[dict[str, object]] = []
    for component_index in range(min(2, pca.components_.shape[0])):
        component = pca.components_[component_index]
        ranked = np.argsort(np.abs(component))[::-1][:50]
        for rank, feature_index in enumerate(ranked, start=1):
            loading_rows.append(
                {
                    "component": f"PC{component_index + 1}",
                    "rank": rank,
                    "feature": feature_names[feature_index],
                    "loading": float(component[feature_index]),
                    "absolute_loading": float(abs(component[feature_index])),
                }
            )
    write_csv(output / "top_pca_distance_loadings.csv", loading_rows)

    population_rows: list[dict[str, object]] = []
    for system in ("Apo", "Z944", "mZ944"):
        system_labels = [int(row["state"]) for row in latent_rows if row["system"] == system]
        for state in range(n_clusters):
            count = system_labels.count(state)
            population_rows.append(
                {
                    "system": system,
                    "state": state,
                    "frames": count,
                    "population": count / len(system_labels),
                }
            )
    write_csv(output / "state_populations_by_system.csv", population_rows)
    plot_latent(latent_rows, output / "common_protein_landscape.png")

    summary = {
        "input_frames": {system: sum(row["system"] == system for row in metadata)
                         for system in ("Apo", "Z944", "mZ944")},
        "frame_interval_ns": args.frame_interval_ns,
        "stride": args.stride,
        "n_distance_features": len(feature_names),
        "pca_components_fitted": n_components,
        "components_used_for_clustering": cluster_dimensions,
        "variance_captured_for_clustering": float(cumulative[cluster_dimensions - 1]),
        "selected_clusters": n_clusters,
        "random_seed": args.seed,
        "scope_note": (
            "Exploratory one-trajectory-per-system baseline; independent replicas "
            "are required before kinetic or mechanistic claims."
        ),
    }
    (output / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote baseline results to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
