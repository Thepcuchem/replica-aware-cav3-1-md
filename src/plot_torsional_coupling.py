#!/usr/bin/env python3
"""Plot replica-recurrent protein-distance/ligand-torsion coupling."""

from __future__ import annotations

import csv
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results/frame_matched_ligand_coupling"
PHASE_MINIMUM = 0.75


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    output_dir = args.output_dir or args.results_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = read_csv(args.results_dir / "reproducible_distance_dihedral_coupling.csv")
    raw = read_csv(args.results_dir / "distance_dihedral_correlations_by_replica.csv")
    raw_index = {
        (row["system"], row["feature"], row["dihedral"], int(row["replica"])): row
        for row in raw
    }

    selected: list[dict[str, str]] = []
    for system in ("Z944", "mZ944"):
        for chi_index in range(1, 8):
            dihedral = f"chi{chi_index}_deg"
            eligible = [
                row for row in summary
                if row["system"] == system
                and row["dihedral"] == dihedral
                and float(row["phase_consistency"]) >= PHASE_MINIMUM
            ]
            if not eligible:
                raise RuntimeError(f"No phase-recurrent relationship: {system} {dihedral}")
            selected.append(
                max(eligible, key=lambda row: float(row["minimum_strength"]))
            )

    labels: list[str] = []
    strengths = np.zeros((len(selected), 3), dtype=float)
    phases = np.zeros((len(selected), 3), dtype=float)
    consistency = np.zeros(len(selected), dtype=float)
    for row_index, row in enumerate(selected):
        pair = row["feature"].replace("ca_dist_", "").replace("_A", "").replace("_", "–")
        chi = row["dihedral"].replace("_deg", "")
        labels.append(f"{row['system']} {chi}  |  {pair}")
        consistency[row_index] = float(row["phase_consistency"])
        for replica in (1, 2, 3):
            detail = raw_index[(row["system"], row["feature"], row["dihedral"], replica)]
            strengths[row_index, replica - 1] = float(detail["circular_linear_strength"])
            phases[row_index, replica - 1] = float(detail["phase_rad"])

    figure, (axis, p_axis) = plt.subplots(
        1, 2, figsize=(11.5, 9.2), gridspec_kw={"width_ratios": (4.8, 1.0)},
    )
    # Reserve a lower band for the marker-size key so it does not cover plotted
    # relationships in panel A after publication-size scaling.
    figure.subplots_adjust(left=0.27, right=0.95, bottom=0.18, top=0.88, wspace=0.30)
    y = np.arange(len(selected))
    for replica in range(3):
        scatter = axis.scatter(
            np.full(len(selected), replica), y,
            s=45 + 950 * strengths[:, replica] ** 2,
            c=phases[:, replica], cmap="twilight", vmin=-np.pi, vmax=np.pi,
            edgecolors="#263746", linewidths=0.45,
        )
    axis.set_xticks((0, 1, 2), ("run1", "run2", "run3"))
    axis.set_yticks(y, labels, fontsize=14)
    axis.set_xlim(-0.55, 2.55)
    axis.set_ylim(len(selected) - 0.4, -0.6)
    axis.grid(axis="x", color="#d8dde2", linewidth=0.8)
    axis.axhline(6.5, color="#607383", linewidth=1.2)
    axis.set_title("A   Circular-linear strength and phase by replica", loc="left", fontweight="bold")
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.045, pad=0.02)
    colorbar.set_label("Phase, φ (rad)")
    colorbar.set_ticks((-np.pi, 0, np.pi))
    colorbar.set_ticklabels(("−π", "0", "+π"))

    p_axis.barh(y, consistency, color="#4c829f", edgecolor="#263746", linewidth=0.5)
    p_axis.axvline(PHASE_MINIMUM, color="#b34845", linestyle="--", linewidth=1.2)
    p_axis.set_xlim(0, 1.02)
    p_axis.set_ylim(len(selected) - 0.4, -0.6)
    p_axis.set_yticks([])
    p_axis.set_xlabel("Phase consistency, P")
    p_axis.set_title("B   Across-replica\nphase recurrence", loc="left", fontweight="bold")
    for row_index, value in enumerate(consistency):
        p_axis.text(min(value + 0.025, 0.96), row_index, f"{value:.2f}", va="center", fontsize=14)

    for strength, label in ((0.10, "R = 0.10"), (0.25, "R = 0.25"), (0.40, "R = 0.40")):
        axis.scatter([], [], s=45 + 950 * strength**2, facecolor="#b7b7b7",
                     edgecolor="#263746", linewidth=0.45, label=label)
    axis.legend(
        title="Marker size", loc="upper center", bbox_to_anchor=(0.5, -0.085),
        ncol=3, frameon=False, fontsize=14, borderaxespad=0,
    )
    figure.suptitle(
        "Replica-recurrent coupling of ligand torsions to protein Cα–Cα distances",
        fontsize=18, fontweight="bold",
    )
    figure.savefig(output_dir / "reproducible_distance_dihedral_coupling.png", dpi=300)
    figure.savefig(output_dir / "reproducible_distance_dihedral_coupling.pdf")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
