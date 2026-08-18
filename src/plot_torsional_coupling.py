#!/usr/bin/env python3
"""Plot replica-recurrent protein-distance/ligand-torsion coupling."""

from __future__ import annotations

import csv
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
    summary = read_csv(RESULTS / "reproducible_distance_dihedral_coupling.csv")
    raw = read_csv(RESULTS / "distance_dihedral_correlations_by_replica.csv")
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
        1, 2, figsize=(11.5, 8.2), gridspec_kw={"width_ratios": (4.8, 1.0)},
        constrained_layout=True,
    )
    y = np.arange(len(selected))
    for replica in range(3):
        scatter = axis.scatter(
            np.full(len(selected), replica), y,
            s=45 + 950 * strengths[:, replica] ** 2,
            c=phases[:, replica], cmap="twilight", vmin=-np.pi, vmax=np.pi,
            edgecolors="#263746", linewidths=0.45,
        )
    axis.set_xticks((0, 1, 2), ("run1", "run2", "run3"))
    axis.set_yticks(y, labels, fontsize=9)
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
        p_axis.text(min(value + 0.025, 0.96), row_index, f"{value:.2f}", va="center", fontsize=8)

    for strength, label in ((0.10, "R = 0.10"), (0.25, "R = 0.25"), (0.40, "R = 0.40")):
        axis.scatter([], [], s=45 + 950 * strength**2, facecolor="#b7b7b7",
                     edgecolor="#263746", linewidth=0.45, label=label)
    axis.legend(title="Marker size", loc="lower right", frameon=False, fontsize=8)
    figure.suptitle(
        "Replica-recurrent coupling of ligand torsions to protein Cα–Cα distances",
        fontsize=14, fontweight="bold",
    )
    figure.savefig(RESULTS / "reproducible_distance_dihedral_coupling.png", dpi=300)
    figure.savefig(RESULTS / "reproducible_distance_dihedral_coupling.pdf")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
