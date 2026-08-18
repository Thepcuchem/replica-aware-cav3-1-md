#!/usr/bin/env python3
"""Rank and plot per-residue MM/GBSA decomposition results."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--system", required=True)
    parser.add_argument("--run", required=True, type=int)
    args = parser.parse_args()
    root = args.directory.resolve()
    table = pd.read_csv(root / "per_residue" / "per_residue_statistics.csv")
    table["residue"] = table.resname.str.title() + table.resid.astype(str)
    ranked = table.sort_values("mean_TOTAL")
    ranked.to_csv(root / "per_residue" / "per_residue_ranked.csv", index=False)

    top = ranked.head(20).sort_values("mean_TOTAL", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(top.residue, top.mean_VDW, color="#E07A1F", label="van der Waals")
    ax.barh(top.residue, top.mean_ELEC, left=top.mean_VDW,
            color="#2878B5", label="Electrostatic-related")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Mean residue–ligand contribution (kcal/mol)")
    ax.set_ylabel("Residue")
    ax.set_title(f"{args.system} run{args.run}: top favorable per-residue contributions")
    ax.grid(axis="x", alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(root / "per_residue" / "top20_per_residue_decomposition.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(root / "per_residue" / "top20_summary.txt", "w") as handle:
        handle.write(f"{args.system} run{args.run}: top favorable residue contributions\n")
        handle.write("=" * 64 + "\n")
        handle.write("Values are mean +/- temporal SD in kcal/mol.\n\n")
        for _, row in ranked.head(20).iterrows():
            handle.write(
                f"{row['residue']:10s} TOTAL={row.mean_TOTAL:8.3f} +/- {row.sd_TOTAL:7.3f} "
                f"VDW={row.mean_VDW:8.3f} ELEC={row.mean_ELEC:8.3f}\n"
            )


if __name__ == "__main__":
    main()
