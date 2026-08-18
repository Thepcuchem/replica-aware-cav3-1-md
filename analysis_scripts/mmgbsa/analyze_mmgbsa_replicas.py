#!/usr/bin/env python3
"""Extract NAMD energies, calculate binding terms, and summarize MM/GBSA."""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ENERGY_RE = re.compile(r"^ENERGY:\s+(\d+)\s+")


def read_namd(path):
    records = []
    startup_seen = False
    with open(path) as handle:
        for line in handle:
            if "Info: Finished startup at" in line:
                startup_seen = True
                continue
            if not startup_seen or not ENERGY_RE.match(line):
                continue
            fields = line.split()
            records.append((int(fields[1]), float(fields[6]), float(fields[7])))
    if not records:
        raise RuntimeError(f"No post-startup ENERGY records in {path}")
    # NAMD prints frame 0 before reading the first DCD coordinate. Keep only
    # positive frame indices, matching the run1 extraction scripts.
    records = [row for row in records if row[0] > 0]
    return pd.DataFrame(records, columns=["frame", "ELECT", "VDW"]).drop_duplicates("frame")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--system", required=True)
    parser.add_argument("--run", required=True, type=int)
    parser.add_argument("--start-ns", default=200.0, type=float)
    parser.add_argument("--snapshot-ns", default=0.5, type=float)
    args = parser.parse_args()
    root = args.directory.resolve()

    energies = {}
    for component in ("cpx", "prot", "lig"):
        energies[component] = read_namd(root / component / f"mmgbsa_{component}.out")

    frames = set(energies["cpx"].frame)
    for component in ("prot", "lig"):
        frames &= set(energies[component].frame)
    frames = sorted(frames)
    if not frames:
        raise RuntimeError("No common energy frames across complex, protein, and ligand")

    tables = {key: value.set_index("frame").loc[frames] for key, value in energies.items()}
    result = pd.DataFrame(index=frames)
    result.index.name = "frame"
    result["time_ns"] = args.start_ns + (np.arange(len(frames)) + 1) * args.snapshot_ns
    result["Ebind_ELEC"] = tables["cpx"].ELECT - tables["prot"].ELECT - tables["lig"].ELECT
    result["Ebind_VDW"] = tables["cpx"].VDW - tables["prot"].VDW - tables["lig"].VDW
    result["Ebind_TOTAL"] = result.Ebind_ELEC + result.Ebind_VDW
    result.reset_index().to_csv(root / "mmgbsa_binding_energy.csv", index=False)

    for column, filename in (
        ("Ebind_ELEC", "Ebind_ELEC.dat"),
        ("Ebind_VDW", "Ebind_VDW.dat"),
        ("Ebind_TOTAL", "Ebind_TOTAL.dat"),
    ):
        np.savetxt(root / filename, result.reset_index()[["frame", column]].to_numpy(),
                   fmt=["%d", "%.8f"])

    stats = result[["Ebind_ELEC", "Ebind_VDW", "Ebind_TOTAL"]].agg(["mean", "std", "min", "max"]).T
    stats.to_csv(root / "mmgbsa_statistics.csv")
    with open(root / "summary.txt", "w") as handle:
        handle.write(f"{args.system} run{args.run} NAMD/GBIS MM/GBSA summary\n")
        handle.write("=" * 52 + "\n\n")
        handle.write(f"Frames analyzed: {len(result)}\n")
        handle.write(f"Time interval: {result.time_ns.iloc[0]:.3f}-{result.time_ns.iloc[-1]:.3f} ns\n")
        handle.write(f"Snapshot interval: {args.snapshot_ns:.3f} ns\n")
        handle.write("Energy units: kcal/mol\n\n")
        for name, row in stats.iterrows():
            handle.write(
                f"{name:12s} mean={row['mean']:10.3f} SD={row['std']:9.3f} "
                f"min={row['min']:10.3f} max={row['max']:10.3f}\n"
            )
        handle.write("\nEbind_TOTAL = Ebind_ELEC + Ebind_VDW\n")
        handle.write("Entropy was not calculated.\n")

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    colors = {"Ebind_ELEC": "#4472C4", "Ebind_VDW": "#ED7D31", "Ebind_TOTAL": "#2E8B57"}
    labels = {"Ebind_ELEC": "Electrostatic-related", "Ebind_VDW": "van der Waals", "Ebind_TOTAL": "Total"}
    for ax, column in zip(axes, ("Ebind_ELEC", "Ebind_VDW", "Ebind_TOTAL")):
        ax.plot(result.time_ns, result[column], lw=0.8, color=colors[column])
        ax.axhline(result[column].mean(), ls="--", lw=1, color="black")
        ax.set_ylabel(f"{labels[column]}\n(kcal/mol)")
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("Time (ns)")
    fig.suptitle(f"{args.system} run{args.run} MM/GBSA binding-energy components")
    fig.tight_layout()
    fig.savefig(root / "mmgbsa_binding_energy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
