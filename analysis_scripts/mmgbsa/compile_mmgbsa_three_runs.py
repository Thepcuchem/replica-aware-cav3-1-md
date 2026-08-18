#!/usr/bin/env python3
"""Compile run1-run3 MM/GBSA component means into comparison files."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
RUN1 = {
    "Z944": Path("/work1/ted/June2023/WithLigand/md_run1/bombam_analyse/mmgbsa"),
    "mZ944": Path("/work1/ted/June2023/WithmodLigand/md_run1/bombam_analyse/mmgbsa"),
}
COMPONENTS = ["Ebind_ELEC", "Ebind_VDW", "Ebind_TOTAL"]


def load_run1(system):
    root = RUN1[system]
    data = {}
    for component in COMPONENTS:
        suffix = component.removeprefix("Ebind_")
        data[component] = np.loadtxt(root / f"Ebind_{suffix}.dat")[:, 1]
    return data


def load_replica(system, run):
    frame = pd.read_csv(HERE / "mmgbsa" / system / f"run{run}" / "mmgbsa_binding_energy.csv")
    return {component: frame[component].to_numpy() for component in COMPONENTS}


records = []
for system in ("Z944", "mZ944"):
    for run in (1, 2, 3):
        data = load_run1(system) if run == 1 else load_replica(system, run)
        for component in COMPONENTS:
            values = data[component]
            records.append({
                "system": system,
                "run": run,
                "component": component,
                "frames": len(values),
                "mean_kcal_mol": values.mean(),
                "sd_kcal_mol": values.std(ddof=1),
            })

out = HERE / "mmgbsa" / "three_run_comparison"
out.mkdir(parents=True, exist_ok=True)
table = pd.DataFrame(records)
table.to_csv(out / "mmgbsa_run_statistics.csv", index=False)

replica = (
    table.groupby(["system", "component"])["mean_kcal_mol"]
    .agg(["mean", "std", "count"])
    .reset_index()
    .rename(columns={"mean": "mean_of_run_means", "std": "sd_across_run_means"})
)
replica.to_csv(out / "mmgbsa_three_run_statistics.csv", index=False)

with open(out / "summary.txt", "w") as handle:
    handle.write("Three-run MM/GBSA comparison\n")
    handle.write("============================\n\n")
    handle.write("Values are kcal/mol. Per-run uncertainty is the temporal SD.\n")
    handle.write("Three-run uncertainty is the SD among the three run means.\n")
    handle.write("Entropy was not calculated.\n\n")
    for system in ("Z944", "mZ944"):
        handle.write(f"{system}\n")
        for run in (1, 2, 3):
            handle.write(f"  run{run}: ")
            subset = table[(table.system == system) & (table.run == run)].set_index("component")
            handle.write(", ".join(
                f"{name}={subset.loc[name, 'mean_kcal_mol']:.3f} +/- "
                f"{subset.loc[name, 'sd_kcal_mol']:.3f}"
                for name in COMPONENTS
            ))
            handle.write("\n")
        handle.write("  Across-run means: ")
        subset = replica[replica.system == system].set_index("component")
        handle.write(", ".join(
            f"{name}={subset.loc[name, 'mean_of_run_means']:.3f} +/- "
            f"{subset.loc[name, 'sd_across_run_means']:.3f}"
            for name in COMPONENTS
        ))
        handle.write("\n\n")

fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
colors = {"Z944": "#2878B5", "mZ944": "#E07A1F"}
x = np.arange(3)
width = 0.36
for ax, component in zip(axes, COMPONENTS):
    for offset, system in zip((-width / 2, width / 2), ("Z944", "mZ944")):
        subset = table[(table.system == system) & (table.component == component)].sort_values("run")
        ax.bar(x + offset, subset.mean_kcal_mol, width, yerr=subset.sd_kcal_mol,
               capsize=3, color=colors[system], label=system)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x, ["run1", "run2", "run3"])
    ax.set_title(component.replace("Ebind_", ""))
    ax.grid(axis="y", alpha=0.2)
axes[0].set_ylabel("Binding-energy component (kcal/mol)")
axes[-1].legend(frameon=False)
fig.suptitle("NAMD/GBIS MM/GBSA across three MD replicas")
fig.tight_layout()
fig.savefig(out / "mmgbsa_three_run_comparison.png", dpi=300, bbox_inches="tight")
plt.close(fig)
