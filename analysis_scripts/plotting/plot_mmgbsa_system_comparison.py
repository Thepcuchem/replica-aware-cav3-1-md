#!/usr/bin/env python3
"""Create publication-ready Z944 versus mZ944 MM/GBSA comparison plots."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
MMGBSA = ROOT / "mmgbsa"
OUT = MMGBSA / "Z944_vs_mZ944_comparison"
RUN1 = {
    "Z944": Path("/work1/ted/June2023/WithLigand/md_run1/bombam_analyse/"
                 "mmgbsa_perres/FINAL_RESULT_DECOMPOSITION_ENERGY.dat"),
    "mZ944": Path("/work1/ted/June2023/WithmodLigand/md_run1/bombam_analyse/"
                  "mmgbsa_perres/FINAL_RESULT_DECOMPOSITION_ENERGY.dat"),
}
SYSTEMS = ["Z944", "mZ944"]
COMPONENTS = ["VDW", "ELEC", "TOTAL"]
COLORS = {"Z944": "#2878B5", "mZ944": "#D9534F"}


def residue_name_map():
    maps = []
    for system in SYSTEMS:
        for run in (2, 3):
            path = MMGBSA / system / f"run{run}" / "per_residue" / "per_residue_statistics.csv"
            frame = pd.read_csv(path, dtype={"resid": str})
            maps.append(frame[["resid", "resname"]])
    names = pd.concat(maps).drop_duplicates("resid").set_index("resid")["resname"]
    return names.to_dict()


def load_per_residue():
    names = residue_name_map()
    tables = []
    for system in SYSTEMS:
        old = pd.read_csv(RUN1[system], dtype={"Residue": str})
        old = old.rename(columns={
            "Residue": "resid", "AvgVDW": "VDW", "AvgELEC": "ELEC",
            "AvgTOTAL": "TOTAL", "FramesUsed": "frames",
        })
        old["system"] = system
        old["run"] = 1
        old["resname"] = old["resid"].map(names).fillna("UNK")
        tables.append(old[["system", "run", "resid", "resname", "frames"] + COMPONENTS])

        for run in (2, 3):
            path = MMGBSA / system / f"run{run}" / "per_residue" / "per_residue_statistics.csv"
            new = pd.read_csv(path, dtype={"resid": str}).rename(columns={
                "mean_VDW": "VDW", "mean_ELEC": "ELEC", "mean_TOTAL": "TOTAL",
            })
            new["system"] = system
            new["run"] = run
            tables.append(new[["system", "run", "resid", "resname", "frames"] + COMPONENTS])
    return pd.concat(tables, ignore_index=True)


def save_per_residue_tables(data):
    long_path = OUT / "per_residue_all_runs.csv"
    data.to_csv(long_path, index=False)

    summary = (
        data.groupby(["system", "resid", "resname"])[COMPONENTS]
        .agg(["mean", "std"]).reset_index()
    )
    summary.columns = [
        "_".join(str(x) for x in col if x).rstrip("_") if isinstance(col, tuple) else col
        for col in summary.columns
    ]
    summary.to_csv(OUT / "per_residue_three_run_summary.csv", index=False)
    return summary


def choose_residues(summary, number_per_system=12):
    selected = set()
    for system in SYSTEMS:
        part = summary[summary["system"] == system].nsmallest(number_per_system, "TOTAL_mean")
        selected.update(part["resid"])
    numeric = lambda x: int(x) if str(x).lstrip("-").isdigit() else 10**9
    return sorted(selected, key=numeric)


def plot_per_residue(summary, residues):
    names = (
        summary.drop_duplicates("resid").set_index("resid")["resname"].to_dict()
    )
    labels = [f"{names.get(r, 'UNK').title()}{r}" for r in residues]
    x = np.arange(len(residues))
    width = 0.37

    fig, axes = plt.subplots(3, 1, figsize=(14.5, 11.5), sharex=True)
    for ax, component in zip(axes, COMPONENTS):
        for offset, system in zip((-width / 2, width / 2), SYSTEMS):
            part = summary[summary["system"] == system].set_index("resid").reindex(residues)
            means = part[f"{component}_mean"].to_numpy(float)
            sds = part[f"{component}_std"].to_numpy(float)
            ax.bar(
                x + offset, means, width, yerr=sds, capsize=2.5,
                color=COLORS[system], edgecolor="black", linewidth=0.45,
                label=system, error_kw={"elinewidth": 0.8},
            )
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_ylabel(f"$\\Delta E_{{\\mathrm{{{component}}}}}$\n(kcal/mol)")
        ax.grid(axis="y", alpha=0.22, linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, ncol=2, loc="lower right")
    axes[0].set_title("Per-residue MM/GBSA decomposition: Z944 versus mZ944")
    axes[-1].set_xticks(x, labels, rotation=48, ha="right")
    axes[-1].set_xlabel("Protein residue")
    fig.text(
        0.5, 0.005,
        "Bars: mean of run1–run3; error bars: SD among replica means. "
        "Residues are the union of the 12 most favorable total-energy contributors in either system.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"per_residue_energy_comparison.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def load_overall():
    path = MMGBSA / "three_run_comparison" / "mmgbsa_run_statistics.csv"
    data = pd.read_csv(path)
    data["component"] = data["component"].str.replace("Ebind_", "", regex=False)
    return data


def plot_overall(data):
    rows = []
    for system in SYSTEMS:
        for component in COMPONENTS:
            values = data[(data.system == system) & (data.component == component)]["mean_kcal_mol"]
            rows.append({
                "system": system, "run": "Mean", "component": component,
                "mean_kcal_mol": values.mean(), "sd_kcal_mol": values.std(ddof=1),
                "error_definition": "SD among three replica means",
            })
    aggregate = pd.DataFrame(rows)
    individual = data.copy()
    individual["run"] = individual["run"].map(lambda n: f"Run {int(n)}")
    individual["error_definition"] = "Temporal SD within replica"
    complete = pd.concat([individual, aggregate], ignore_index=True)
    complete.to_csv(OUT / "overall_energy_plot_data.csv", index=False)

    categories = ["Run 1", "Run 2", "Run 3", "Mean"]
    x = np.arange(len(categories))
    width = 0.36
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
    for ax, component in zip(axes, COMPONENTS):
        for offset, system in zip((-width / 2, width / 2), SYSTEMS):
            part = complete[
                (complete.system == system) & (complete.component == component)
            ].set_index("run").reindex(categories)
            ax.bar(
                x + offset, part["mean_kcal_mol"], width,
                yerr=part["sd_kcal_mol"], capsize=3,
                color=COLORS[system], edgecolor="black", linewidth=0.55,
                label=system, error_kw={"elinewidth": 0.9},
            )
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xticks(x, categories, rotation=25, ha="right")
        ax.set_title(component)
        ax.grid(axis="y", alpha=0.22, linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Binding-energy contribution (kcal/mol)")
    axes[0].legend(frameon=False)
    fig.suptitle("Overall MM/GBSA binding energy: Z944 versus mZ944", y=1.02)
    fig.text(
        0.5, -0.02,
        "Run 1–3 error bars: temporal SD; Mean error bars: SD among the three replica means.",
        ha="center", fontsize=9,
    )
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"overall_energy_comparison.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_caption(residues, summary):
    lookup = summary.drop_duplicates("resid").set_index("resid")["resname"].to_dict()
    labels = ", ".join(f"{lookup.get(r, 'UNK').title()}{r}" for r in residues)
    text = (
        "Figure 1. Overall MM/GBSA comparison of Z944 and mZ944. Bars show the "
        "electrostatic, van der Waals, and total binding-energy components. For "
        "individual runs, error bars denote temporal standard deviations; for the "
        "three-run mean, error bars denote the standard deviation among replica means.\n\n"
        "Figure 2. Per-residue MM/GBSA decomposition for Z944 and mZ944, averaged "
        "over three independent replicas. Error bars denote the standard deviation "
        "among replica means. Negative energies indicate favorable contributions. "
        "The plotted residues are the union of the 12 most favorable total-energy "
        f"contributors in either system: {labels}.\n\n"
        "Caution: the pairwise GBIS per-residue terms are not strictly additive, so "
        "their sum is not expected to reproduce the whole-complex binding energy.\n"
    )
    (OUT / "figure_captions.txt").write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    per_residue = load_per_residue()
    summary = save_per_residue_tables(per_residue)
    residues = choose_residues(summary)
    plot_per_residue(summary, residues)
    overall = load_overall()
    plot_overall(overall)
    write_caption(residues, summary)
    with (OUT / "selected_residues.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["resid"])
        writer.writerows([[r] for r in residues])
    print(f"Wrote comparison results to {OUT}")
    print(f"Selected {len(residues)} residues: {', '.join(residues)}")


if __name__ == "__main__":
    main()
