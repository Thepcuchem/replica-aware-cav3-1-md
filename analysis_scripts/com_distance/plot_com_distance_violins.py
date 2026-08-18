#!/usr/bin/env python3
"""Compile three-replica COM distances and create manuscript violin plots."""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "com_distance" / "Z944_vs_mZ944"
COLORS = {"Z944": "#2537D8", "mZ944": "#E63737"}
RUN1 = {
    ("Z944", "pocket"): Path(
        "/work1/ted/June2023/WithLigand/md_run1/bombam_analyse/"
        "distance_pocket2lig/distance_pocket2lig.out"
    ),
    ("Z944", "filter"): Path(
        "/work1/ted/June2023/WithLigand/md_run1/bombam_analyse/"
        "distance_sf2lig/distance_sf2lig.out"
    ),
    ("mZ944", "pocket"): Path(
        "/work1/ted/June2023/WithmodLigand/md_run1/bombam_analyse/"
        "distance_pocket2lig/distance_pocket2lig.out"
    ),
    ("mZ944", "filter"): Path(
        "/work1/ted/June2023/WithmodLigand/md_run1/bombam_analyse/"
        "distance_sf2lig/distance_sf2lig.out"
    ),
}


def load_run1():
    records = []
    for system in ("Z944", "mZ944"):
        series = {}
        for metric in ("pocket", "filter"):
            data = pd.read_csv(
                RUN1[(system, metric)], sep=r"\s+", header=None,
                names=["source_frame", "run_label", "distance", "source_system"],
            )
            # The files contain 31,000 frames at 0.01 ns. Retain the final
            # 30,000 frames (300 ns) and sample every 50 frames (0.5 ns).
            data = data.tail(30000).iloc[::50].reset_index(drop=True)
            if len(data) != 600:
                raise ValueError(f"Expected 600 run1 samples, got {len(data)}")
            series[metric] = data["distance"].to_numpy(float)
        for frame in range(600):
            records.append({
                "system": system, "run": 1, "frame": frame,
                "time_ns": frame * 0.5,
                "pocket_distance_A": series["pocket"][frame],
                "filter_distance_A": series["filter"][frame],
            })
    return pd.DataFrame(records)


def load_new_runs():
    tables = []
    for system in ("Z944", "mZ944"):
        for run in (2, 3):
            path = ROOT / f"com_distance/{system}/run{run}/com_distances.csv"
            table = pd.read_csv(path)
            if len(table) != 600:
                raise ValueError(f"{path} contains {len(table)} rather than 600 frames")
            tables.append(table)
    return pd.concat(tables, ignore_index=True)


def compile_data():
    wide = pd.concat([load_run1(), load_new_runs()], ignore_index=True)
    wide = wide.sort_values(["system", "run", "frame"]).reset_index(drop=True)
    long = wide.melt(
        id_vars=["system", "run", "frame", "time_ns"],
        value_vars=["pocket_distance_A", "filter_distance_A"],
        var_name="distance_type", value_name="distance_A",
    )
    long["distance_type"] = long["distance_type"].map({
        "pocket_distance_A": "Blocker-binding pocket",
        "filter_distance_A": "Selectivity filter",
    })
    return wide, long


def violin(ax, values, positions, colors, widths=0.28):
    result = ax.violinplot(
        values, positions=positions, widths=widths,
        showmeans=False, showmedians=False, showextrema=False,
        bw_method="scott", points=250,
    )
    for body, color in zip(result["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_linewidth(0.75)
        body.set_alpha(0.95)
    for vals, pos in zip(values, positions):
        q1, median, q3 = np.percentile(vals, [25, 50, 75])
        mean = np.mean(vals)
        ax.vlines(pos, q1, q3, color="black", linewidth=2.0, zorder=4)
        ax.scatter(pos, median, color="black", s=15, zorder=5)
        ax.scatter(pos, mean, facecolor="white", edgecolor="black",
                   linewidth=0.55, s=27, zorder=6)


def pooled_plot(long):
    metrics = ["Blocker-binding pocket", "Selectivity filter"]
    systems = ["Z944", "mZ944"]
    positions = [0.84, 1.16, 1.84, 2.16]
    values, colors = [], []
    for metric in metrics:
        for system in systems:
            values.append(long[
                (long.distance_type == metric) & (long.system == system)
            ].distance_A.to_numpy())
            colors.append(COLORS[system])

    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    violin(ax, values, positions, colors, widths=0.29)
    ax.set_xticks([1, 2], ["Pocket-to-ligand", "Selectivity filter-to-ligand"])
    ax.set_ylabel("Center-of-mass distance (Å)")
    ax.set_title("Ligand position relative to the blocker-binding cavity")
    ax.grid(axis="y", alpha=0.25, linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        handles=[Patch(facecolor=COLORS[x], edgecolor="black", label=x)
                 for x in systems],
        frameon=False, loc="upper left",
    )
    ax.text(
        0.99, 0.01,
        "Violin: pooled 0.5-ns samples from three replicas\n"
        "White circle: mean; black dot/line: median and interquartile range",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
    )
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"com_distance_violin_pooled.{suffix}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def replica_plot(long):
    metrics = ["Blocker-binding pocket", "Selectivity filter"]
    systems = ["Z944", "mZ944"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), sharey=False)
    for ax, metric in zip(axes, metrics):
        positions, values, colors = [], [], []
        for run in (1, 2, 3):
            for system, offset in zip(systems, (-0.16, 0.16)):
                positions.append(run + offset)
                values.append(long[
                    (long.distance_type == metric)
                    & (long.system == system)
                    & (long.run == run)
                ].distance_A.to_numpy())
                colors.append(COLORS[system])
        violin(ax, values, positions, colors, widths=0.28)
        ax.set_xticks([1, 2, 3], ["Run 1", "Run 2", "Run 3"])
        ax.set_title(metric)
        ax.set_ylabel("Center-of-mass distance (Å)")
        ax.grid(axis="y", alpha=0.25, linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(
        handles=[Patch(facecolor=COLORS[x], edgecolor="black", label=x)
                 for x in systems],
        frameon=False, loc="best",
    )
    fig.suptitle("Replica-resolved ligand COM-distance distributions", y=1.01)
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"com_distance_violin_by_replica.{suffix}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def statistics(long):
    run_stats = (
        long.groupby(["system", "run", "distance_type"])["distance_A"]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
    )
    pooled_stats = (
        long.groupby(["system", "distance_type"])["distance_A"]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
    )
    run_stats.to_csv(OUT / "com_distance_statistics_by_run.csv", index=False)
    pooled_stats.to_csv(OUT / "com_distance_statistics_pooled.csv", index=False)
    return run_stats, pooled_stats


def write_caption():
    text = (
        "Center-of-mass distance distributions between each ligand and the "
        "blocker-binding pocket or selectivity filter. Distances were calculated "
        "from mass-weighted centers of mass using the same selections as the "
        "original run1 analysis. Each replica contributes 600 snapshots sampled "
        "at 0.5-ns intervals over its final 300 ns. The pooled violins therefore "
        "contain equal contributions from all three replicas. White circles show "
        "means, black dots show medians, and black vertical lines show "
        "interquartile ranges.\n"
    )
    (OUT / "figure_caption.txt").write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wide, long = compile_data()
    wide.to_csv(OUT / "com_distances_all_replicas_wide.csv", index=False)
    long.to_csv(OUT / "com_distances_all_replicas_long.csv", index=False)
    statistics(long)
    pooled_plot(long)
    replica_plot(long)
    write_caption()
    print(f"Wrote COM-distance comparison to {OUT}")


if __name__ == "__main__":
    main()
