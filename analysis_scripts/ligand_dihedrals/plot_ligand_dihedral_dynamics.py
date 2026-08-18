#!/usr/bin/env python3
"""Compile and plot the seven DZR dihedrals across three MD replicas."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "ligand_dihedrals" / "Z944_vs_mZ944"
COLORS = {"Z944": "#2138D8", "mZ944": "#ED3030"}
CHI = [f"chi{i}_deg" for i in range(1, 8)]
DEFINITIONS = {
    "chi1_deg": "C9–N3–C8–C7",
    "chi2_deg": "N3–C8–C7–N1",
    "chi3_deg": "C8–C7–N1–C4",
    "chi4_deg": "C3–C1–C6–N2",
    "chi5_deg": "C1–C6–N2–C10",
    "chi6_deg": "C6–N2–C10–C14",
    "chi7_deg": "N2–C10–C14–C16",
}
RUN1 = {
    "Z944": Path(
        "/work1/ted/June2023/WithLigand/md_run1/bombam_analyse/"
        "dihedral_angle/dihedral_angle_lig_all.out"
    ),
    "mZ944": Path(
        "/work1/ted/June2023/WithmodLigand/md_run1/bombam_analyse/"
        "dihedral_angle/dihedral_angle_lig_all.out"
    ),
}


def load_run1():
    tables = []
    columns = ["source_frame", "source_time"] + CHI
    for system, path in RUN1.items():
        data = pd.read_csv(path, sep=r"\s+", header=None, names=columns)
        # The 31,000 values cover 310 ns at 0.01 ns/frame; use the final
        # 30,000 frames and retain every 50th frame for 0.5-ns sampling.
        data = data.tail(30000).iloc[::50].reset_index(drop=True)
        if len(data) != 600:
            raise ValueError(f"Expected 600 run1 samples in {path}, found {len(data)}")
        data["system"] = system
        data["run"] = 1
        data["frame"] = np.arange(600)
        data["time_ns"] = np.arange(600) * 0.5
        tables.append(data[["system", "run", "frame", "time_ns"] + CHI])
    return pd.concat(tables, ignore_index=True)


def load_new_runs():
    tables = []
    for system in ("Z944", "mZ944"):
        for run in (2, 3):
            path = ROOT / f"ligand_dihedrals/{system}/run{run}/dihedrals.csv"
            table = pd.read_csv(path)
            if len(table) != 600:
                raise ValueError(f"Expected 600 frames in {path}, found {len(table)}")
            tables.append(table)
    return pd.concat(tables, ignore_index=True)


def circular_stats(values):
    radians = np.deg2rad(values)
    vector = np.mean(np.exp(1j * radians))
    mean = np.rad2deg(np.angle(vector))
    resultant = np.abs(vector)
    circ_sd = np.rad2deg(np.sqrt(max(0, -2 * np.log(max(resultant, 1e-12)))))
    return mean, circ_sd, resultant


def compile_stats(data):
    rows = []
    for (system, run), group in data.groupby(["system", "run"]):
        for chi in CHI:
            mean, sd, resultant = circular_stats(group[chi].to_numpy())
            rows.append({
                "system": system, "run": run, "dihedral": chi.replace("_deg", ""),
                "atom_definition": DEFINITIONS[chi], "frames": len(group),
                "circular_mean_deg": mean, "circular_sd_deg": sd,
                "mean_resultant_length": resultant,
            })
    by_run = pd.DataFrame(rows)
    pooled_rows = []
    for system, group in data.groupby("system"):
        for chi in CHI:
            mean, sd, resultant = circular_stats(group[chi].to_numpy())
            pooled_rows.append({
                "system": system, "dihedral": chi.replace("_deg", ""),
                "atom_definition": DEFINITIONS[chi], "frames": len(group),
                "circular_mean_deg": mean, "circular_sd_deg": sd,
                "mean_resultant_length": resultant,
            })
    pooled = pd.DataFrame(pooled_rows)
    by_run.to_csv(OUT / "dihedral_circular_statistics_by_run.csv", index=False)
    pooled.to_csv(OUT / "dihedral_circular_statistics_pooled.csv", index=False)


def timeline_plot(data):
    fig, axes = plt.subplots(
        7, 3, figsize=(16, 17.5), sharex=True, sharey=True,
        gridspec_kw={"hspace": 0.13, "wspace": 0.10},
    )
    for row, chi in enumerate(CHI):
        for col, run in enumerate((1, 2, 3)):
            ax = axes[row, col]
            for system in ("Z944", "mZ944"):
                part = data[(data.system == system) & (data.run == run)]
                angles = part[chi].to_numpy(float).copy()
                # Do not draw artificial vertical lines when a continuous
                # torsion crosses the conventional -180/180-degree boundary.
                wrap = np.where(np.abs(np.diff(angles)) > 180)[0] + 1
                angles[wrap] = np.nan
                ax.plot(
                    part.time_ns, angles, color=COLORS[system],
                    linewidth=0.72, alpha=0.88, label=system,
                )
            ax.set_ylim(-185, 185)
            ax.set_yticks([-180, -90, 0, 90, 180])
            ax.grid(alpha=0.22, linewidth=0.5)
            ax.spines[["top", "right"]].set_visible(False)
            if row == 0:
                ax.set_title(f"Run {run}", fontsize=12)
            if col == 0:
                ax.set_ylabel(
                    f"$\\chi_{{{row + 1}}}$ (°)\n{DEFINITIONS[chi]}",
                    fontsize=9,
                )
            if row == 6:
                ax.set_xlabel("Time within final 300 ns (ns)")
    axes[0, 0].legend(frameon=False, ncol=2, loc="lower left")
    fig.suptitle(
        "Ligand dihedral-angle dynamics during the final 300 ns",
        fontsize=17, y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.987))
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"dihedral_timelines_all_replicas.{suffix}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def circular_density(values, bin_width=5):
    edges = np.arange(-180, 180 + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    # Circular Gaussian smoothing with sigma=2 bins.
    radius, sigma = 8, 2.0
    offsets = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel /= kernel.sum()
    padded = np.concatenate([counts[-radius:], counts, counts[:radius]])
    smooth = np.convolve(padded, kernel, mode="same")[radius:-radius]
    return centers, smooth


def distribution_plot(data, by_replica=False):
    fig, axes = plt.subplots(4, 2, figsize=(11.5, 13), sharex=True)
    axes = axes.ravel()
    linestyles = {1: "-", 2: "--", 3: ":"}
    for i, chi in enumerate(CHI):
        ax = axes[i]
        if by_replica:
            for system in ("Z944", "mZ944"):
                for run in (1, 2, 3):
                    vals = data[(data.system == system) & (data.run == run)][chi]
                    x, y = circular_density(vals.to_numpy())
                    ax.plot(
                        x, y, color=COLORS[system], linestyle=linestyles[run],
                        linewidth=1.5, alpha=0.90,
                        label=f"{system} run{run}",
                    )
        else:
            for system in ("Z944", "mZ944"):
                vals = data[data.system == system][chi]
                x, y = circular_density(vals.to_numpy())
                ax.plot(x, y, color=COLORS[system], linewidth=2.0, label=system)
                ax.fill_between(x, 0, y, color=COLORS[system], alpha=0.10)
        ax.set_xlim(-180, 180)
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_title(f"$\\chi_{{{i + 1}}}$: {DEFINITIONS[chi]}", fontsize=11)
        ax.set_ylabel("Probability density")
        ax.grid(alpha=0.22, linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        if i >= 5:
            ax.set_xlabel("Dihedral angle (°)")
    axes[7].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[7].legend(handles, labels, loc="center", frameon=False, ncol=1)
    title = (
        "Replica-resolved ligand dihedral distributions"
        if by_replica else
        "Three-replica pooled ligand dihedral distributions"
    )
    fig.suptitle(title, fontsize=17, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    stem = "dihedral_distributions_by_replica" if by_replica else "dihedral_distributions_pooled"
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_caption():
    text = (
        "Ligand conformational dynamics during the final 300 ns of three "
        "independent MD replicas. Seven dihedral angles were defined exactly as "
        "in the original Figure 5 analysis. Each replica contributes 600 snapshots "
        "at 0.5-ns intervals. Timeline panels compare Z944 and mZ944 within each "
        "replica. Distribution curves use periodic 5-degree histograms with "
        "circular Gaussian smoothing; pooled distributions contain equal numbers "
        "of observations from each replica. Angular summary values are circular "
        "means and circular standard deviations because arithmetic statistics are "
        "inappropriate at the -180/180-degree boundary.\n"
    )
    (OUT / "figure_caption.txt").write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.concat([load_run1(), load_new_runs()], ignore_index=True)
    data = data.sort_values(["system", "run", "frame"]).reset_index(drop=True)
    data.to_csv(OUT / "dihedral_angles_all_replicas.csv", index=False)
    compile_stats(data)
    timeline_plot(data)
    distribution_plot(data, by_replica=False)
    distribution_plot(data, by_replica=True)
    write_caption()
    print(f"Wrote ligand-dihedral analysis to {OUT}")


if __name__ == "__main__":
    main()
