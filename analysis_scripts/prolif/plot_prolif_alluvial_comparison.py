#!/usr/bin/env python3
"""Three-replica ProLIF residue-rank alluvial comparison for Z944 and mZ944."""

from pathlib import Path
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "prolif" / "Z944_vs_mZ944_alluvial"
INPUTS = {
    "Z944": {
        1: Path("/work1/ted/June2023/WithLigand/md_run1/analyze/ProLif/output/"
                "interaction_fingerprint.csv"),
        2: ROOT / "prolif/Z944/run2/interaction_fingerprint.csv",
        3: ROOT / "prolif/Z944/run3/interaction_fingerprint.csv",
    },
    "mZ944": {
        1: Path("/work1/ted/June2023/WithmodLigand/md_run1/analyze/ProLif/output/"
                "interaction_fingerprint.csv"),
        2: ROOT / "prolif/mZ944/run2/interaction_fingerprint.csv",
        3: ROOT / "prolif/mZ944/run3/interaction_fingerprint.csv",
    },
}
DOMAIN_COLORS = {
    "DI": "#2ca02c", "DII": "#2878B5", "DIII": "#9467bd", "DIV": "#ff7f0e",
}
INTERACTION_LABELS = {
    "Hydrophobic": "Hydrophobic",
    "VdWContact": "VdW",
    "HBDonor": "H-bond donor",
    "HBAcceptor": "H-bond acceptor",
    "PiCation": "Cation–π",
    "CationPi": "Cation–π",
    "PiStacking": "π-stacking",
    "Anionic": "Ionic",
    "Cationic": "Ionic",
}


def parse_residue(protein_label):
    token = protein_label.split(".")[0]
    match = re.fullmatch(r"([A-Za-z]+)(-?\d+)", token)
    if not match:
        raise ValueError(f"Cannot parse residue label: {protein_label}")
    return match.group(1).upper(), int(match.group(2))


def domain_for(resid):
    if 210 <= resid <= 391:
        return "DI"
    if 823 <= resid <= 956:
        return "DII"
    if 1374 <= resid <= 1509:
        return "DIII"
    if 1689 <= resid <= 1833:
        return "DIV"
    return "Other"


def load_run(path):
    fp = pd.read_csv(path, header=[0, 1, 2], index_col=0)
    records = []
    protein_labels = list(dict.fromkeys(fp.columns.get_level_values("protein")))
    for protein in protein_labels:
        resname, resid = parse_residue(protein)
        block = fp.xs(protein, axis=1, level="protein")
        any_contact = block.any(axis=1).mean()
        for interaction in list(dict.fromkeys(block.columns.get_level_values("interaction"))):
            cols = block.xs(interaction, axis=1, level="interaction")
            if isinstance(cols, pd.Series):
                occupancy = cols.mean()
            else:
                occupancy = cols.any(axis=1).mean()
            records.append({
                "resid": resid, "resname": resname, "interaction": interaction,
                "interaction_occupancy": float(occupancy),
                "any_contact": float(any_contact), "frames": len(fp),
            })
    return pd.DataFrame(records)


def aggregate():
    all_runs = []
    for system, runs in INPUTS.items():
        for run, path in runs.items():
            frame = load_run(path)
            frame["system"] = system
            frame["run"] = run
            all_runs.append(frame)
    detail = pd.concat(all_runs, ignore_index=True)

    residue_run = (
        detail.groupby(["system", "run", "resid", "resname"], as_index=False)
        .agg(any_contact=("any_contact", "first"), frames=("frames", "first"))
    )
    # Reindex absent residues to zero so every replica contributes equally.
    residue_rows = []
    interaction_rows = []
    for system in INPUTS:
        residues = (
            detail.loc[detail.system == system, ["resid", "resname"]]
            .drop_duplicates().sort_values("resid")
        )
        interactions = sorted(detail.loc[detail.system == system, "interaction"].unique())
        for run in (1, 2, 3):
            sub = residue_run[(residue_run.system == system) & (residue_run.run == run)]
            merged = residues.merge(sub, on=["resid", "resname"], how="left")
            merged["system"] = system
            merged["run"] = run
            merged["any_contact"] = merged["any_contact"].fillna(0)
            residue_rows.append(merged)

            idata = detail[(detail.system == system) & (detail.run == run)]
            grid = (
                residues.assign(key=1)
                .merge(pd.DataFrame({"interaction": interactions, "key": 1}), on="key")
                .drop(columns="key")
                .merge(
                    idata[["resid", "resname", "interaction", "interaction_occupancy"]],
                    on=["resid", "resname", "interaction"], how="left",
                )
            )
            grid["system"] = system
            grid["run"] = run
            grid["interaction_occupancy"] = grid["interaction_occupancy"].fillna(0)
            interaction_rows.append(grid)

    residue_complete = pd.concat(residue_rows, ignore_index=True)
    interaction_complete = pd.concat(interaction_rows, ignore_index=True)
    summary = (
        residue_complete.groupby(["system", "resid", "resname"], as_index=False)
        .agg(mean_any_contact=("any_contact", "mean"),
             sd_any_contact=("any_contact", "std"))
    )
    itype = (
        interaction_complete.groupby(
            ["system", "resid", "resname", "interaction"], as_index=False
        )["interaction_occupancy"].mean()
    )
    dominant = (
        itype.sort_values(
            ["system", "resid", "interaction_occupancy"],
            ascending=[True, True, False],
        ).drop_duplicates(["system", "resid"])
        .rename(columns={
            "interaction": "dominant_interaction",
            "interaction_occupancy": "dominant_occupancy",
        })
    )
    summary = summary.merge(
        dominant[["system", "resid", "dominant_interaction", "dominant_occupancy"]],
        on=["system", "resid"], how="left",
    )
    summary["domain"] = summary["resid"].map(domain_for)
    summary["residue"] = (
        summary["resname"].str.title() + summary["resid"].astype(str)
    )
    return detail, residue_complete, interaction_complete, summary, itype


def bezier(ax, x0, y0, x1, y1, color, width, linestyle="-", alpha=0.85):
    delta = (x1 - x0) * 0.48
    path = MplPath(
        [(x0, y0), (x0 + delta, y0), (x1 - delta, y1), (x1, y1)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(PathPatch(
        path, facecolor="none", edgecolor=color, linewidth=width,
        linestyle=linestyle, alpha=alpha, capstyle="round",
    ))


def plot_alluvial(summary, top_n=15):
    tops = {}
    for system in ("Z944", "mZ944"):
        tops[system] = (
            summary[summary.system == system]
            .sort_values(["mean_any_contact", "resid"], ascending=[False, True])
            .head(top_n).reset_index(drop=True)
        )
        tops[system]["rank"] = np.arange(1, top_n + 1)

    selected = pd.concat(tops.values(), ignore_index=True)
    selected.to_csv(OUT / "top15_three_replica_residue_contacts.csv", index=False)

    fig, ax = plt.subplots(figsize=(15, 10))
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.8, top_n + 1.2)
    ax.axis("off")
    x_left, x_right = 0.34, 0.66
    node_w, node_h = 0.025, 0.68
    ys = np.arange(top_n - 1, -1, -1)

    positions = {}
    for system, x, side in (("Z944", x_left, "left"), ("mZ944", x_right, "right")):
        table = tops[system]
        for i, row in table.iterrows():
            y = ys[i]
            positions[(system, int(row.resid))] = y
            color = DOMAIN_COLORS.get(row.domain, "#777777")
            ax.add_patch(Rectangle(
                (x - node_w / 2, y - node_h / 2), node_w, node_h,
                facecolor=color, edgecolor="black", linewidth=0.65, zorder=4,
            ))
            interaction = INTERACTION_LABELS.get(
                row.dominant_interaction, row.dominant_interaction
            )
            label = (
                f"{i + 1}. {row.residue}  {100 * row.mean_any_contact:.1f}%"
                f"  [{interaction} {100 * row.dominant_occupancy:.1f}%]"
            )
            if side == "left":
                ax.text(x - 0.022, y, label, ha="right", va="center", fontsize=8.3)
            else:
                ax.text(x + 0.022, y, label, ha="left", va="center", fontsize=8.3)

    left_ids = set(tops["Z944"].resid.astype(int))
    right_ids = set(tops["mZ944"].resid.astype(int))
    shared = left_ids & right_ids
    for resid in sorted(shared):
        left = tops["Z944"].set_index("resid").loc[resid]
        right = tops["mZ944"].set_index("resid").loc[resid]
        strength = (left.mean_any_contact + right.mean_any_contact) / 2
        bezier(
            ax, x_left + node_w / 2, positions[("Z944", resid)],
            x_right - node_w / 2, positions[("mZ944", resid)],
            DOMAIN_COLORS.get(left.domain, "#777777"),
            width=1.2 + 5.2 * strength, alpha=0.78,
        )

    # Dashed stubs identify residues occurring in only one top-15 set.
    for system, ids, x, direction in (
        ("Z944", left_ids - shared, x_left + node_w / 2, 1),
        ("mZ944", right_ids - shared, x_right - node_w / 2, -1),
    ):
        table = tops[system].set_index("resid")
        for resid in sorted(ids):
            row = table.loc[resid]
            y = positions[(system, resid)]
            target_x = 0.50 + direction * 0.015
            target_y = -0.70
            bezier(
                ax, x, y, target_x, target_y,
                DOMAIN_COLORS.get(row.domain, "#777777"),
                width=1.0 + 2.8 * row.mean_any_contact,
                linestyle=(0, (2, 2)), alpha=0.42,
            )

    ax.text(x_left, top_n + 0.35, "Z944", ha="center", va="bottom",
            fontsize=16, fontweight="bold")
    ax.text(x_right, top_n + 0.35, "mZ944", ha="center", va="bottom",
            fontsize=16, fontweight="bold")
    ax.text(
        0.5, top_n + 0.88,
        "Three-replica comparison of residue–ligand interaction networks",
        ha="center", va="bottom", fontsize=17, fontweight="bold",
    )
    ax.text(
        0.5, top_n + 0.48,
        "Ranked by mean ProLIF any-contact occupancy",
        ha="center", va="bottom", fontsize=11,
    )

    legend_y = -1.35
    start_x = 0.25
    for i, (domain, color) in enumerate(DOMAIN_COLORS.items()):
        lx = start_x + i * 0.105
        ax.add_patch(Rectangle((lx, legend_y), 0.015, 0.28,
                               facecolor=color, edgecolor="black", linewidth=0.5))
        ax.text(lx + 0.020, legend_y + 0.14, domain, va="center", fontsize=9)
    ax.plot([0.68, 0.72], [legend_y + 0.14] * 2, color="#666", linewidth=4)
    ax.text(0.725, legend_y + 0.14, "Present in both top-15 lists",
            va="center", fontsize=9)
    ax.plot([0.68, 0.72], [legend_y - 0.28] * 2, color="#888",
            linewidth=2, linestyle=(0, (2, 2)))
    ax.text(0.725, legend_y - 0.28, "Only in one top-15 list",
            va="center", fontsize=9)

    fig.text(
        0.5, 0.012,
        "Node text: three-run mean any-contact occupancy; brackets: dominant "
        "interaction type and its three-run mean occupancy. Link width scales with occupancy.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"prolif_alluvial_Z944_vs_mZ944.{suffix}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_caption():
    caption = (
        "Alluvial comparison of the residue–ligand interaction networks of Z944 "
        "and mZ944. Residues are independently ranked by their mean any-contact "
        "occupancy across three MD replicas. Solid curves connect residues present "
        "in both top-15 lists, whereas dashed curves indicate residues unique to "
        "one top-15 list. Node and curve colors denote channel domains (DI–DIV), "
        "and curve width scales with mean occupancy. The dominant ProLIF interaction "
        "class and its mean occupancy are shown in brackets. Because interaction "
        "classes can coexist in one frame, their occupancies are not additive.\n"
    )
    (OUT / "figure_caption.txt").write_text(caption)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    detail, residue_runs, interaction_runs, summary, interaction_summary = aggregate()
    detail.to_csv(OUT / "interaction_occupancy_by_run.csv", index=False)
    residue_runs.to_csv(OUT / "residue_any_contact_by_run.csv", index=False)
    interaction_runs.to_csv(OUT / "interaction_complete_by_run.csv", index=False)
    summary.to_csv(OUT / "three_replica_residue_summary.csv", index=False)
    interaction_summary.to_csv(OUT / "three_replica_interaction_summary.csv", index=False)
    plot_alluvial(summary)
    write_caption()
    print(f"Wrote alluvial comparison to {OUT}")


if __name__ == "__main__":
    main()
