#!/usr/bin/env python3
"""Map reproducible distance determinants to PDB B-factors and network figures."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / ".deps"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

COMPARISONS = ("Z944_vs_Apo", "mZ944_vs_Apo", "mZ944_vs_Z944")
REFERENCE_PDB = {
    "Z944_vs_Apo": PROJECT / "z944/z944_dl.pdb",
    "mZ944_vs_Apo": PROJECT / "mz944/mz944_dl.pdb",
    "mZ944_vs_Z944": PROJECT / "mz944/mz944_dl.pdb",
}
DOMAIN_COLORS = {
    "DI": "#4C78A8",
    "DII": "#F58518",
    "DIII": "#54A24B",
    "DIV": "#B279A2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--determinants",
        type=Path,
        default=PROJECT
        / "results/reproducible_distance_determinants/"
        "top_reproducible_distance_determinants.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "structural_determinant_mapping",
    )
    parser.add_argument("--top-pairs", type=int, default=25)
    return parser.parse_args()


def domain(resid: int) -> str:
    if 210 <= resid <= 397:
        return "DI"
    if 863 <= resid <= 968:
        return "DII"
    if 1388 <= resid <= 1516:
        return "DIII"
    if 1717 <= resid <= 1832:
        return "DIV"
    return "Other"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_annotated_pdb(
    source: Path, destination: Path, scores: dict[int, float]
) -> None:
    output: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 66:
            try:
                resid = int(line[22:26])
            except ValueError:
                output.append(line)
                continue
            score = max(-99.99, min(99.99, scores.get(resid, 0.0)))
            line = f"{line[:60]}{score:6.2f}{line[66:]}"
        output.append(line)
    destination.write_text("".join(output), encoding="utf-8")


def write_pymol_script(
    path: Path,
    pdb_name: str,
    comparison: str,
    pairs: list[dict[str, object]],
    residues: list[int],
) -> None:
    commands = [
        f"load {pdb_name}, mapped",
        "hide everything, mapped",
        "show cartoon, mapped and protein",
        "color gray80, mapped and protein",
        "spectrum b, blue_white_red, mapped and protein, minimum=-100, maximum=100",
        f"select determinants, mapped and resid {'+'.join(map(str, residues))}",
        "show sticks, determinants",
        "set stick_radius, 0.18, determinants",
        "show spheres, mapped and resname DZR",
        "set sphere_scale, 0.3, mapped and resname DZR",
        "color yellow, mapped and resname DZR",
    ]
    for index, row in enumerate(pairs[:10], start=1):
        first = int(row["residue_1"])
        second = int(row["residue_2"])
        name = f"pair_{index}_{first}_{second}"
        commands.append(
            f"distance {name}, mapped and name CA and resid {first}, "
            f"mapped and name CA and resid {second}"
        )
        commands.append(f"set dash_width, 2.5, {name}")
        commands.append(
            f"color {'red' if float(row['mean_difference_A']) > 0 else 'blue'}, {name}"
        )
    commands.extend(
        [
            "bg_color white",
            "set ray_opaque_background, off",
            "orient determinants",
            f"set_name mapped, {comparison}",
            f"save {comparison}_determinants.pse",
        ]
    )
    path.write_text("\n".join(commands) + "\n", encoding="utf-8")


def write_vmd_script(
    path: Path,
    pdb_name: str,
    pairs: list[dict[str, object]],
    residues: list[int],
) -> None:
    commands = [
        f"mol new {pdb_name} type pdb waitfor all",
        "set molid [molinfo top]",
        "mol delrep 0 $molid",
        "color scale method BWR",
        "mol representation NewCartoon",
        "mol color Beta",
        "mol selection {protein}",
        "mol material Opaque",
        "mol addrep $molid",
        "mol scaleminmax $molid 0 -100 100",
        "mol representation Licorice 0.18 12.0 12.0",
        "mol color Beta",
        f"mol selection {{protein and resid {' '.join(map(str, residues))}}}",
        "mol material Opaque",
        "mol addrep $molid",
        "mol scaleminmax $molid 1 -100 100",
        "mol representation Licorice 0.25 12.0 12.0",
        "mol color ColorID 4",
        "mol selection {resname DZR}",
        "mol material Opaque",
        "mol addrep $molid",
        "graphics $molid materials on",
    ]
    for index, row in enumerate(pairs[:10], start=1):
        first = int(row["residue_1"])
        second = int(row["residue_2"])
        color = "red" if float(row["mean_difference_A"]) > 0 else "blue"
        commands.extend(
            [
                f"set sel1_{index} [atomselect $molid {{protein and name CA and resid {first}}}]",
                f"set sel2_{index} [atomselect $molid {{protein and name CA and resid {second}}}]",
                f"set p1_{index} [lindex [$sel1_{index} get {{x y z}}] 0]",
                f"set p2_{index} [lindex [$sel2_{index} get {{x y z}}] 0]",
                f"graphics $molid color {color}",
                f"graphics $molid line $p1_{index} $p2_{index} width 3 style dashed",
                f"$sel1_{index} delete",
                f"$sel2_{index} delete",
            ]
        )
    commands.extend(
        [
            "display projection Orthographic",
            "color Display Background white",
            "display resetview",
            "puts {Blue = reproducible contraction; red = reproducible expansion}",
            "quit",
        ]
    )
    path.write_text("\n".join(commands) + "\n", encoding="utf-8")


def network_positions(residues: list[int]) -> dict[int, tuple[float, float]]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for residue in residues:
        grouped[domain(residue)].append(residue)
    positions: dict[int, tuple[float, float]] = {}
    domain_x = {"DI": 0.0, "DII": 1.0, "DIII": 2.0, "DIV": 3.0, "Other": 4.0}
    for label, members in grouped.items():
        members = sorted(members)
        offsets = np.linspace(-1.0, 1.0, len(members)) if len(members) > 1 else [0.0]
        for residue, y in zip(members, offsets):
            positions[residue] = (domain_x[label], float(y))
    return positions


def plot_network(
    comparison: str,
    pairs: list[dict[str, object]],
    residue_strength: dict[int, float],
    output: Path,
) -> None:
    residues = sorted(residue_strength)
    positions = network_positions(residues)
    figure, axis = plt.subplots(figsize=(12, 8), constrained_layout=True)
    maximum_effect = max(abs(float(row["mean_standardized_effect"])) for row in pairs)
    for row in reversed(pairs):
        first, second = int(row["residue_1"]), int(row["residue_2"])
        x1, y1 = positions[first]
        x2, y2 = positions[second]
        effect = float(row["mean_standardized_effect"])
        axis.plot(
            (x1, x2),
            (y1, y2),
            color="#D62728" if effect > 0 else "#1F77B4",
            alpha=0.25 + 0.6 * abs(effect) / maximum_effect,
            linewidth=0.7 + 3.0 * abs(effect) / maximum_effect,
            zorder=1,
        )
    max_strength = max(residue_strength.values())
    for residue in residues:
        x, y = positions[residue]
        label = domain(residue)
        axis.scatter(
            x,
            y,
            s=100 + 650 * residue_strength[residue] / max_strength,
            color=DOMAIN_COLORS.get(label, "#888888"),
            edgecolor="white",
            linewidth=0.8,
            zorder=2,
        )
        axis.text(x + 0.035, y, str(residue), va="center", fontsize=8, zorder=3)
    axis.set_xticks((0, 1, 2, 3), ("DI", "DII", "DIII", "DIV"))
    axis.set_xlim(-0.25, 3.45)
    axis.set_ylim(-1.2, 1.2)
    axis.set_ylabel("Residues arranged within domain")
    axis.set_title(
        f"{comparison.replace('_', ' ')} reproducible distance network\n"
        "red = expansion; blue = contraction"
    )
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", left=False, labelleft=False)
    axis.legend(
        handles=[
            Line2D([0], [0], color="#D62728", lw=3, label="Expanded distance"),
            Line2D([0], [0], color="#1F77B4", lw=3, label="Contracted distance"),
        ],
        frameon=False,
        loc="upper right",
    )
    figure.savefig(output, dpi=300)
    figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.determinants.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    pattern = re.compile(r"ca_dist_(\d+)_(\d+)_A")
    summary: dict[str, object] = {}
    for comparison in COMPARISONS:
        selected = [
            row
            for row in rows
            if row["comparison"] == comparison and int(row["rank"]) <= args.top_pairs
        ]
        pairs: list[dict[str, object]] = []
        signed_scores: dict[int, float] = defaultdict(float)
        strengths: dict[int, float] = defaultdict(float)
        for row in selected:
            match = pattern.fullmatch(row["feature"])
            if not match:
                raise ValueError(row["feature"])
            first, second = map(int, match.groups())
            signed = float(row["mean_standardized_effect"]) / math.sqrt(int(row["rank"]))
            strength = abs(signed)
            for residue in (first, second):
                signed_scores[residue] += signed
                strengths[residue] += strength
            pairs.append({**row, "residue_1": first, "residue_2": second})
        scale = max(abs(value) for value in signed_scores.values())
        normalized = {residue: 100.0 * value / scale for residue, value in signed_scores.items()}
        residue_rows = [
            {
                "comparison": comparison,
                "resid": residue,
                "domain": domain(residue),
                "signed_mapping_score": signed_scores[residue],
                "normalized_B_factor": normalized[residue],
                "structural_strength": strengths[residue],
            }
            for residue in sorted(signed_scores)
        ]
        write_csv(args.output_dir / f"{comparison}_residue_mapping.csv", residue_rows)
        pdb_path = args.output_dir / f"{comparison}_determinants.pdb"
        write_annotated_pdb(REFERENCE_PDB[comparison], pdb_path, normalized)
        write_pymol_script(
            args.output_dir / f"{comparison}_determinants.pml",
            pdb_path.name,
            comparison,
            pairs,
            sorted(signed_scores),
        )
        write_vmd_script(
            args.output_dir / f"{comparison}_determinants.vmd.tcl",
            pdb_path.name,
            pairs,
            sorted(signed_scores),
        )
        plot_network(
            comparison,
            pairs,
            strengths,
            args.output_dir / f"{comparison}_distance_network.png",
        )
        summary[comparison] = {
            "pairs_mapped": len(pairs),
            "residues_mapped": len(signed_scores),
            "reference_pdb": str(REFERENCE_PDB[comparison]),
            "positive_B_factor": "net expansion",
            "negative_B_factor": "net contraction",
        }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
