#!/usr/bin/env python3
"""Integrate reproducible distances with ProLIF, water bridges, and MM/GBSA."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / ".deps"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANALYSIS_ROOT = PROJECT / "processed_data" / "auxiliary_inputs" / "analysis_replicas"
DISTANCE_RESULTS = PROJECT / "results" / "reproducible_distance_determinants"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, default=ANALYSIS_ROOT)
    parser.add_argument("--distance-results", type=Path, default=DISTANCE_RESULTS)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT / "results" / "mechanistic_evidence_integration",
    )
    parser.add_argument("--pairs-per-comparison", type=int, default=25)
    return parser.parse_args()


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


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


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    determinant_rows = read_csv(
        args.distance_results / "top_reproducible_distance_determinants.csv"
    )
    selected_pairs = [
        row for row in determinant_rows if int(row["rank"]) <= args.pairs_per_comparison
    ]

    prolif_rows = read_csv(
        args.analysis_root
        / "prolif/Z944_vs_mZ944_alluvial/three_replica_residue_summary.csv"
    )
    prolif = {
        (row["system"], int(row["resid"])): row
        for row in prolif_rows
    }
    mmgbsa_rows = read_csv(
        args.analysis_root
        / "mmgbsa/Z944_vs_mZ944_comparison/per_residue_three_run_summary.csv"
    )
    mmgbsa = {
        (row["system"], int(row["resid"])): row
        for row in mmgbsa_rows
    }
    water_rows = read_csv(
        args.analysis_root
        / "water_bridges/Z944/three_run_tables/leading_bridging_residues.tsv",
        delimiter="\t",
    )
    water: dict[int, float] = {}
    residue_names: dict[int, str] = {}
    for row in water_rows:
        _, resid, resname = row["residue"].split(":")
        residue = int(resid)
        water[residue] = float(
            np.mean(
                [
                    float(row["run1_percent"]),
                    float(row["run2_percent"]),
                    float(row["run3_percent"]),
                ]
            )
        )
        residue_names[residue] = resname
    for row in prolif_rows:
        residue_names[int(row["resid"])] = row["resname"]
    for row in mmgbsa_rows:
        residue_names[int(row["resid"])] = row["resname"]

    residue_structural_score: dict[int, float] = defaultdict(float)
    residue_comparisons: dict[int, set[str]] = defaultdict(set)
    residue_pair_count: dict[int, int] = defaultdict(int)
    pair_output: list[dict[str, object]] = []
    pattern = re.compile(r"ca_dist_(\d+)_(\d+)_A")
    for row in selected_pairs:
        match = pattern.fullmatch(row["feature"])
        if not match:
            raise ValueError(f"Unexpected feature name: {row['feature']}")
        first, second = map(int, match.groups())
        robustness = float(row["minimum_absolute_standardized_effect"])
        rank_weighted = robustness / np.sqrt(int(row["rank"]))
        for residue in (first, second):
            residue_structural_score[residue] += rank_weighted
            residue_comparisons[residue].add(row["comparison"])
            residue_pair_count[residue] += 1
        pair_output.append(
            {
                **row,
                "residue_1": first,
                "residue_1_name": residue_names.get(first, ""),
                "residue_1_domain": domain(first),
                "residue_2": second,
                "residue_2_name": residue_names.get(second, ""),
                "residue_2_domain": domain(second),
                "interdomain": domain(first) != domain(second),
            }
        )
    write_csv(args.output_dir / "integrated_pair_evidence.csv", pair_output)

    residue_output: list[dict[str, object]] = []
    for residue in sorted(residue_structural_score):
        z_prolif = float(prolif.get(("Z944", residue), {}).get("mean_any_contact", 0))
        m_prolif = float(prolif.get(("mZ944", residue), {}).get("mean_any_contact", 0))
        z_energy = float(mmgbsa.get(("Z944", residue), {}).get("TOTAL_mean", 0))
        m_energy = float(mmgbsa.get(("mZ944", residue), {}).get("TOTAL_mean", 0))
        water_occupancy = water.get(residue, 0.0)
        contact_supported = max(z_prolif, m_prolif) >= 0.10
        energy_supported = min(z_energy, m_energy) <= -0.50
        water_supported = water_occupancy >= 5.0
        direct_layers = sum((contact_supported, energy_supported, water_supported))
        role = "direct/contact-coupled" if direct_layers else "structural-network"
        residue_output.append(
            {
                "resid": residue,
                "resname": residue_names.get(residue, ""),
                "domain": domain(residue),
                "structural_score": residue_structural_score[residue],
                "top_pair_occurrences": residue_pair_count[residue],
                "comparisons_supported": ";".join(sorted(residue_comparisons[residue])),
                "comparison_count": len(residue_comparisons[residue]),
                "Z944_ProLIF_contact": z_prolif,
                "mZ944_ProLIF_contact": m_prolif,
                "Z944_MMGBSA_total_kcal_mol": z_energy,
                "mZ944_MMGBSA_total_kcal_mol": m_energy,
                "Z944_water_bridge_mean_percent": water_occupancy,
                "contact_supported": contact_supported,
                "energy_supported": energy_supported,
                "water_supported": water_supported,
                "direct_evidence_layers": direct_layers,
                "mechanistic_role": role,
            }
        )
    residue_output.sort(
        key=lambda row: (
            int(row["direct_evidence_layers"]),
            int(row["comparison_count"]),
            float(row["structural_score"]),
        ),
        reverse=True,
    )
    write_csv(args.output_dir / "integrated_residue_evidence.csv", residue_output)

    shown = residue_output[:20]
    metrics = np.asarray(
        [
            [
                float(row["structural_score"]),
                max(float(row["Z944_ProLIF_contact"]), float(row["mZ944_ProLIF_contact"])),
                max(0.0, -min(float(row["Z944_MMGBSA_total_kcal_mol"]),
                              float(row["mZ944_MMGBSA_total_kcal_mol"]))),
                float(row["Z944_water_bridge_mean_percent"]) / 100.0,
            ]
            for row in shown
        ]
    )
    normalized = metrics / np.maximum(metrics.max(axis=0, keepdims=True), 1e-12)
    figure, axis = plt.subplots(figsize=(9, 9), constrained_layout=True)
    image = axis.imshow(normalized, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    axis.set_xticks(
        range(4),
        ("Structural\nreproducibility", "ProLIF\ncontact", "Favorable\nMM/GBSA", "Z944 water\nbridge"),
    )
    axis.set_yticks(
        range(len(shown)),
        [
            f"{row['resname'].title()}{row['resid']} ({row['domain']})"
            if row["resname"]
            else f"{row['resid']} ({row['domain']})"
            for row in shown
        ],
    )
    axis.set_title("Integrated evidence for reproducible structural determinants")
    figure.colorbar(image, ax=axis, label="Within-layer normalized evidence")
    figure.savefig(args.output_dir / "integrated_mechanistic_evidence.png", dpi=300)
    figure.savefig(args.output_dir / "integrated_mechanistic_evidence.pdf")
    plt.close(figure)

    summary = {
        "top_pairs_integrated_per_comparison": args.pairs_per_comparison,
        "unique_residues": len(residue_output),
        "direct_contact_coupled_residues": sum(
            row["mechanistic_role"] == "direct/contact-coupled"
            for row in residue_output
        ),
        "structural_network_residues": sum(
            row["mechanistic_role"] == "structural-network"
            for row in residue_output
        ),
        "evidence_thresholds": {
            "ProLIF_mean_contact": 0.10,
            "MMGBSA_total_kcal_mol": -0.50,
            "Z944_water_bridge_mean_percent": 5.0,
        },
        "scope": (
            "Water-bridge evidence is available only for Z944. MM/GBSA and ProLIF "
            "apply only to ligand-bound systems."
        ),
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
