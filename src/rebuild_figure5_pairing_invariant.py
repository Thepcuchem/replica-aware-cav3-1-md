#!/usr/bin/env python3
"""Rebuild Figure 5A from all nine cross-replica standardized effects."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / "processed_data/common_ca_distances"
PAIRING = ROOT / "results/pairing"
OUT = ROOT / "results/figure5"
CONTRASTS = (("Z944", "Apo"), ("mZ944", "Apo"), ("mZ944", "Z944"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    means, sds = {}, {}
    names = None
    for system in ("Apo", "Z944", "mZ944"):
        for run in (1, 2, 3):
            with np.load(FEATURES / f"{system.lower()}_run{run}_common_ca_distances.npz") as data:
                current = data["feature_names"].astype(str)
                values = data["distances"].astype(float)
            if names is None:
                names = current
            elif not np.array_equal(names, current):
                raise ValueError("Feature-name mismatch")
            means[(system, run)] = values.mean(axis=0)
            sds[(system, run)] = values.std(axis=0, ddof=1)
    index = {name: i for i, name in enumerate(names)}
    columns = [f"{i}–{j}" for i in (1, 2, 3) for j in (1, 2, 3)]
    fig, axes = plt.subplots(1, 3, figsize=(18, 8.4), constrained_layout=True)
    output_rows = []
    for axis, (numerator, denominator) in zip(axes, CONTRASTS):
        label = f"{numerator}_vs_{denominator}"
        selected = read_csv(PAIRING / f"{label}_top_pairing_invariant_features.csv")[:15]
        effects = []
        for rank, row in enumerate(selected, start=1):
            k = index[row["feature"]]
            values = []
            for i in (1, 2, 3):
                for j in (1, 2, 3):
                    pooled = np.sqrt((sds[(numerator, i)][k] ** 2 + sds[(denominator, j)][k] ** 2) / 2)
                    value = (means[(numerator, i)][k] - means[(denominator, j)][k]) / pooled
                    values.append(float(value))
                    output_rows.append({"comparison": label, "rank": rank,
                                        "feature": row["feature"],
                                        "numerator_run": i, "denominator_run": j,
                                        "standardized_effect": float(value)})
            effects.append(values)
        matrix = np.asarray(effects)
        limit = max(1.0, float(np.max(np.abs(matrix))))
        image = axis.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
        axis.set_xticks(range(9), columns, rotation=45, ha="right")
        axis.set_yticks(range(len(selected)),
                        [r["feature"].replace("ca_dist_", "").replace("_A", "").replace("_", "–")
                         for r in selected])
        axis.set_xlabel("Numerator run–denominator run")
        axis.set_title(label.replace("_vs_", " vs "))
        fig.colorbar(image, ax=axis, shrink=.72, label="Standardized effect")
    fig.suptitle(
        "Top pairing-invariant determinants across all nine cross-replica comparisons",
        fontweight="bold",
    )
    fig.savefig(OUT / "pairing_invariant_all9_effects.png", dpi=300)
    fig.savefig(OUT / "pairing_invariant_all9_effects.pdf")
    plt.close(fig)
    with (OUT / "panelA_all9_standardized_effects.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader(); writer.writerows(output_rows)


if __name__ == "__main__":
    main()
