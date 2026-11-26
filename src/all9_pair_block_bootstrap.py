#!/usr/bin/env python3
"""Block-bootstrap all nine cross-system trajectory pairs for revised determinants."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / "processed_data/common_ca_distances"
PAIRING = ROOT / "results/pairing"
OUT = ROOT / "results/original_all9"
CONTRASTS = (("Z944", "Apo"), ("mZ944", "Apo"), ("mZ944", "Z944"))
RUNS = (1, 2, 3)
BLOCK_FRAMES = 50  # 10 ns at 0.2 ns/frame
ITERATIONS = 2000
SEED = 2026


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def block_means(values: np.ndarray) -> np.ndarray:
    complete = len(values) // BLOCK_FRAMES
    return values[:complete * BLOCK_FRAMES].reshape(complete, BLOCK_FRAMES, -1).mean(axis=1)


def bootstrap_pair(a: np.ndarray, b: np.ndarray, rng: np.random.Generator):
    ia = rng.integers(0, len(a), size=(ITERATIONS, len(a)))
    ib = rng.integers(0, len(b), size=(ITERATIONS, len(b)))
    diff = a[ia].mean(axis=1) - b[ib].mean(axis=1)
    return diff.mean(axis=0), np.quantile(diff, .025, axis=0), np.quantile(diff, .975, axis=0)


def welch(first: np.ndarray, second: np.ndarray):
    n1, n2 = len(first), len(second)
    v1, v2 = first.var(axis=0, ddof=1), second.var(axis=0, ddof=1)
    se2 = v1/n1 + v2/n2
    df = se2**2 / ((v1/n1)**2/(n1-1) + (v2/n2)**2/(n2-1))
    mean = first.mean(axis=0) - second.mean(axis=0)
    half = t.ppf(.975, df) * np.sqrt(se2)
    return mean, mean-half, mean+half, df


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    detail, summary = [], []
    for numerator, denominator in CONTRASTS:
        label = f"{numerator}_vs_{denominator}"
        selected = read_csv(PAIRING / f"{label}_top_pairing_invariant_features.csv")[:25]
        selected_names = [r["feature"] for r in selected]
        blocks, run_means = {}, {}
        for system in (numerator, denominator):
            for run in RUNS:
                with np.load(FEATURES / f"{system.lower()}_run{run}_common_ca_distances.npz") as data:
                    names = data["feature_names"].astype(str)
                    idx = np.asarray([np.where(names == f)[0][0] for f in selected_names])
                    values = data["distances"][:, idx].astype(float)
                blocks[(system, run)] = block_means(values)
                run_means[(system, run)] = values.mean(axis=0)
        pair_results = {}
        for i in RUNS:
            for j in RUNS:
                estimate, low, high = bootstrap_pair(blocks[(numerator, i)], blocks[(denominator, j)], rng)
                pair_results[(i, j)] = (estimate, low, high)
                for k, feature in enumerate(selected_names):
                    detail.append({"comparison": label, "rank": k+1, "feature": feature,
                                   "numerator_run": i, "denominator_run": j,
                                   "mean_difference_A": run_means[(numerator, i)][k]-run_means[(denominator, j)][k],
                                   "bootstrap_mean_difference_A": estimate[k],
                                   "block_bootstrap_95_CI_low_A": low[k],
                                   "block_bootstrap_95_CI_high_A": high[k],
                                   "CI_excludes_zero": bool(low[k] > 0 or high[k] < 0)})
        first = np.vstack([run_means[(numerator, i)] for i in RUNS])
        second = np.vstack([run_means[(denominator, j)] for j in RUNS])
        mean, wlow, whigh, df = welch(first, second)
        for k, source in enumerate(selected):
            lows = np.asarray([pair_results[p][1][k] for p in pair_results])
            highs = np.asarray([pair_results[p][2][k] for p in pair_results])
            point_positive = float(source["minimum_pair_difference_A"]) > 0
            all_supported = bool(np.all(lows > 0) if point_positive else np.all(highs < 0))
            welch_supported = bool(wlow[k] > 0 or whigh[k] < 0)
            summary.append({"comparison": label, "rank": k+1, "feature": source["feature"],
                            "system_mean_difference_A": mean[k],
                            "unpaired_Welch_95_CI_low_A": wlow[k],
                            "unpaired_Welch_95_CI_high_A": whigh[k], "Welch_df": df[k],
                            "all_9_point_estimates_same_direction": True,
                            "all_9_block_bootstrap_CIs_exclude_zero_same_direction": all_supported,
                            "number_of_pair_CIs_excluding_zero": int(sum((pair_results[p][1][k] > 0 or pair_results[p][2][k] < 0) for p in pair_results)),
                            "Welch_CI_excludes_zero": welch_supported,
                            "all9_bootstrap_and_Welch_supported": all_supported and welch_supported,
                            "minimum_pair_difference_A": source["minimum_pair_difference_A"],
                            "maximum_pair_difference_A": source["maximum_pair_difference_A"],
                            "minimum_absolute_standardized_effect_all_9_pairs": source["minimum_absolute_standardized_effect_all_9_pairs"]})
    write_csv(OUT / "all9_pair_block_bootstrap_intervals.csv", detail)
    write_csv(OUT / "all9_pair_block_bootstrap_summary.csv", summary)

    # Reserve a dedicated lower band for the shared legend so it cannot overlap
    # the two-line x-axis titles after the panel is scaled for publication.
    fig, axes = plt.subplots(1, 3, figsize=(17, 9))
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.22, top=0.86, wspace=0.26)
    for ax, (numerator, denominator) in zip(axes, CONTRASTS):
        label = f"{numerator}_vs_{denominator}"
        rows = [r for r in summary if r["comparison"] == label][:15][::-1]
        y = np.arange(len(rows)); means = np.array([r["system_mean_difference_A"] for r in rows])
        lo = np.array([r["unpaired_Welch_95_CI_low_A"] for r in rows]); hi = np.array([r["unpaired_Welch_95_CI_high_A"] for r in rows])
        colors = ["#20854E" if r["all9_bootstrap_and_Welch_supported"] else
                  ("#E18727" if r["all_9_block_bootstrap_CIs_exclude_zero_same_direction"] else "#7F7F7F") for r in rows]
        for yy, mm, ll, hh, color in zip(y, means, lo, hi, colors):
            ax.errorbar(mm, yy, xerr=([mm-ll], [hh-mm]), fmt="none",
                        ecolor=color, elinewidth=2, capsize=3)
        ax.scatter(means, y, c=colors, s=42, zorder=3); ax.axvline(0, color="black", lw=.8)
        ax.set_yticks(y, [r["feature"].replace("ca_dist_", "").replace("_A", "").replace("_", "–") for r in rows])
        ax.set_title(label.replace("_vs_", " vs ")); ax.set_xlabel("Mean distance difference (Å)\nunpaired 95% Welch interval")
    fig.suptitle("Uncertainty of pairing-invariant determinants", fontweight="bold")
    fig.legend(handles=[
        Line2D([0], [0], marker="o", color="#20854E", lw=2,
               label="All 9 block-bootstrap CIs and Welch CI exclude zero"),
        Line2D([0], [0], marker="o", color="#E18727", lw=2,
               label="All 9 block-bootstrap CIs exclude zero; Welch CI includes zero"),
        Line2D([0], [0], marker="o", color="#7F7F7F", lw=2,
               label="At least one pairwise block-bootstrap CI includes zero"),
    ], loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=3,
       frameon=False, fontsize=11)
    fig.savefig(OUT / "pairing_invariant_uncertainty_forest.png", dpi=300)
    fig.savefig(OUT / "pairing_invariant_uncertainty_forest.pdf"); plt.close(fig)
    counts = {f"{a}_vs_{b}": {"all9_block_CIs": sum(r["all_9_block_bootstrap_CIs_exclude_zero_same_direction"] for r in summary if r["comparison"] == f"{a}_vs_{b}"),
                              "all9_block_and_Welch": sum(r["all9_bootstrap_and_Welch_supported"] for r in summary if r["comparison"] == f"{a}_vs_{b}")} for a,b in CONTRASTS}
    (OUT / "analysis_summary.json").write_text(json.dumps({"block_length_ns":10,"iterations":ITERATIONS,"features_per_contrast":25,"counts":counts}, indent=2)+"\n")
    print(json.dumps(counts, indent=2))

if __name__ == "__main__": main()
