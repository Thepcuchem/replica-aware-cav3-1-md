#!/usr/bin/env python3
"""Compile Z944 run1/run2/run3 water-bridge results into tables."""
import csv
from pathlib import Path


WORK = Path(__file__).resolve().parent
RUNS = {
    "run1": Path("/work1/ted/June2023/WithLigand/md_run1/analyze/water-mediated/results"),
    "run2": WORK / "water_bridges/Z944/run2",
    "run3": WORK / "water_bridges/Z944/run3",
}
META = {
    "run1": ("302.0-612.0", "md-59--md-89", 0.1),
    "run2": ("200.0-497.0", "md-16--md-30 (partial)", 0.2),
    "run3": ("200.0-497.0", "md-16--md-30 (partial)", 0.2),
}
OUT = WORK / "water_bridges/Z944/three_run_tables"
OUT.mkdir(parents=True, exist_ok=True)


def summary_values(path):
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            values[key.strip()] = value.strip()
    return values


overall = []
residues = {}
for run, directory in RUNS.items():
    summary = summary_values(directory / "run_summary.txt")
    with (directory / "bridge_events.tsv").open(encoding="utf-8") as source:
        event_count = sum(1 for _ in source) - 1
    interval, dcds, sampling = META[run]
    overall.append({
        "run": run,
        "interval": interval,
        "dcds": dcds,
        "sampling": sampling,
        "frames": int(summary["Sampled frames"]),
        "bridge_frames": int(summary["Frames with >=1 bridge"]),
        "occupancy": float(summary["Any-bridge occupancy (%)"]),
        "events": event_count,
    })
    with (directory / "residue_bridge_occupancy.tsv").open(encoding="utf-8") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            key = f'{row["protein_segname"]}:{row["resid"] if "resid" in row else row["protein_resid"]}:{row["resname"] if "resname" in row else row["protein_resname"]}'
            residues.setdefault(key, {})[run] = {
                "frames": int(row["frames"]),
                "occupancy": float(row["occupancy_percent"]),
            }

with (OUT / "overall_water_bridge_results.tsv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["run", "interval_ns", "dcd_files", "sampling_ns", "sampled_frames", "bridge_frames", "any_bridge_occupancy_percent", "atom_level_events"])
    for row in overall:
        writer.writerow([row["run"], row["interval"], row["dcds"], row["sampling"], row["frames"], row["bridge_frames"], f'{row["occupancy"]:.6f}', row["events"]])

ranked = sorted(residues, key=lambda key: max(value["occupancy"] for value in residues[key].values()), reverse=True)
with (OUT / "leading_bridging_residues.tsv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["residue", "run1_frames", "run1_percent", "run2_frames", "run2_percent", "run3_frames", "run3_percent"])
    for key in ranked:
        values = [key]
        for run in ("run1", "run2", "run3"):
            entry = residues[key].get(run, {"frames": 0, "occupancy": 0.0})
            values.extend([entry["frames"], f'{entry["occupancy"]:.6f}'])
        writer.writerow(values)

with (OUT / "water_bridge_three_run_summary.md").open("w", encoding="utf-8") as out:
    out.write("# Z944 water-mediated bridge analysis: three-run summary\n\n")
    out.write("## Overall results\n\n")
    out.write("| Run | Interval (ns) | Sampled frames | Bridge frames | Any-bridge occupancy | Atom-level events |\n")
    out.write("|---|---:|---:|---:|---:|---:|\n")
    for row in overall:
        out.write(f'| {row["run"]} | {row["interval"]} | {row["frames"]:,} | {row["bridge_frames"]:,} | {row["occupancy"]:.2f}% | {row["events"]:,} |\n')
    out.write("\n## Leading bridging residues\n\n")
    out.write("| Residue | Run1 | Run2 | Run3 |\n|---|---:|---:|---:|\n")
    for key in ranked[:15]:
        vals = [residues[key].get(run, {"occupancy": 0.0})["occupancy"] for run in ("run1", "run2", "run3")]
        out.write(f"| {key} | {vals[0]:.2f}% | {vals[1]:.2f}% | {vals[2]:.2f}% |\n")

print(OUT)
