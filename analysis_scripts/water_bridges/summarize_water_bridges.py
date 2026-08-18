#!/usr/bin/env python3
"""Create dynamic per-DCD and time-binned summaries for a water-bridge run."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("result_dir", type=Path)
parser.add_argument("--first-dcd", type=int, required=True)
parser.add_argument("--last-dcd", type=int, required=True)
parser.add_argument("--start-ns", type=float, required=True)
parser.add_argument("--dcd-ns", type=float, default=20.0)
parser.add_argument("--frame-dt-ns", type=float, default=0.02)
parser.add_argument("--stride", type=int, default=10)
args = parser.parse_args()

events_path = args.result_dir / "bridge_events.tsv"
events = []
with events_path.open(encoding="utf-8") as source:
    for row in csv.DictReader(source, delimiter="\t"):
        events.append(row)

# Unique bridge observations per frame, residue, and water.
frame_events = defaultdict(int)
frame_residues = defaultdict(set)
frame_waters = defaultdict(set)
for row in events:
    frame = (int(row["dcd"]), int(row["dcd_frame"]))
    residue = (row["protein_segname"], row["protein_resid"], row["protein_resname"])
    water = (row["water_index"], row["water_segname"], row["water_resid"])
    frame_events[frame] += 1
    frame_residues[frame].add(residue)
    frame_waters[frame].add(water)

# Infer all sampled frames from the last observed/readable DCD frame. The Tcl
# summary supplies the authoritative total, including frames without bridges.
summary_text = (args.result_dir / "run_summary.txt").read_text(encoding="utf-8")
total_frames = None
for line in summary_text.splitlines():
    if line.startswith("Sampled frames:"):
        total_frames = int(line.split(":", 1)[1])
if total_frames is None:
    raise SystemExit("Sampled-frame count missing from run_summary.txt")

full_per_dcd = round(args.dcd_ns / (args.frame_dt_ns * args.stride))
remaining = total_frames
all_frames = []
per_dcd_total = {}
for dcd in range(args.first_dcd, args.last_dcd + 1):
    count = min(full_per_dcd, remaining)
    per_dcd_total[dcd] = count
    for sampled in range(count):
        original = sampled * args.stride
        time = args.start_ns + (dcd - args.first_dcd) * args.dcd_ns + original * args.frame_dt_ns
        all_frames.append((dcd, original, time))
    remaining -= count
if remaining != 0:
    raise SystemExit(f"Could not distribute {total_frames} frames across DCD files")

with (args.result_dir / "bridge_frame_timeseries.tsv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["time_ns", "dcd", "dcd_frame", "has_bridge", "event_count", "residue_count", "water_count"])
    for dcd, original, time in all_frames:
        key = (dcd, original)
        has_bridge = key in frame_events
        writer.writerow([
            f"{time:.3f}", dcd, original, int(has_bridge), frame_events.get(key, 0),
            len(frame_residues.get(key, set())), len(frame_waters.get(key, set())),
        ])

with (args.result_dir / "bridge_by_dcd.tsv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["dcd", "sampled_frames", "bridge_frames", "occupancy_percent"])
    for dcd in range(args.first_dcd, args.last_dcd + 1):
        bridge = sum(1 for key in frame_events if key[0] == dcd)
        total = per_dcd_total[dcd]
        writer.writerow([dcd, total, bridge, f"{100.0 * bridge / total:.6f}"])

# Ten-ns bins facilitate replica comparison even with a partial final DCD.
bins = defaultdict(lambda: [0, 0])
for dcd, original, time in all_frames:
    start = int(time // 10) * 10
    bins[start][0] += 1
    bins[start][1] += int((dcd, original) in frame_events)
with (args.result_dir / "bridge_10ns_occupancy.tsv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["bin_start_ns", "bin_end_ns", "sampled_frames", "bridge_frames", "occupancy_percent"])
    for start in sorted(bins):
        total, bridge = bins[start]
        writer.writerow([start, start + 10, total, bridge, f"{100.0 * bridge / total:.6f}"])

print(f"summarized {total_frames} frames and {len(events)} atom-level bridge events")
