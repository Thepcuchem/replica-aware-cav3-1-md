#!/usr/bin/env python3
"""Create the run2 RMSD plot and summary using Python's standard library.

The plot itself is rendered by the system gnuplot executable, avoiding a
dependency on numpy/matplotlib.
"""
import math
import statistics
import subprocess
import argparse
from pathlib import Path


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Plot and summarize an RMSD table")
parser.add_argument("--run", choices=("run1", "run2", "run3"), default="run2")
parser.add_argument("--system", choices=("Z944", "mZ944", "apo"), default="Z944")
parser.add_argument("--window", choices=("300_640", "0_500", "50_650", "0_600"), default="300_640")
args = parser.parse_args()
prefix = {"Z944": "rmsd", "mZ944": "rmsd_mZ944", "apo": "rmsd_apo"}[args.system]
stem = f"{prefix}_{args.run}_{args.window}ns"
data_path = HERE / f"{stem}.dat"
plot_path = HERE / f"{stem}.png"
summary_path = HERE / f"{stem}_summary.txt"
gnuplot_path = HERE / "plot_rmsd_run2.gnuplot"

rows = []
with data_path.open(encoding="utf-8") as source:
    for line in source:
        if line.strip() and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])

expected_columns = 8 if args.system == "apo" else 9
if not rows or any(len(row) != expected_columns for row in rows):
    raise SystemExit(f"Invalid or empty RMSD table: {data_path}")
if any(not math.isfinite(value) for row in rows for value in row):
    raise SystemExit(f"Non-finite value found in {data_path}")

series = {
    "Transmembrane domain": [row[3] for row in rows],
    "S6": [row[4] for row in rows],
    "Selectivity filter": [row[5] for row in rows],
    "Binding pocket": [row[6] for row in rows],
    "Protein": [row[7] for row in rows],
}
if args.system != "apo":
    series["DZR"] = [row[8] for row in rows]

with summary_path.open("w", encoding="utf-8") as out:
    out.write(f"Samples: {len(rows)}\n")
    out.write(f"Time range: {rows[0][0]:.3f}--{rows[-1][0]:.3f} ns\n")
    out.write("RMSD statistics (angstrom):\n")
    for name, values in series.items():
        out.write(
            f"{name:24s} mean={statistics.fmean(values):8.3f}  "
            f"SD={statistics.stdev(values):8.3f}  "
            f"min={min(values):8.3f}  max={max(values):8.3f}\n"
        )

xmax = {"0_500": 500, "300_640": 640, "50_650": 650, "0_600": 600}[args.window]
settings = (
    f"datafile='{data_path.name}'; outputfile='{plot_path.name}'; "
    f"plot_title='{args.system} {args.run} RMSD (common reference at {rows[0][0]:g} ns)'; "
    f"xmin={rows[0][0]:g}; xmax={xmax}"
)
if args.system == "apo":
    gnuplot_path = HERE / "plot_rmsd_apo.gnuplot"
subprocess.run(["gnuplot", "-e", settings, str(gnuplot_path)], cwd=HERE, check=True)
print(plot_path)
print(summary_path)
