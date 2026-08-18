#!/usr/bin/env python3
"""Compare Z944 water-mediated bridge occupancies between two replicas."""
import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent / "water_bridges" / "Z944"
OUT = ROOT / "comparison"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with path.open(encoding="utf-8") as source:
        return list(csv.DictReader(source, delimiter="\t"))


occupancy = {}
for run in ("run2", "run3"):
    for row in read_tsv(ROOT / run / "residue_bridge_occupancy.tsv"):
        key = f'{row["protein_segname"]}:{row["protein_resid"]}:{row["protein_resname"]}'
        occupancy.setdefault(key, {})[run] = float(row["occupancy_percent"])

ranked = sorted(occupancy, key=lambda key: max(occupancy[key].values()), reverse=True)
with (OUT / "residue_occupancy_comparison.tsv").open("w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["residue", "run2_percent", "run3_percent"])
    for key in ranked:
        writer.writerow([key, f'{occupancy[key].get("run2", 0):.6f}', f'{occupancy[key].get("run3", 0):.6f}'])

with (OUT / "top12_residue_occupancy.dat").open("w", encoding="utf-8") as out:
    for key in ranked[:12]:
        out.write(f'{key} {occupancy[key].get("run2", 0):.6f} {occupancy[key].get("run3", 0):.6f}\n')

# Combine 10-ns time bins for a two-line comparison plot.
bins = {}
for run in ("run2", "run3"):
    for row in read_tsv(ROOT / run / "bridge_10ns_occupancy.tsv"):
        start = float(row["bin_start_ns"])
        bins.setdefault(start, {})[run] = float(row["occupancy_percent"])
with (OUT / "bridge_10ns_comparison.dat").open("w", encoding="utf-8") as out:
    for start in sorted(bins):
        out.write(f'{start + 5:g} {bins[start].get("run2", float("nan")):g} {bins[start].get("run3", float("nan")):g}\n')

gnuplot = OUT / "water_bridge_comparison.gnuplot"
gnuplot.write_text(
    '''set terminal pngcairo size 3300,2400 enhanced font "Sans,30"
set output "water_bridge_comparison.png"
set multiplot layout 2,1 title "Z944 protein-water-DZR bridges (200-500 ns)" font ",34"
set grid back lc rgb "#d0d0d0"
set key top right
set xrange [200:500]
set yrange [0:105]
set ylabel "Frames with bridge (%)"
set xlabel "Time (ns)"
plot "bridge_10ns_comparison.dat" using 1:2 with linespoints lw 3 pt 7 title "run2", \\
     "" using 1:3 with linespoints lw 3 pt 5 title "run3"
unset grid
set style data histograms
set style histogram clustered gap 1
set style fill solid 0.85 border -1
set boxwidth 0.9
set auto x
set yrange [0:*]
set ylabel "Residue bridge occupancy (%)"
set xlabel "Protein residue"
set xtics rotate by -35
plot "top12_residue_occupancy.dat" using 2:xtic(1) title "run2", \\
     "" using 3 title "run3"
unset multiplot
''',
    encoding="utf-8",
)
subprocess.run(["gnuplot", gnuplot.name], cwd=OUT, check=True)
print(OUT / "water_bridge_comparison.png")
