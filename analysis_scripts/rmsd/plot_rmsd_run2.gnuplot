set terminal pngcairo size 3300,2400 enhanced font "Sans,30"
if (!exists("datafile")) datafile = "rmsd_run2_300_640ns.dat"
if (!exists("outputfile")) outputfile = "rmsd_run2_300_640ns.png"
if (!exists("plot_title")) plot_title = "Z944 run2 RMSD (common reference at 300 ns)"
if (!exists("xmin")) xmin = 300
if (!exists("xmax")) xmax = 640
set output outputfile
set multiplot layout 2,1 title plot_title font ",34"
set xrange [xmin:xmax]
set grid back lc rgb "#d0d0d0"
set key outside right center spacing 1.15
set ylabel "Backbone RMSD (Å)"
set format x ""
plot datafile using 1:4 with lines lw 2 title "Transmembrane domain", \
     "" using 1:5 with lines lw 2 title "S6", \
     "" using 1:6 with lines lw 2 title "Selectivity filter", \
     "" using 1:7 with lines lw 2 title "Binding pocket", \
     "" using 1:8 with lines lw 2 title "Protein"
unset title
set format x "%g"
set xlabel "Time (ns)"
set ylabel "DZR heavy-atom RMSD (Å)"
unset key
plot datafile using 1:9 with lines lw 2 lc rgb "#9c2f45"
unset multiplot
