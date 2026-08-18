set terminal pngcairo size 3300,1800 enhanced font "Sans,30"
set output outputfile
set title plot_title font ",34"
set xrange [xmin:xmax]
set grid back lc rgb "#d0d0d0"
set key outside right center spacing 1.15
set xlabel "Time (ns)"
set ylabel "Backbone RMSD (Å)"
plot datafile using 1:4 with lines lw 2 title "Transmembrane domain", \
     "" using 1:5 with lines lw 2 title "S6", \
     "" using 1:6 with lines lw 2 title "Selectivity filter", \
     "" using 1:7 with lines lw 2 title "Binding pocket", \
     "" using 1:8 with lines lw 2 title "Protein"
