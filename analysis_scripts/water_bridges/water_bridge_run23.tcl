# Adapter for the validated Z944 run1 water-bridge analysis.
# Run2/run3 store frames every 0.02 ns and each full DCD spans 20 ns.

set source_script "/work1/ted/June2023/WithLigand/md_run1/analyze/water-mediated/water_bridge_analysis.tcl"
set fh [open $source_script r]
set analysis [read $fh]
close $fh

set analysis [string map [list \
    {set frame_dt_ns 0.01} \
    {set frame_dt_ns 0.02} \
    {($dcdnum-$firstdcd)*10.0 + ($origframe+1)*$frame_dt_ns} \
    {($dcdnum-$firstdcd)*20.0 + $origframe*$frame_dt_ns}] $analysis]

eval $analysis
