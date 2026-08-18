# Calculate mass-weighted ligand-to-pocket and ligand-to-filter COM distances.
# Usage: vmd -dispdev text -e com_distance_prepared.tcl -args \
#        SYSTEM RUN PSF PDB DCD OUTPUT_CSV START_NS INTERVAL_NS

if {$argc == 8} {
    lassign $argv system run psf pdb dcd output start_ns interval_ns
} elseif {[info exists env(COM_SYSTEM)]} {
    set system $env(COM_SYSTEM)
    set run $env(COM_RUN)
    set psf $env(COM_PSF)
    set pdb $env(COM_PDB)
    set dcd $env(COM_DCD)
    set output $env(COM_OUTPUT)
    set start_ns $env(COM_START_NS)
    set interval_ns $env(COM_INTERVAL_NS)
} else {
    puts stderr "Expected 8 arguments or COM_* environment variables"
    exit 2
}

set pocket_resids "384 387 388 391 868 872 875 876 916 917 918 920 921 922 948 950 951 952 953 955 956 957 959 960 1816 1819 1820 1823 1824 1462 1495 1498 1499 1502 1505 1506 1509 1510"
set filter_resids "351 to 357 919 to 926 1459 to 1466 1776 to 1782"

proc vector_distance {a b} {
    set dx [expr {[lindex $a 0] - [lindex $b 0]}]
    set dy [expr {[lindex $a 1] - [lindex $b 1]}]
    set dz [expr {[lindex $a 2] - [lindex $b 2]}]
    return [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
}

set molid [mol new $psf type psf waitfor all]
mol addfile $pdb type pdb waitfor all molid $molid
mol addfile $dcd type dcd waitfor all molid $molid

set ligand [atomselect $molid "resname DZR"]
set pocket [atomselect $molid "protein and backbone and resid $pocket_resids"]
set filter [atomselect $molid "protein and backbone and resid $filter_resids"]

if {[$ligand num] == 0 || [$pocket num] == 0 || [$filter num] == 0} {
    puts stderr "Empty selection: ligand=[$ligand num], pocket=[$pocket num], filter=[$filter num]"
    exit 3
}

set handle [open $output w]
puts $handle "system,run,frame,time_ns,pocket_distance_A,filter_distance_A"
set frames [molinfo $molid get numframes]
# Frame 0 is the coordinate PDB loaded before the DCD; analyze DCD frames only.
for {set frame 1} {$frame < $frames} {incr frame} {
    $ligand frame $frame
    $pocket frame $frame
    $filter frame $frame
    set ligand_com [measure center $ligand weight mass]
    set pocket_com [measure center $pocket weight mass]
    set filter_com [measure center $filter weight mass]
    set pocket_distance [vector_distance $ligand_com $pocket_com]
    set filter_distance [vector_distance $ligand_com $filter_com]
    set dcd_frame [expr {$frame - 1}]
    set time_ns [expr {$start_ns + $dcd_frame * $interval_ns}]
    puts $handle [format "%s,%d,%d,%.3f,%.6f,%.6f" \
        $system $run $dcd_frame $time_ns $pocket_distance $filter_distance]
}
close $handle

puts "Completed $system run$run: [expr {$frames - 1}] DCD frames; ligand=[$ligand num], pocket=[$pocket num], filter=[$filter num]"
$ligand delete
$pocket delete
$filter delete
mol delete $molid
exit
