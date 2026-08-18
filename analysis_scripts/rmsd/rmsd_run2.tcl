# RMSD analysis for Z944 run2, final completed 300--640 ns interval.
# Run with: vmd -dispdev text -e rmsd_run2.tcl

set psf "/work1/ted/June2023/WithLigand/ionized.psf"
set trajdir "/work1/ted/June2023/WithLigand/md_run2"
set output_file "rmsd_run2_300_640ns.dat"
set include_ligand 1
set variable_run1_timing 0
if {[info exists env(RMSD_PSF)]} { set psf $env(RMSD_PSF) }
if {[info exists env(RMSD_TRAJDIR)]} { set trajdir $env(RMSD_TRAJDIR) }
if {[info exists env(RMSD_OUTPUT)]} { set output_file $env(RMSD_OUTPUT) }
if {[info exists env(RMSD_INCLUDE_LIGAND)]} { set include_ligand $env(RMSD_INCLUDE_LIGAND) }
if {[info exists env(RMSD_VARIABLE_RUN1_TIMING)]} { set variable_run1_timing $env(RMSD_VARIABLE_RUN1_TIMING) }
set outfile [open $output_file w]

set firstdcd 13
set lastdcd 29
set stride 10
set frame_dt_ns 0.02
set start_time_ns 300.0
if {[info exists env(RMSD_FIRSTDCD)]} { set firstdcd $env(RMSD_FIRSTDCD) }
if {[info exists env(RMSD_LASTDCD)]} { set lastdcd $env(RMSD_LASTDCD) }
if {[info exists env(RMSD_START_NS)]} { set start_time_ns $env(RMSD_START_NS) }

set s5 "resid 210 to 233 864 to 886 1389 to 1410 1716 to 1742"
set s6 "resid 370 to 397 939 to 968 1489 to 1516 1802 to 1832"
set s6e "resid 340 to 397 908 to 968 1448 to 1516 1764 to 1832"
set sf "resid 351 to 357 919 to 926 1459 to 1466 1776 to 1782"
set pk_resid "384 387 388 391 868 872 875 876 916 917 918 920 921 922 948 950 951 952 953 955 956 957 959 960 1816 1819 1820 1823 1824 1462 1495 1498 1499 1502 1505 1506 1509 1510"

set aligntext "protein and alpha and ($s5 or $s6e)"
set pdtext "protein and backbone and ($s5 or $s6e)"
set s6text "protein and backbone and ($s6)"
set sftext "protein and backbone and ($sf)"
set pockettext "protein and backbone and resid $pk_resid"
set proteintext "protein and backbone"
set ligtext "noh and resname DZR"

if {$include_ligand} {
    puts $outfile "# Time_ns DCD Frame PD_RMSD_A S6_RMSD_A SF_RMSD_A Pocket_RMSD_A Protein_RMSD_A DZR_RMSD_A"
} else {
    puts $outfile "# Time_ns DCD Frame PD_RMSD_A S6_RMSD_A SF_RMSD_A Pocket_RMSD_A Protein_RMSD_A"
}

# Use one common reference: frame 0 of the first selected DCD.
set refmol [mol new $psf type psf waitfor all]
set reffile [format "%s/md-%02d.dcd" $trajdir $firstdcd]
mol addfile $reffile type dcd first 0 last 0 waitfor all molid $refmol
set ref_align [atomselect $refmol $aligntext frame 0]
set ref_pd [atomselect $refmol $pdtext frame 0]
set ref_s6 [atomselect $refmol $s6text frame 0]
set ref_sf [atomselect $refmol $sftext frame 0]
set ref_pocket [atomselect $refmol $pockettext frame 0]
set ref_protein [atomselect $refmol $proteintext frame 0]
if {$include_ligand} { set ref_lig [atomselect $refmol $ligtext frame 0] }

set global_frame 0
set current_time_ns $start_time_ns
for {set j $firstdcd} {$j <= $lastdcd} {incr j} {
    set dcdfile [format "%s/md-%02d.dcd" $trajdir $j]
    puts "Processing $dcdfile"

    set current_stride $stride
    set current_frame_dt_ns $frame_dt_ns
    if {$variable_run1_timing} {
        if {$variable_run1_timing == 3 && $j <= 61} {
            # Apo run1 md-03--61: 2 ps per stored frame.
            set current_stride 50
            set current_frame_dt_ns 0.002
        } elseif {$variable_run1_timing == 3 && $j <= 89} {
            # Apo run1 md-62--89: 10 ps per stored frame.
            set current_stride 10
            set current_frame_dt_ns 0.010
        } elseif {$variable_run1_timing == 3} {
            # Apo run1 md-90 onward: 100 ps per stored frame.
            set current_stride 1
            set current_frame_dt_ns 0.100
        } elseif {$variable_run1_timing == 2 && $j <= 5} {
            # mZ944 run1 early files: 2 ps per stored frame.
            set current_stride 50
            set current_frame_dt_ns 0.002
        } elseif {$variable_run1_timing == 1 && $j <= 48} {
            # Z944 run1 early files: 5 ps per stored frame.
            set current_stride 20
            set current_frame_dt_ns 0.005
        } else {
            set current_stride 10
            set current_frame_dt_ns 0.010
        }
    }

    set mol [mol new $psf type psf waitfor all]
    mol addfile $dcdfile type dcd first 0 last -1 step $current_stride waitfor all molid $mol
    set nf [molinfo $mol get numframes]

    # Apo run1 absolute times derived from DCD ISTART values (2 fs/step).
    # These branches preserve interruptions and split restart files.
    if {$variable_run1_timing == 3} {
        if {$j <= 61} {
            set current_time_ns [expr {($j - 3) * 1.0}]
        } elseif {$j <= 82} {
            set current_time_ns [expr {59.006 + ($j - 62) * 10.0}]
        } elseif {$j == 83} {
            set current_time_ns 269.006
        } elseif {$j == 84} {
            set current_time_ns 273.276
        } elseif {$j <= 89} {
            set current_time_ns [expr {273.286 + ($j - 85) * 10.0}]
        } else {
            set current_time_ns [expr {323.376 + ($j - 90) * 20.0}]
        }
        if {$j == 73} { puts $outfile "" }
    }

    set all [atomselect $mol all]
    set mob_align [atomselect $mol $aligntext]
    set mob_pd [atomselect $mol $pdtext]
    set mob_s6 [atomselect $mol $s6text]
    set mob_sf [atomselect $mol $sftext]
    set mob_pocket [atomselect $mol $pockettext]
    set mob_protein [atomselect $mol $proteintext]
    if {$include_ligand} { set mob_lig [atomselect $mol $ligtext] }

    for {set i 0} {$i < $nf} {incr i} {
        set mobile_sels [list $all $mob_align $mob_pd $mob_s6 $mob_sf $mob_pocket $mob_protein]
        if {$include_ligand} { lappend mobile_sels $mob_lig }
        foreach sel $mobile_sels {
            $sel frame $i
        }
        $all move [measure fit $mob_align $ref_align]

        set time_ns $current_time_ns
        set values [list $time_ns $j [expr {$i * $current_stride}] \
            [measure rmsd $mob_pd $ref_pd] [measure rmsd $mob_s6 $ref_s6] \
            [measure rmsd $mob_sf $ref_sf] [measure rmsd $mob_pocket $ref_pocket] \
            [measure rmsd $mob_protein $ref_protein]]
        if {$include_ligand} {
            lappend values [measure rmsd $mob_lig $ref_lig]
            puts $outfile [format "%10.3f %3d %6d %10.4f %10.4f %10.4f %13.4f %14.4f %10.4f" {*}$values]
        } else {
            puts $outfile [format "%10.3f %3d %6d %10.4f %10.4f %10.4f %13.4f %14.4f" {*}$values]
        }
        incr global_frame
        set current_time_ns [expr {$current_time_ns + $current_frame_dt_ns * $current_stride}]
    }

    foreach sel $mobile_sels {
        $sel delete
    }
    mol delete $mol

    # Preserve exact continuation times when the source frame count is not an
    # integer multiple of the sampling stride (mZ944 run1 only).
    if {$variable_run1_timing == 2 && $j == 3} {
        set current_time_ns [expr {$current_time_ns - 0.096}]
    }
    if {$variable_run1_timing == 2 && $j == 31} {
        set current_time_ns [expr {$current_time_ns - 0.020}]
    }
}

set reference_sels [list $ref_align $ref_pd $ref_s6 $ref_sf $ref_pocket $ref_protein]
if {$include_ligand} { lappend reference_sels $ref_lig }
foreach sel $reference_sels {
    $sel delete
}
mol delete $refmol
close $outfile
puts "Wrote $output_file ($global_frame sampled frames)"
exit
