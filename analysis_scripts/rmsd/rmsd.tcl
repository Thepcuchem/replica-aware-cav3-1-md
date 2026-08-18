#===========================================================
# RMSD analysis
# Cav3.1 channel pore domain + Ca2+ ions
#===========================================================

set outfile [open rmsd.dat w]

puts $outfile "#Time(ns) DCD Frame PD_RMSD(A) S6_RMSD(A) SF_RMSD(A) POCKET_RMSD(A) PROTEIN_RMSD(A) LOOP_RMSD(A) LIG_RMSD(A)  Ion1_RMSD(A) Ion2_RMSD(A)"

set prefix md

#-----------------------------------------------------------
# Selections
#-----------------------------------------------------------

set s5 "resid 210 to 233 864 to 886 1389 to 1410 1716 to 1742"
set s6 "resid 370 to 397 939 to 968 1489 to 1516 1802 to 1832"
set s6e "resid 340 to 397 908 to 968 1448 to 1516 1764 to 1832"
set sf "resid 351 to 357 919 to 926 1459 to 1466 1776 to 1782"
# pocket "(noh and protein) and same residue as within 7.0 of (resname DZR and noh)"
# residue in the pocket "384 387 388 391 868 872 875 876 916 917 918 920 921 922 948 950 951 952 953 955 956 957 959 960 1816 1819 1820 1823 1824 1462 1495 1498 1499 1502 1505 1506 1509 1510"
set pk_resid "384 387 388 391 868 872 875 876 916 917 918 920 921 922 948 950 951 952 953 955 956 957 959 960 1816 1819 1820 1823 1824 1462 1495 1498 1499 1502 1505 1506 1509 1510"
set pocket "resid $pk_resid "

set seltext "protein and alpha"
set sellig "noh and resname DZR"

set selion1 "resname CAL and resid 2301"
set selion2 "resname CAL and resid 2302"

#-----------------------------------------------------------
# Global simulation time
#-----------------------------------------------------------

set time_ps 0.0
set stride  10

#-----------------------------------------------------------
# Loop over DCD files
#-----------------------------------------------------------

for {set j 58} {$j <= 89} {incr j} {

    puts "Processing $prefix-$j.dcd"

    mol load psf ../../../ionized.psf pdb ../../../ionized.pdb

    if {$j < 10} {
        mol addfile ../../$prefix-0$j.dcd \
            type dcd first 0 last -1 step $stride waitfor all
    } else {
        mol addfile ../../$prefix-$j.dcd \
            type dcd first 0 last -1 step $stride waitfor all
    }

    set nf [molinfo top get numframes]

    #-------------------------------------------------------
    # Reference structures (frame 0)
    #-------------------------------------------------------

    set frame0_pd  [atomselect top "$seltext and ($s5 or $s6e)" frame 0]
    set frame0_sf  [atomselect top "$seltext and ($sf)" frame 0]
    set frame0_s6  [atomselect top "$seltext and ($s6)" frame 0]
    set frame0_pk  [atomselect top "$seltext and $pocket" frame 0]
    set frame0_lig [atomselect top "$sellig" frame 0]
    set frame0_prot  [atomselect top "$seltext" frame 0]
    set frame0_loop  [atomselect top "$seltext and not ($s5 or $s6e)" frame 0]

    set framei1 [atomselect top "$selion1" frame 0]
    set framei2 [atomselect top "$selion2" frame 0]

    #-------------------------------------------------------
    # Mobile selections
    #-------------------------------------------------------

    set all   [atomselect top all]
    set sel_pd   [atomselect top "$seltext and ($s5 or $s6e)"]
    set sel_sf   [atomselect top "$seltext and ($sf)"]
    set sel_s6   [atomselect top "$seltext and ($s6)"]
    set sel_pk   [atomselect top "$seltext and $pocket"]
    set sel_lig  [atomselect top "$sellig"]
    set sel_prot  [atomselect top "$seltext"]
    set sel_loop  [atomselect top "$seltext and not ($s5 or $s6e)"]

    set seli1 [atomselect top "$selion1"]
    set seli2 [atomselect top "$selion2"]

    #-------------------------------------------------------
    # Frame spacing
    #-------------------------------------------------------

    if {$j <= 58} {
        set save_interval_ps 2.0
    } else {
        set save_interval_ps 10.0
    }

    set read_stride $stride

    set dt_ps [expr {$save_interval_ps * $read_stride}]

    #-------------------------------------------------------
    # RMSD calculation
    #-------------------------------------------------------

    for {set i 1} {$i < $nf} {incr i} {

        $all frame $i
        $sel_pd frame $i
        $sel_sf frame $i
        $sel_s6 frame $i
        $sel_pk frame $i
        $sel_lig frame $i
        $sel_prot frame $i
        $sel_loop frame $i
        $seli1 frame $i
        $seli2 frame $i

        # Align protein
        $all move [measure fit $sel_pd $frame0_pd]

        # Update simulation time
        set time_ps [expr {$time_ps + $dt_ps}]
        set time_ns [expr {$time_ps / 1000.0}]

        # RMSD calculations
        set prot_pd [measure rmsd $sel_pd   $frame0_pd]
        set prot_sf [measure rmsd $sel_sf   $frame0_sf]
        set prot_s6 [measure rmsd $sel_s6   $frame0_s6]
        set prot_pk [measure rmsd $sel_pk   $frame0_pk]
        set prot    [measure rmsd $sel_prot   $frame0_prot]
        set prot_loop [measure rmsd $sel_loop   $frame0_loop]
        set lig [measure rmsd $sel_lig $frame0_lig]
        set ion1 [measure rmsd $seli1 $framei1]
        set ion2 [measure rmsd $seli2 $framei2]

        puts $outfile \
            [format "%12.3f %4d %8d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f" \
                $time_ns $j $i $prot_pd $prot_sf $prot_s6 $prot_pk $prot $prot_loop $lig $ion1 $ion2]
    }

    #-------------------------------------------------------
    # Cleanup
    #-------------------------------------------------------

    $frame0_pd delete
    $frame0_sf delete
    $frame0_s6 delete
    $framei1 delete
    $framei2 delete

    $all delete
    $sel_pd delete
    $sel_sf delete
    $sel_s6 delete
    $seli1 delete
    $seli2 delete

    mol delete all
}

close $outfile

exit
