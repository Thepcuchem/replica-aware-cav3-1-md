# Protein-water-DZR single-water bridge analysis for Cav3.1-z944.
# Run: vmd -dispdev text -e water_bridge_analysis.tcl

set psf "/work1/ted/June2023/WithLigand/ionized.psf"
set trajdir "/work1/ted/June2023/WithLigand/md_run1"
set outdir "/work1/ted/June2023/WithLigand/md_run1/analyze/water-mediated/results"
set firstdcd 59
set lastdcd 89
set stride 10
set frame_dt_ns 0.01
set production_start_ns 302.0
set dist_cutoff 3.5
set angle_cutoff 60.0
set ligand_seg "L"
if {[info exists env(WB_PSF)]} { set psf $env(WB_PSF) }
if {[info exists env(WB_TRAJDIR)]} { set trajdir $env(WB_TRAJDIR) }
if {[info exists env(WB_OUTDIR)]} { set outdir $env(WB_OUTDIR) }
if {[info exists env(WB_FIRSTDCD)]} { set firstdcd $env(WB_FIRSTDCD) }
if {[info exists env(WB_LASTDCD)]} { set lastdcd $env(WB_LASTDCD) }
if {[info exists env(WB_STRIDE)]} { set stride $env(WB_STRIDE) }
if {[info exists env(WB_START_NS)]} { set production_start_ns $env(WB_START_NS) }
if {[info exists env(WB_LIGAND_SEG)]} { set ligand_seg $env(WB_LIGAND_SEG) }

file mkdir $outdir
set eventfh [open "$outdir/bridge_events.tsv" w]
puts $eventfh "time_ns\tdcd\tdcd_frame\twater_index\twater_segname\twater_resid\tprotein_segname\tprotein_resid\tprotein_resname\tprotein_atom\tligand_atom\tligand_water_direction\tprotein_water_direction"

mol new $psf type psf waitfor all
set mid [molinfo top]
set ligdon [atomselect $mid "segname $ligand_seg and resname DZR and name N2 N3"]
set ligacc [atomselect $mid "segname $ligand_seg and resname DZR and name N1 O1 O2"]
# PSF-only loading leaves the VMD element field as X for protein atoms, so do
# not filter on element here.  measure hbonds uses topology/bonding to identify
# eligible protein donors and acceptors.
set wat [atomselect $mid "water"]
puts "Ligand role selections: [$ligdon num] donors, [$ligacc num] acceptors"

# Build chemically explicit protein donor/acceptor selections. The PSF element
# field is unset, so CHARMM atom types, bonds, and residue chemistry are used.
# This is wrapped in a procedure to prevent VMD's console from echoing the
# large intermediate atom-index lists.
proc setup_protein_roles {mid} {
    global protdon protacc
    set protall [atomselect $mid "segname A B C D and protein"]
    set prot_indices [$protall get index]
    set prot_names [$protall get name]
    set prot_resnames [$protall get resname]
    set prot_types [$protall get type]
    set prot_bonds [$protall getbonds]
    set masssel [atomselect $mid all]
    set all_masses [$masssel get mass]
    $masssel delete
    set prot_donor_indices {}
    set prot_acceptor_indices {}
    foreach idx $prot_indices nm $prot_names rn $prot_resnames tp $prot_types bonds $prot_bonds {
        set first [string index $tp 0]
        if {$first ni {N O S}} { continue }
        set has_h 0
        foreach b $bonds {
            if {[lindex $all_masses $b] < 2.0} { set has_h 1; break }
        }
        if {$has_h} { lappend prot_donor_indices $idx }
        if {$first eq "O" || $first eq "S"} {
            lappend prot_acceptor_indices $idx
        } elseif {$first eq "N" && $rn in {HIS HSD HSE HSP} && $nm in {ND1 NE2} && !$has_h} {
            lappend prot_acceptor_indices $idx
        }
    }
    $protall delete
    set protdon [atomselect $mid "index [join $prot_donor_indices { }]"]
    set protacc [atomselect $mid "index [join $prot_acceptor_indices { }]"]
    $protdon global
    $protacc global
    puts "Protein role selections: [$protdon num] donors, [$protacc num] acceptors"
    return
}
setup_protein_roles $mid

# Static atom metadata, indexed by VMD's zero-based atom index.
array set aname {}
array set aseg {}
array set aresid {}
array set aresname {}
set allmeta [atomselect $mid "segname A B C D $ligand_seg or water"]
foreach idx [$allmeta get index] nm [$allmeta get name] sg [$allmeta get segname] ri [$allmeta get resid] rn [$allmeta get resname] {
    set aname($idx) $nm
    set aseg($idx) $sg
    set aresid($idx) $ri
    set aresname($idx) $rn
}
$allmeta delete

array set residue_frames {}
array set residue_events {}
array set atompair_frames {}
array set atompair_events {}
array set water_frames {}
array set motif_seen {}
set total_frames 0
set bridge_frames 0

proc add_hbonds {hb direction mapname} {
    upvar 1 $mapname hmap
    lassign $hb donors acceptors hydrogens
    foreach d $donors a $acceptors h $hydrogens {
        # The water atom is whichever endpoint belongs to a water residue.
        upvar 1 aresname aresname aname aname
        if {[info exists aresname($d)] && ($aresname($d) eq "TIP3" || $aresname($d) eq "HOH" || $aresname($d) eq "WAT")} {
            set w $d
            set other $a
        } else {
            set w $a
            set other $d
        }
        lappend hmap($w) [list $other $direction]
    }
}

for {set dcdnum $firstdcd} {$dcdnum <= $lastdcd} {incr dcdnum} {
    set dcd [format "%s/md-%02d.dcd" $trajdir $dcdnum]
    puts "Loading $dcd (stride $stride)"
    mol addfile $dcd type dcd first 0 last -1 step $stride waitfor all $mid
    set nf [molinfo $mid get numframes]
    puts "Processing $nf sampled frames from md-[format %02d $dcdnum]"

    for {set f 0} {$f < $nf} {incr f} {
        incr total_frames
        foreach s [list $ligdon $ligacc $protdon $protacc $wat] {$s frame $f}
        array unset lw
        array unset pw
        array set lw {}
        array set pw {}

        # Both H-bond directions are evaluated. VMD's angle cutoff is the
        # deviation from linearity: 60 degrees corresponds to D-H...A >= 120.
        add_hbonds [measure hbonds $dist_cutoff $angle_cutoff $ligdon $wat] "ligand_donor" lw
        add_hbonds [measure hbonds $dist_cutoff $angle_cutoff $wat $ligacc] "ligand_acceptor" lw
        if {[array size lw] == 0} { continue }

        # Limit the protein-water calculation to waters already H-bonded to DZR.
        set windices [array names lw]
        set bridgewat [atomselect $mid "water and index [join $windices { }]" frame $f]
        add_hbonds [measure hbonds $dist_cutoff $angle_cutoff $protdon $bridgewat] "protein_donor" pw
        add_hbonds [measure hbonds $dist_cutoff $angle_cutoff $bridgewat $protacc] "protein_acceptor" pw
        $bridgewat delete
        if {[array size pw] == 0} { continue }

        set frame_has_bridge 0
        array unset frame_res_seen
        array unset frame_pair_seen
        array unset frame_water_seen
        array set frame_res_seen {}
        array set frame_pair_seen {}
        array set frame_water_seen {}

        set origframe [expr {$f * $stride}]
        set time_ns [expr {$production_start_ns + ($dcdnum-$firstdcd)*10.0 + ($origframe+1)*$frame_dt_ns}]
        foreach w [array names lw] {
            if {![info exists pw($w)]} { continue }
            set frame_has_bridge 1
            set frame_water_seen($w) 1
            foreach le $lw($w) {
                lassign $le latom ldir
                foreach pe $pw($w) {
                    lassign $pe patom pdir
                    set rkey "$aseg($patom):$aresid($patom):$aresname($patom)"
                    set pkey "$rkey:$aname($patom)-$aname($latom)"
                    set frame_res_seen($rkey) 1
                    set frame_pair_seen($pkey) 1
                    incr residue_events($rkey)
                    incr atompair_events($pkey)
                    puts $eventfh [join [list [format %.3f $time_ns] $dcdnum $origframe $w $aseg($w) $aresid($w) $aseg($patom) $aresid($patom) $aresname($patom) $aname($patom) $aname($latom) $ldir $pdir] "\t"]
                }
            }
        }
        if {$frame_has_bridge} { incr bridge_frames }
        foreach k [array names frame_res_seen] { incr residue_frames($k) }
        foreach k [array names frame_pair_seen] { incr atompair_frames($k) }
        foreach k [array names frame_water_seen] { incr water_frames($k) }
    }
    animate delete all $mid
    flush $eventfh
    puts "Completed md-[format %02d $dcdnum]; total sampled frames=$total_frames, bridge frames=$bridge_frames"
}
close $eventfh

set sfh [open "$outdir/residue_bridge_occupancy.tsv" w]
puts $sfh "protein_segname\tprotein_resid\tprotein_resname\tframes\toccupancy_percent\tevent_count"
foreach key [lsort [array names residue_frames]] {
    lassign [split $key :] sg ri rn
    puts $sfh [join [list $sg $ri $rn $residue_frames($key) [format %.6f [expr {100.0*$residue_frames($key)/$total_frames}]] $residue_events($key)] "\t"]
}
close $sfh

set pfh [open "$outdir/atom_pair_bridge_occupancy.tsv" w]
puts $pfh "protein_segname\tprotein_resid\tprotein_resname\tprotein_atom\tligand_atom\tframes\toccupancy_percent\tevent_count"
foreach key [lsort [array names atompair_frames]] {
    regexp {^([^:]+):([^:]+):([^:]+):([^-]+)-(.+)$} $key -> sg ri rn pa la
    puts $pfh [join [list $sg $ri $rn $pa $la $atompair_frames($key) [format %.6f [expr {100.0*$atompair_frames($key)/$total_frames}]] $atompair_events($key)] "\t"]
}
close $pfh

set wfh [open "$outdir/water_participation.tsv" w]
puts $wfh "water_index\twater_segname\twater_resid\tframes\toccupancy_percent"
foreach w [lsort -integer [array names water_frames]] {
    puts $wfh [join [list $w $aseg($w) $aresid($w) $water_frames($w) [format %.6f [expr {100.0*$water_frames($w)/$total_frames}]]] "\t"]
}
close $wfh

set mfh [open "$outdir/run_summary.txt" w]
puts $mfh "Cav3.1-z944 protein-water-DZR bridge analysis"
puts $mfh "DCD files: md-[format %02d $firstdcd].dcd through md-[format %02d $lastdcd].dcd"
puts $mfh "Production-time mapping: $production_start_ns ns + trajectory time"
puts $mfh "Native frame interval: $frame_dt_ns ns; sampling stride: $stride ($stride*$frame_dt_ns ns)"
puts $mfh "Hydrogen bond criterion: donor-acceptor <= $dist_cutoff A; D-H...A >= [expr {180.0-$angle_cutoff}] degrees"
puts $mfh "Single-water bridge definition: same water H-bonded to DZR and a protein residue in the same sampled frame"
puts $mfh "Ligand donors: N2, N3; ligand acceptors: N1, O1, O2"
puts $mfh "Sampled frames: $total_frames"
puts $mfh "Frames with >=1 bridge: $bridge_frames"
puts $mfh "Any-bridge occupancy (%): [format %.6f [expr {100.0*$bridge_frames/$total_frames}]]"
close $mfh

puts "DONE: sampled frames=$total_frames bridge frames=$bridge_frames"
quit
