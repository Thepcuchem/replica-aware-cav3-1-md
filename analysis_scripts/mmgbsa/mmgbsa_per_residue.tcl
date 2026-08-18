# Full NAMD/GBIS residue-ligand energy decomposition.
# Required environment:
#   PERRES_ROOT, PERRES_PARAM_DIR, PERRES_NAMD
# Optional:
#   PERRES_MAX (0=all), PERRES_KEEP (0=delete intermediates)

foreach variable {PERRES_ROOT PERRES_PARAM_DIR PERRES_NAMD} {
    if {![info exists env($variable)]} {
        puts stderr "Missing environment variable: $variable"
        exit 2
    }
}
set root [file normalize $env(PERRES_ROOT)]
set param_dir [file normalize $env(PERRES_PARAM_DIR)]
set namd_exe $env(PERRES_NAMD)
set max_residues 0
set keep_intermediate 0
if {[info exists env(PERRES_MAX)]} { set max_residues $env(PERRES_MAX) }
if {[info exists env(PERRES_KEEP)]} { set keep_intermediate $env(PERRES_KEEP) }

set work [file join $root per_residue]
file mkdir $work
cd $work
set psf_all [file join $root cpx cpx.psf]
set dcd_all [file join $root cpx sumtraj.cpx.dcd]
set ligand_log [file join $root lig mmgbsa_lig.out]

proc write_conf {path psf pdb dcd prefix param_dir} {
    set handle [open $path w]
    puts $handle "set path2par {$param_dir}"
    puts $handle "paraTypeCharmm on"
    puts $handle "structure {$psf}"
    puts $handle "coordinates {$pdb}"
    foreach parameter {
        par_all36m_prot.prm par_all36_lipid.prm par_all36_carb.prm
        par_all36_cgenff.prm par_all36_na.prm
        toppar_water_ions_cufix.str dzr.prm
    } {
        puts $handle "parameters \[file join \$path2par $parameter\]"
    }
    puts $handle "outputname {$prefix}"
    puts $handle "numsteps 0"
    puts $handle "temperature 310"
    puts $handle "GBIS on"
    puts $handle "solventDielectric 78.5"
    puts $handle "ionConcentration 0.3"
    puts $handle "GBISBeta 0.8"
    puts $handle "alphaCutoff 15"
    puts $handle "sasa on"
    puts $handle "surfaceTension 0.00542"
    puts $handle "switching on"
    puts $handle "switchdist 15"
    puts $handle "cutoff 16"
    puts $handle "pairlistdist 18"
    puts $handle "exclude scaled1-4"
    puts $handle "timestep 2"
    puts $handle "nonbondedFreq 2"
    puts $handle "fullElectFrequency 4"
    puts $handle "set ts 0"
    puts $handle "coorfile open dcd {$dcd}"
    puts $handle {while {![coorfile read]} {
    firstTimestep $ts
    run 0
    incr ts
}}
    puts $handle "coorfile close"
    close $handle
}

proc parse_energy {path} {
    set handle [open $path r]
    set started 0
    set result [dict create]
    while {[gets $handle line] >= 0} {
        if {[string match "*Info: Finished startup at*" $line]} {
            set started 1
            continue
        }
        if {!$started || ![string match "ENERGY:*" $line]} { continue }
        set fields [regexp -inline -all {\S+} $line]
        if {[llength $fields] < 8} { continue }
        set frame [lindex $fields 1]
        # Match the total MM/GBSA convention: omit NAMD timestep 0.
        if {![string is integer -strict $frame] || $frame <= 0} { continue }
        dict set result $frame [list [lindex $fields 6] [lindex $fields 7]]
    }
    close $handle
    return $result
}

proc mean_sd {values} {
    set n [llength $values]
    if {$n == 0} { return [list nan nan] }
    set sum 0.0
    foreach value $values { set sum [expr {$sum + $value}] }
    set mean [expr {$sum / double($n)}]
    if {$n == 1} { return [list $mean 0.0] }
    set ss 0.0
    foreach value $values { set ss [expr {$ss + ($value-$mean)*($value-$mean)}] }
    return [list $mean [expr {sqrt($ss / double($n-1))}]]
}

proc safe_delete {paths} {
    foreach path $paths { catch {file delete -force $path} }
}

puts "Loading $psf_all and $dcd_all"
set molid [mol new $psf_all type psf waitfor all]
mol addfile $dcd_all type dcd waitfor all molid $molid
set nframes [molinfo $molid get numframes]
set ligand [atomselect $molid "resname DZR" frame 0]
set protein [atomselect $molid "protein" frame 0]
set residues [lsort -integer -unique [$protein get resid]]
if {$max_residues > 0 && $max_residues < [llength $residues]} {
    set residues [lrange $residues 0 [expr {$max_residues - 1}]]
}
puts "Trajectory frames: $nframes"
puts "Residues selected: [llength $residues]"

if {![file exists $ligand_log]} {
    puts stderr "Missing ligand energy log: $ligand_log"
    exit 3
}
set lig_energy [parse_energy $ligand_log]
if {[dict size $lig_energy] == 0} {
    puts stderr "No ligand energies parsed from $ligand_log"
    exit 4
}

set frame_out [open per_residue_energy.csv w]
puts $frame_out "segid,resid,resname,frame,time_ns,Ebind_VDW,Ebind_ELEC,Ebind_TOTAL"
set average_out [open per_residue_statistics.csv w]
puts $average_out "segid,resid,resname,frames,mean_VDW,sd_VDW,mean_ELEC,sd_ELEC,mean_TOTAL,sd_TOTAL"

set completed 0
foreach resid $residues {
    set residue [atomselect $molid "protein and resid $resid" frame 0]
    set pair [atomselect $molid "(protein and resid $resid) or resname DZR" frame 0]
    set segid [lindex [lsort -unique [$residue get segid]] 0]
    set resname [lindex [lsort -unique [$residue get resname]] 0]
    set tag [format "%s_%04d_%s" $segid $resid $resname]
    puts "Residue [expr {$completed + 1}]/[llength $residues]: $tag"

    set generated {}
    foreach item [list [list res $residue] [list pair $pair]] {
        lassign $item name selection
        set psf "${tag}_${name}.psf"
        set pdb "${tag}_${name}.pdb"
        set dcd "${tag}_${name}.dcd"
        set conf "${tag}_${name}.conf"
        set log "${tag}_${name}.log"
        $selection frame 0
        $selection writepsf $psf
        $selection writepdb $pdb
        animate write dcd $dcd beg 0 end [expr {$nframes - 1}] skip 1 \
            waitfor all sel $selection $molid
        write_conf $conf $psf $pdb $dcd "${tag}_${name}_out" $param_dir
        if {[catch {
            exec $namd_exe +p1 $conf > $log 2>@1
        } message]} {
            puts stderr "NAMD failed for $tag $name: $message"
            exit 5
        }
        set energy($name) [parse_energy $log]
        lappend generated $psf $pdb $dcd $conf $log \
            "${tag}_${name}_out.coor" "${tag}_${name}_out.vel" "${tag}_${name}_out.xsc"
    }

    set common_frames {}
    foreach frame [lsort -integer [dict keys $energy(pair)]] {
        if {[dict exists $energy(res) $frame] && [dict exists $lig_energy $frame]} {
            lappend common_frames $frame
        }
    }
    if {[llength $common_frames] == 0} {
        puts stderr "No common frames for $tag"
        exit 6
    }

    set vdw_values {}
    set elec_values {}
    set total_values {}
    foreach frame $common_frames {
        lassign [dict get $energy(pair) $frame] pair_elec pair_vdw
        lassign [dict get $energy(res) $frame] res_elec res_vdw
        lassign [dict get $lig_energy $frame] lig_elec lig_vdw
        set bind_vdw [expr {$pair_vdw - $res_vdw - $lig_vdw}]
        set bind_elec [expr {$pair_elec - $res_elec - $lig_elec}]
        set bind_total [expr {$bind_vdw + $bind_elec}]
        set time_ns [expr {200.0 + $frame * 0.5}]
        puts $frame_out "$segid,$resid,$resname,$frame,$time_ns,$bind_vdw,$bind_elec,$bind_total"
        lappend vdw_values $bind_vdw
        lappend elec_values $bind_elec
        lappend total_values $bind_total
    }
    lassign [mean_sd $vdw_values] mean_vdw sd_vdw
    lassign [mean_sd $elec_values] mean_elec sd_elec
    lassign [mean_sd $total_values] mean_total sd_total
    puts $average_out "$segid,$resid,$resname,[llength $common_frames],$mean_vdw,$sd_vdw,$mean_elec,$sd_elec,$mean_total,$sd_total"
    flush $frame_out
    flush $average_out

    if {!$keep_intermediate} { safe_delete $generated }
    $residue delete
    $pair delete
    incr completed
}

close $frame_out
close $average_out
set summary [open summary.txt w]
puts $summary "Full protein-ligand per-residue NAMD/GBIS decomposition"
puts $summary "Protein residues: $completed"
puts $summary "Prepared trajectory frames: $nframes"
puts $summary "Reported frames per residue: [dict size $lig_energy]"
puts $summary "DCD source: $dcd_all"
puts $summary "Parameter directory: $param_dir"
puts $summary "NAMD executable: $namd_exe"
puts $summary "Frame 0 excluded to match total MM/GBSA analysis"
close $summary
$ligand delete
$protein delete
mol delete $molid
exit
