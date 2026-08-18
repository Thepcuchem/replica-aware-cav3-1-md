# Prepare stripped complex, protein, and ligand trajectories for NAMD MM/GBSA.
# Usage:
#   vmd -dispdev text -e this_script.tcl -args PSF PDB TRAJDIR FIRST LAST STRIDE OUTDIR

if {$argc == 7} {
    lassign $argv psf_file pdb_file traj_dir first_dcd last_dcd stride out_dir
} elseif {
    [info exists env(MMGBSA_PSF)] &&
    [info exists env(MMGBSA_PDB)] &&
    [info exists env(MMGBSA_TRAJDIR)] &&
    [info exists env(MMGBSA_FIRST)] &&
    [info exists env(MMGBSA_LAST)] &&
    [info exists env(MMGBSA_STRIDE)] &&
    [info exists env(MMGBSA_OUTDIR)]
} {
    set psf_file $env(MMGBSA_PSF)
    set pdb_file $env(MMGBSA_PDB)
    set traj_dir $env(MMGBSA_TRAJDIR)
    set first_dcd $env(MMGBSA_FIRST)
    set last_dcd $env(MMGBSA_LAST)
    set stride $env(MMGBSA_STRIDE)
    set out_dir $env(MMGBSA_OUTDIR)
} else {
    puts stderr "Provide 7 arguments or the MMGBSA_* environment variables"
    exit 2
}
file mkdir $out_dir
foreach component {cpx prot lig} {
    file mkdir [file join $out_dir $component]
}

puts "Loading topology: $psf_file"
set molid [mol new $psf_file type psf waitfor all]
mol addfile $pdb_file type pdb waitfor all molid $molid

for {set index $first_dcd} {$index <= $last_dcd} {incr index} {
    set dcd_file [file join $traj_dir [format "md-%02d.dcd" $index]]
    if {![file exists $dcd_file]} {
        puts stderr "Missing trajectory: $dcd_file"
        exit 3
    }
    puts "Loading $dcd_file with stride $stride"
    mol addfile $dcd_file type dcd first 0 last -1 step $stride waitfor all molid $molid
}

# Frame 0 is the input PDB; trajectory frames begin at frame 1.
set first_frame 1
set last_frame [expr {[molinfo $molid get numframes] - 1}]
set trajectory_frames $last_frame
if {$trajectory_frames < 1} {
    puts stderr "No trajectory frames were loaded"
    exit 4
}

set cpx [atomselect $molid "protein or resname DZR" frame $first_frame]
set prot [atomselect $molid "protein" frame $first_frame]
set lig [atomselect $molid "resname DZR" frame $first_frame]

if {[$lig num] == 0} {
    puts stderr "The DZR ligand selection is empty"
    exit 5
}

foreach item [list [list cpx $cpx] [list prot $prot] [list lig $lig]] {
    lassign $item name selection
    set component_dir [file join $out_dir $name]
    $selection writepdb [file join $component_dir "${name}.pdb"]
    $selection writepsf [file join $component_dir "${name}.psf"]
    animate write dcd [file join $component_dir "sumtraj.${name}.dcd"] \
        beg $first_frame end $last_frame skip 1 waitfor all \
        sel $selection $molid
    puts "$name: [$selection num] atoms; $trajectory_frames frames"
}

set info [open [file join $out_dir preparation_summary.txt] w]
puts $info "PSF: $psf_file"
puts $info "PDB: $pdb_file"
puts $info "Trajectory directory: $traj_dir"
puts $info "DCD files: md-[format %02d $first_dcd].dcd through md-[format %02d $last_dcd].dcd"
puts $info "Native trajectory interval: 0.020 ns"
puts $info "Input stride: $stride"
puts $info "Snapshot interval: [format %.3f [expr {$stride * 0.020}]] ns"
puts $info "Prepared snapshots: $trajectory_frames"
puts $info "Complex atoms: [$cpx num]"
puts $info "Protein atoms: [$prot num]"
puts $info "Ligand atoms: [$lig num]"
close $info

$cpx delete
$prot delete
$lig delete
mol delete $molid
exit
