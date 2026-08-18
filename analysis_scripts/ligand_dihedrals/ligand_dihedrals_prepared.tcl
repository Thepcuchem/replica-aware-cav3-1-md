# Calculate the seven DZR dihedrals used in manuscript Figure 5.
# Environment variables: DIH_SYSTEM, DIH_RUN, DIH_PSF, DIH_PDB, DIH_DCD,
# DIH_OUTPUT, DIH_START_NS, and DIH_INTERVAL_NS.

foreach variable {DIH_SYSTEM DIH_RUN DIH_PSF DIH_PDB DIH_DCD DIH_OUTPUT DIH_START_NS DIH_INTERVAL_NS} {
    if {![info exists env($variable)]} {
        puts stderr "Missing environment variable $variable"
        exit 2
    }
}

set system $env(DIH_SYSTEM)
set run $env(DIH_RUN)
set psf $env(DIH_PSF)
set pdb $env(DIH_PDB)
set dcd $env(DIH_DCD)
set output $env(DIH_OUTPUT)
set start_ns $env(DIH_START_NS)
set interval_ns $env(DIH_INTERVAL_NS)

set definitions {
    {C9 N3 C8 C7}
    {N3 C8 C7 N1}
    {C8 C7 N1 C4}
    {C3 C1 C6 N2}
    {C1 C6 N2 C10}
    {C6 N2 C10 C14}
    {N2 C10 C14 C16}
}

set molid [mol new $psf type psf waitfor all]
mol addfile $pdb type pdb waitfor all molid $molid
mol addfile $dcd type dcd waitfor all molid $molid

set dihedral_indices {}
foreach definition $definitions {
    set indices {}
    foreach atom_name $definition {
        set selection [atomselect $molid "resname DZR and name $atom_name"]
        if {[$selection num] != 1} {
            puts stderr "Expected one DZR atom named $atom_name, found [$selection num]"
            exit 3
        }
        lappend indices [lindex [$selection get index] 0]
        $selection delete
    }
    lappend dihedral_indices $indices
}

set handle [open $output w]
puts $handle "system,run,frame,time_ns,chi1_deg,chi2_deg,chi3_deg,chi4_deg,chi5_deg,chi6_deg,chi7_deg"
set frames [molinfo $molid get numframes]
# Frame 0 is the coordinate PDB; frames 1 onward are the 600 DCD snapshots.
for {set frame 1} {$frame < $frames} {incr frame} {
    set dcd_frame [expr {$frame - 1}]
    set time_ns [expr {$start_ns + $dcd_frame * $interval_ns}]
    set values {}
    foreach indices $dihedral_indices {
        lappend values [measure dihed $indices frame $frame]
    }
    puts $handle [format "%s,%d,%d,%.3f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f" \
        $system $run $dcd_frame $time_ns \
        [lindex $values 0] [lindex $values 1] [lindex $values 2] \
        [lindex $values 3] [lindex $values 4] [lindex $values 5] \
        [lindex $values 6]]
}
close $handle
puts "Completed $system run$run: [expr {$frames - 1}] DCD frames"
mol delete $molid
exit
