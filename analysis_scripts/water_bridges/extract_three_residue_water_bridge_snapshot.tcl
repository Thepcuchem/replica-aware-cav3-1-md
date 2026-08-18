# Extract a representative simultaneous Ser383/Asn957/Lys1462 water-bridge
# snapshot from Z944 run3, md-19 frame 690 (273.8 ns).

set psf "/work1/ted/June2023/WithLigand/ionized.psf"
set dcd "/work1/ted/June2023/WithLigand/md_run3/md-19.dcd"
set outdir "/work1/ted/June2023/WithLigand/analysis_replicas/water_bridges/Z944/representative_snapshot_three_residues"
set frame 690
set water_oxygen_indices "74 47 145232"

file mkdir $outdir
mol new $psf type psf waitfor all
mol addfile $dcd type dcd first $frame last $frame waitfor all

set selection "(segname A and resid 383) or (segname B and resid 957) or (segname C and resid 1462) or (segname L and resname DZR) or (same residue as index $water_oxygen_indices)"
set snapshot [atomselect top $selection frame 0]
puts "Selected [$snapshot num] atoms"
puts "Residues: [lsort -unique [$snapshot get {segname resid resname}]]"
$snapshot writepdb "$outdir/Z944_run3_273.8ns_three_residue_water_bridges.pdb"
$snapshot writepsf "$outdir/Z944_run3_273.8ns_three_residue_water_bridges.psf"

set info [open "$outdir/Z944_run3_273.8ns_three_residue_water_bridges.txt" w]
puts $info "Representative simultaneous protein-water-DZR bridge snapshot"
puts $info "System: Z944 run3"
puts $info "Trajectory: md-19.dcd"
puts $info "Original DCD frame: 690"
puts $info "Simulation time: 273.8 ns"
puts $info "Hydrogen-bond criterion: donor-acceptor <= 3.5 A; D-H...A >= 120 degrees"
puts $info ""
puts $info "Detected bridge paths from bridge_events.tsv:"
puts $info "A:SER383 O -- water WS1:894 (oxygen index 74) -- DZR O1"
puts $info "B:ASN957 OD1 -- water WT3:3977 (oxygen index 145232) -- DZR N3"
puts $info "B:ASN957 OD1 -- water WT3:3977 (oxygen index 145232) -- DZR N1"
puts $info "C:LYS1462 NZ -- water WS1:495 (oxygen index 47) -- DZR N2"
puts $info ""
puts $info "Included atoms: [$snapshot num]"
puts $info "Selection: $selection"
close $info

$snapshot delete
quit
