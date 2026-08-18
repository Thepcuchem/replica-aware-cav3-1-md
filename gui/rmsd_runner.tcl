# ReplicaLab generic RMSD adapter.
# Required environment variables:
# RL_PSF, RL_PDB, RL_MANIFEST, RL_OUTPUT, RL_ALIGNMENT, RL_SELECTIONS

foreach variable {RL_PSF RL_PDB RL_MANIFEST RL_OUTPUT RL_ALIGNMENT RL_SELECTIONS} {
    if {![info exists env($variable)]} {
        puts stderr "Missing environment variable $variable"
        exit 2
    }
}

set psf $env(RL_PSF)
set manifest_path $env(RL_MANIFEST)
set output_path $env(RL_OUTPUT)
set alignment_text $env(RL_ALIGNMENT)
set selection_spec $env(RL_SELECTIONS)

set manifest [open $manifest_path r]
gets $manifest header
set segments {}
while {[gets $manifest line] >= 0} {
    if {[string trim $line] eq ""} { continue }
    # Paths containing commas are not supported by this compact manifest format.
    set fields [split $line ","]
    lappend segments [list [lindex $fields 0] [lindex $fields 1] [lindex $fields 2] [lindex $fields 3]]
}
close $manifest

if {[llength $segments] == 0} {
    puts stderr "Trajectory manifest is empty"
    exit 3
}

set first_path [lindex [lindex $segments 0] 0]
set refmol [mol new $psf type psf waitfor all]
mol addfile $first_path type dcd first 0 last 0 step 1 waitfor all molid $refmol
set refalign [atomselect $refmol $alignment_text frame 0]
if {[$refalign num] == 0} {
    puts stderr "Alignment selection is empty: $alignment_text"
    exit 4
}

set definitions {}
foreach item [split $selection_spec ";"] {
    if {[string trim $item] eq ""} { continue }
    set separator [string first "|" $item]
    if {$separator < 1} { continue }
    set name [string range $item 0 [expr {$separator - 1}]]
    set text [string range $item [expr {$separator + 1}] end]
    lappend definitions [list $name $text]
}
if {[llength $definitions] == 0} {
    puts stderr "No RMSD measurement selections were supplied"
    exit 5
}

set output [open $output_path w]
puts -nonewline $output "time_ns,segment,frame"
foreach definition $definitions {
    puts -nonewline $output ",[lindex $definition 0]"
}
puts $output ""

set molid [mol new $psf type psf waitfor all]
set segment_number 0
foreach segment $segments {
    incr segment_number
    lassign $segment dcd_path start_ns interval_ns stride
    mol addfile $dcd_path type dcd step $stride waitfor all molid $molid
    set frames [molinfo $molid get numframes]
    set mobile [atomselect $molid $alignment_text]
    if {[$mobile num] != [$refalign num]} {
        puts stderr "Alignment atom count mismatch in $dcd_path"
        exit 6
    }
    set all [atomselect $molid all]
    set measurement_selections {}
    foreach definition $definitions {
        set selection [atomselect $molid [lindex $definition 1]]
        if {[$selection num] == 0} {
            puts stderr "Measurement selection is empty: [lindex $definition 1]"
            exit 7
        }
        lappend measurement_selections $selection
    }
    for {set frame 0} {$frame < $frames} {incr frame} {
        $mobile frame $frame
        $all frame $frame
        set transform [measure fit $mobile $refalign]
        $all move $transform
        set time_ns [expr {double($start_ns) + $frame * double($interval_ns) * $stride}]
        puts -nonewline $output [format "%.6f,%d,%d" $time_ns $segment_number $frame]
        foreach selection $measurement_selections {
            $selection frame $frame
            set reference [atomselect $refmol [$selection text] frame 0]
            if {[$reference num] != [$selection num]} {
                puts stderr "Reference atom count mismatch for [$selection text]"
                exit 8
            }
            set value [measure rmsd $selection $reference]
            puts -nonewline $output [format ",%.6f" $value]
            $reference delete
        }
        puts $output ""
    }
    foreach selection $measurement_selections { $selection delete }
    $mobile delete
    $all delete
    animate delete beg 0 end -1 molid $molid
}

close $output
$refalign delete
mol delete $refmol
mol delete $molid
puts "ReplicaLab RMSD completed: $output_path"
exit
