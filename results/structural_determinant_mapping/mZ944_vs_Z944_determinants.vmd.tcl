mol new mZ944_vs_Z944_determinants.pdb type pdb waitfor all
set molid [molinfo top]
mol delrep 0 $molid
color scale method BWR
mol representation NewCartoon
mol color Beta
mol selection {protein}
mol material Opaque
mol addrep $molid
mol scaleminmax $molid 0 -100 100
mol representation Licorice 0.18 12.0 12.0
mol color Beta
mol selection {protein and resid 354 918 925 926 956 1459 1460 1461 1462 1463 1464 1465 1776 1777 1778 1780 1781 1782 1816 1820}
mol material Opaque
mol addrep $molid
mol scaleminmax $molid 1 -100 100
mol representation Licorice 0.25 12.0 12.0
mol color ColorID 4
mol selection {resname DZR}
mol material Opaque
mol addrep $molid
graphics $molid materials on
set sel1_1 [atomselect $molid {protein and name CA and resid 926}]
set sel2_1 [atomselect $molid {protein and name CA and resid 1463}]
set p1_1 [lindex [$sel1_1 get {x y z}] 0]
set p2_1 [lindex [$sel2_1 get {x y z}] 0]
graphics $molid color red
graphics $molid line $p1_1 $p2_1 width 3 style dashed
$sel1_1 delete
$sel2_1 delete
set sel1_2 [atomselect $molid {protein and name CA and resid 925}]
set sel2_2 [atomselect $molid {protein and name CA and resid 1463}]
set p1_2 [lindex [$sel1_2 get {x y z}] 0]
set p2_2 [lindex [$sel2_2 get {x y z}] 0]
graphics $molid color red
graphics $molid line $p1_2 $p2_2 width 3 style dashed
$sel1_2 delete
$sel2_2 delete
set sel1_3 [atomselect $molid {protein and name CA and resid 926}]
set sel2_3 [atomselect $molid {protein and name CA and resid 1461}]
set p1_3 [lindex [$sel1_3 get {x y z}] 0]
set p2_3 [lindex [$sel2_3 get {x y z}] 0]
graphics $molid color red
graphics $molid line $p1_3 $p2_3 width 3 style dashed
$sel1_3 delete
$sel2_3 delete
set sel1_4 [atomselect $molid {protein and name CA and resid 1778}]
set sel2_4 [atomselect $molid {protein and name CA and resid 1464}]
set p1_4 [lindex [$sel1_4 get {x y z}] 0]
set p2_4 [lindex [$sel2_4 get {x y z}] 0]
graphics $molid color blue
graphics $molid line $p1_4 $p2_4 width 3 style dashed
$sel1_4 delete
$sel2_4 delete
set sel1_5 [atomselect $molid {protein and name CA and resid 1777}]
set sel2_5 [atomselect $molid {protein and name CA and resid 1464}]
set p1_5 [lindex [$sel1_5 get {x y z}] 0]
set p2_5 [lindex [$sel2_5 get {x y z}] 0]
graphics $molid color blue
graphics $molid line $p1_5 $p2_5 width 3 style dashed
$sel1_5 delete
$sel2_5 delete
set sel1_6 [atomselect $molid {protein and name CA and resid 1776}]
set sel2_6 [atomselect $molid {protein and name CA and resid 1464}]
set p1_6 [lindex [$sel1_6 get {x y z}] 0]
set p2_6 [lindex [$sel2_6 get {x y z}] 0]
graphics $molid color blue
graphics $molid line $p1_6 $p2_6 width 3 style dashed
$sel1_6 delete
$sel2_6 delete
set sel1_7 [atomselect $molid {protein and name CA and resid 926}]
set sel2_7 [atomselect $molid {protein and name CA and resid 1459}]
set p1_7 [lindex [$sel1_7 get {x y z}] 0]
set p2_7 [lindex [$sel2_7 get {x y z}] 0]
graphics $molid color red
graphics $molid line $p1_7 $p2_7 width 3 style dashed
$sel1_7 delete
$sel2_7 delete
set sel1_8 [atomselect $molid {protein and name CA and resid 1778}]
set sel2_8 [atomselect $molid {protein and name CA and resid 1465}]
set p1_8 [lindex [$sel1_8 get {x y z}] 0]
set p2_8 [lindex [$sel2_8 get {x y z}] 0]
graphics $molid color blue
graphics $molid line $p1_8 $p2_8 width 3 style dashed
$sel1_8 delete
$sel2_8 delete
set sel1_9 [atomselect $molid {protein and name CA and resid 925}]
set sel2_9 [atomselect $molid {protein and name CA and resid 1461}]
set p1_9 [lindex [$sel1_9 get {x y z}] 0]
set p2_9 [lindex [$sel2_9 get {x y z}] 0]
graphics $molid color red
graphics $molid line $p1_9 $p2_9 width 3 style dashed
$sel1_9 delete
$sel2_9 delete
set sel1_10 [atomselect $molid {protein and name CA and resid 1820}]
set sel2_10 [atomselect $molid {protein and name CA and resid 1463}]
set p1_10 [lindex [$sel1_10 get {x y z}] 0]
set p2_10 [lindex [$sel2_10 get {x y z}] 0]
graphics $molid color blue
graphics $molid line $p1_10 $p2_10 width 3 style dashed
$sel1_10 delete
$sel2_10 delete
display projection Orthographic
color Display Background white
display resetview
puts {Blue = reproducible contraction; red = reproducible expansion}
quit
