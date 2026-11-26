load mZ944_vs_Apo_determinants.pdb, mapped
hide everything, mapped
show cartoon, mapped and protein
color gray80, mapped and protein
spectrum b, blue_white_red, mapped and protein, minimum=-100, maximum=100
select determinants, mapped and resid 353+354+355+356+357+384+916+917+918+919+920+921+923+924+925+1460+1461+1462+1463+1465+1466+1780+1781+1782+1816
show sticks, determinants
set stick_radius, 0.18, determinants
show spheres, mapped and resname DZR
set sphere_scale, 0.3, mapped and resname DZR
color yellow, mapped and resname DZR
distance pair_1_357_1781, mapped and name CA and resid 357, mapped and name CA and resid 1781
set dash_width, 2.5, pair_1_357_1781
color red, pair_1_357_1781
distance pair_2_1781_1816, mapped and name CA and resid 1781, mapped and name CA and resid 1816
set dash_width, 2.5, pair_2_1781_1816
color blue, pair_2_1781_1816
distance pair_3_384_921, mapped and name CA and resid 384, mapped and name CA and resid 921
set dash_width, 2.5, pair_3_384_921
color blue, pair_3_384_921
distance pair_4_384_917, mapped and name CA and resid 384, mapped and name CA and resid 917
set dash_width, 2.5, pair_4_384_917
color blue, pair_4_384_917
distance pair_5_384_1465, mapped and name CA and resid 384, mapped and name CA and resid 1465
set dash_width, 2.5, pair_5_384_1465
color blue, pair_5_384_1465
distance pair_6_354_1461, mapped and name CA and resid 354, mapped and name CA and resid 1461
set dash_width, 2.5, pair_6_354_1461
color red, pair_6_354_1461
distance pair_7_355_924, mapped and name CA and resid 355, mapped and name CA and resid 924
set dash_width, 2.5, pair_7_355_924
color red, pair_7_355_924
distance pair_8_384_918, mapped and name CA and resid 384, mapped and name CA and resid 918
set dash_width, 2.5, pair_8_384_918
color blue, pair_8_384_918
distance pair_9_357_1462, mapped and name CA and resid 357, mapped and name CA and resid 1462
set dash_width, 2.5, pair_9_357_1462
color red, pair_9_357_1462
distance pair_10_357_1461, mapped and name CA and resid 357, mapped and name CA and resid 1461
set dash_width, 2.5, pair_10_357_1461
color red, pair_10_357_1461
bg_color white
set ray_opaque_background, off
orient determinants
set_name mapped, mZ944_vs_Apo
save mZ944_vs_Apo_determinants.pse
