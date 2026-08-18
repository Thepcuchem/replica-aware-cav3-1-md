load mZ944_vs_Z944_determinants.pdb, mapped
hide everything, mapped
show cartoon, mapped and protein
color gray80, mapped and protein
spectrum b, blue_white_red, mapped and protein, minimum=-100, maximum=100
select determinants, mapped and resid 354+918+925+926+956+1459+1460+1461+1462+1463+1464+1465+1776+1777+1778+1780+1781+1782+1816+1820
show sticks, determinants
set stick_radius, 0.18, determinants
show spheres, mapped and resname DZR
set sphere_scale, 0.3, mapped and resname DZR
color yellow, mapped and resname DZR
distance pair_1_926_1463, mapped and name CA and resid 926, mapped and name CA and resid 1463
set dash_width, 2.5, pair_1_926_1463
color red, pair_1_926_1463
distance pair_2_925_1463, mapped and name CA and resid 925, mapped and name CA and resid 1463
set dash_width, 2.5, pair_2_925_1463
color red, pair_2_925_1463
distance pair_3_926_1461, mapped and name CA and resid 926, mapped and name CA and resid 1461
set dash_width, 2.5, pair_3_926_1461
color red, pair_3_926_1461
distance pair_4_1778_1464, mapped and name CA and resid 1778, mapped and name CA and resid 1464
set dash_width, 2.5, pair_4_1778_1464
color blue, pair_4_1778_1464
distance pair_5_1777_1464, mapped and name CA and resid 1777, mapped and name CA and resid 1464
set dash_width, 2.5, pair_5_1777_1464
color blue, pair_5_1777_1464
distance pair_6_1776_1464, mapped and name CA and resid 1776, mapped and name CA and resid 1464
set dash_width, 2.5, pair_6_1776_1464
color blue, pair_6_1776_1464
distance pair_7_926_1459, mapped and name CA and resid 926, mapped and name CA and resid 1459
set dash_width, 2.5, pair_7_926_1459
color red, pair_7_926_1459
distance pair_8_1778_1465, mapped and name CA and resid 1778, mapped and name CA and resid 1465
set dash_width, 2.5, pair_8_1778_1465
color blue, pair_8_1778_1465
distance pair_9_925_1461, mapped and name CA and resid 925, mapped and name CA and resid 1461
set dash_width, 2.5, pair_9_925_1461
color red, pair_9_925_1461
distance pair_10_1820_1463, mapped and name CA and resid 1820, mapped and name CA and resid 1463
set dash_width, 2.5, pair_10_1820_1463
color blue, pair_10_1820_1463
bg_color white
set ray_opaque_background, off
orient determinants
set_name mapped, mZ944_vs_Z944
save mZ944_vs_Z944_determinants.pse
