load Z944_vs_Apo_determinants.pdb, mapped
hide everything, mapped
show cartoon, mapped and protein
color gray80, mapped and protein
spectrum b, blue_white_red, mapped and protein, minimum=-100, maximum=100
select determinants, mapped and resid 351+353+354+355+356+357+384+387+388+872+875+876+916+917+918+919+920+921+925+926+948+950+955+960+1465+1466+1781
show sticks, determinants
set stick_radius, 0.18, determinants
show spheres, mapped and resname DZR
set sphere_scale, 0.3, mapped and resname DZR
color yellow, mapped and resname DZR
distance pair_1_384_920, mapped and name CA and resid 384, mapped and name CA and resid 920
set dash_width, 2.5, pair_1_384_920
color blue, pair_1_384_920
distance pair_2_356_1781, mapped and name CA and resid 356, mapped and name CA and resid 1781
set dash_width, 2.5, pair_2_356_1781
color red, pair_2_356_1781
distance pair_3_384_917, mapped and name CA and resid 384, mapped and name CA and resid 917
set dash_width, 2.5, pair_3_384_917
color blue, pair_3_384_917
distance pair_4_384_876, mapped and name CA and resid 384, mapped and name CA and resid 876
set dash_width, 2.5, pair_4_384_876
color blue, pair_4_384_876
distance pair_5_384_916, mapped and name CA and resid 384, mapped and name CA and resid 916
set dash_width, 2.5, pair_5_384_916
color blue, pair_5_384_916
distance pair_6_384_872, mapped and name CA and resid 384, mapped and name CA and resid 872
set dash_width, 2.5, pair_6_384_872
color blue, pair_6_384_872
distance pair_7_384_875, mapped and name CA and resid 384, mapped and name CA and resid 875
set dash_width, 2.5, pair_7_384_875
color blue, pair_7_384_875
distance pair_8_384_921, mapped and name CA and resid 384, mapped and name CA and resid 921
set dash_width, 2.5, pair_8_384_921
color blue, pair_8_384_921
distance pair_9_354_926, mapped and name CA and resid 354, mapped and name CA and resid 926
set dash_width, 2.5, pair_9_354_926
color blue, pair_9_354_926
distance pair_10_351_356, mapped and name CA and resid 351, mapped and name CA and resid 356
set dash_width, 2.5, pair_10_351_356
color red, pair_10_351_356
bg_color white
set ray_opaque_background, off
orient determinants
set_name mapped, Z944_vs_Apo
save Z944_vs_Apo_determinants.pse
