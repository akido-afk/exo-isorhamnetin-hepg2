reinitialize

bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on
set depth_cue, off
set specular, 0.15
set cartoon_fancy_helices, 1
set cartoon_sampling, 14


# ============================================================
# LOAD CLOSED 1Z83
# ============================================================

load AK1_1Z83.pdb, closed_1Z83

remove closed_1Z83 and not chain A

hide everything, closed_1Z83

show cartoon, closed_1Z83 and polymer.protein

color gray70, closed_1Z83 and polymer.protein


# ============================================================
# SHOW NATIVE AP5 FROM 1Z83
# ============================================================

select native_AP5, closed_1Z83 and chain A and resn AP5

show sticks, native_AP5

color yellow, native_AP5

set stick_radius, 0.18, native_AP5


# ============================================================
# LOAD ALIGNED 7X7S SELECTED MODELS
# ============================================================

load AK1_7X7S_model20_least_open_ALIGNED_1Z83.pdb, apo_least_open
load AK1_7X7S_model8_ensemble_medoid_ALIGNED_1Z83.pdb, apo_medoid
load AK1_7X7S_model1_most_open_ALIGNED_1Z83.pdb, apo_most_open


hide everything, apo_least_open
hide everything, apo_medoid
hide everything, apo_most_open


show cartoon, apo_least_open
show cartoon, apo_medoid
show cartoon, apo_most_open


# ============================================================
# COLORS
# ============================================================

color palegreen, apo_least_open
color cyan, apo_medoid
color magenta, apo_most_open


# ============================================================
# TRANSPARENCY
# ============================================================

set cartoon_transparency, 0.15, closed_1Z83
set cartoon_transparency, 0.35, apo_least_open
set cartoon_transparency, 0.35, apo_medoid
set cartoon_transparency, 0.35, apo_most_open


# ============================================================
# IMPORTANT OPENING RESIDUES
#
# 1Z83 numbering:
# V47, R132, T145
# ============================================================

select closed_markers, closed_1Z83 and chain A and resi 47+132+145 and name CA

show spheres, closed_markers

color gray30, closed_markers

set sphere_scale, 0.45, closed_markers


# ============================================================
# KEY NUCLEOTIDE-REGION RESIDUES IN CLOSED 1Z83
# ============================================================

select closed_pocket, closed_1Z83 and chain A and resi 15+16+17+18+21+22+38+61+64+93+94+97+132+138+140+149+153

show sticks, closed_pocket

color wheat, closed_pocket

set stick_radius, 0.10, closed_pocket


# ============================================================
# VIEW
# ============================================================

orient closed_1Z83 and polymer.protein

zoom closed_1Z83 and polymer.protein, 7

turn y, 15
turn x, -8


# ============================================================
# RENDER
# ============================================================

ray 2600, 2000

png ..\results\AK1_ensemble_sensitivity\7X7S_selected_models\AK1_1Z83_vs_7X7S_selected_overlay.png, dpi=300

quit