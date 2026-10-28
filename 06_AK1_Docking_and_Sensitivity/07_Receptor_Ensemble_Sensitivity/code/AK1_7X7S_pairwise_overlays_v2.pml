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

set cartoon_transparency, 0.10, closed_1Z83


# ============================================================
# AP5 REFERENCE
# ============================================================

select native_AP5, closed_1Z83 and chain A and resn AP5

show sticks, native_AP5

color yellow, native_AP5

set stick_radius, 0.12, native_AP5


# ============================================================
# LOAD 7X7S SELECTED CONFORMERS
# ============================================================

load AK1_7X7S_model20_least_open_ALIGNED_1Z83.pdb, apo_least
load AK1_7X7S_model8_ensemble_medoid_ALIGNED_1Z83.pdb, apo_medoid
load AK1_7X7S_model1_most_open_ALIGNED_1Z83.pdb, apo_most


hide everything, apo_least
hide everything, apo_medoid
hide everything, apo_most


# ============================================================
# COLORS
# ============================================================

color palegreen, apo_least
color cyan, apo_medoid
color magenta, apo_most


set cartoon_transparency, 0.18, apo_least
set cartoon_transparency, 0.18, apo_medoid
set cartoon_transparency, 0.18, apo_most


# ============================================================
# ESTABLISH ONE COMMON CAMERA
# ============================================================

orient closed_1Z83 and polymer.protein

zoom closed_1Z83 and polymer.protein, 9

turn y, 18
turn x, -10


# ============================================================
# PANEL 1
# CLOSED 1Z83 VS LEAST-OPEN 7X7S MODEL 20
# ============================================================

hide cartoon, apo_medoid
hide cartoon, apo_most

show cartoon, apo_least

ray 2400, 1800

png ..\results\AK1_ensemble_sensitivity\7X7S_selected_models\AK1_7X7S_V2_closed_vs_least_open.png, dpi=300


# ============================================================
# PANEL 2
# CLOSED 1Z83 VS ENSEMBLE MEDOID MODEL 8
# ============================================================

hide cartoon, apo_least
hide cartoon, apo_most

show cartoon, apo_medoid

ray 2400, 1800

png ..\results\AK1_ensemble_sensitivity\7X7S_selected_models\AK1_7X7S_V2_closed_vs_medoid.png, dpi=300


# ============================================================
# PANEL 3
# CLOSED 1Z83 VS MOST-OPEN MODEL 1
# ============================================================

hide cartoon, apo_least
hide cartoon, apo_medoid

show cartoon, apo_most

ray 2400, 1800

png ..\results\AK1_ensemble_sensitivity\7X7S_selected_models\AK1_7X7S_V2_closed_vs_most_open.png, dpi=300


quit