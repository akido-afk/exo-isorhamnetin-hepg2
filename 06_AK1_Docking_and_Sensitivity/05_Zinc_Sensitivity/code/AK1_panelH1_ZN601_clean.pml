reinitialize

# ============================================================
# LOAD ORIGINAL 1Z83 — SAME FILE THAT WORKED BEFORE
# ============================================================

load AK1_1Z83.pdb, complex

hide everything

# ============================================================
# SELECTIONS
# ============================================================

select AK1A, complex and chain A and polymer.protein
select AP5A, complex and chain A and resn AP5 and resi 631
select ZN601, complex and chain A and resn ZN and resi 601

# ============================================================
# PROTEIN
# ============================================================

show cartoon, AK1A
color gray80, AK1A
set cartoon_transparency, 0.68, AK1A

# ============================================================
# AP5
# ============================================================

show sticks, AP5A

color cyan, AP5A and elem C
color red, AP5A and elem O
color blue, AP5A and elem N
color orange, AP5A and elem P

set stick_radius, 0.18, AP5A

# ============================================================
# ZN601
# ============================================================

show spheres, ZN601
color purple, ZN601
set sphere_scale, 0.52, ZN601

label ZN601, "ZN601"

# ============================================================
# FOUR COORDINATING WATERS
# Show them, but DO NOT draw water-distance labels
# ============================================================

select ZN601_waters, complex and chain A and resn HOH and resi 845+838+703+704

show spheres, ZN601_waters
color marine, ZN601_waters
set sphere_scale, 0.22, ZN601_waters

# ============================================================
# NEARBY PROTEIN RESIDUES
# ============================================================

select ZN601_nearby_res, AK1A and resi 22+38+93+138

show sticks, ZN601_nearby_res

color palegreen, ZN601_nearby_res and elem C
color red, ZN601_nearby_res and elem O
color blue, ZN601_nearby_res and elem N

set stick_radius, 0.13, ZN601_nearby_res

# ============================================================
# ONLY THE TWO IMPORTANT AP5–ZN601 DISTANCES
# ============================================================

distance ZN601_AP5_O2D, ZN601, AP5A and name O2D
distance ZN601_AP5_O1G, ZN601, AP5A and name O1G

set dash_color, black, ZN601_AP5_O2D
set dash_color, black, ZN601_AP5_O1G

set dash_width, 2.5
set dash_gap, 0.25

# ============================================================
# LABEL STYLE
# ============================================================

set label_color, black
set label_font_id, 7
set label_size, 18
set label_outline_color, white

# ============================================================
# APPEARANCE
# ============================================================

bg_color white

set opaque_background, on
set ray_opaque_background, on

set antialias, 2
set orthoscopic, on
set depth_cue, 0
set ambient, 0.45
set specular, 0.20

# ============================================================
# CAMERA
# ============================================================

orient ZN601 or AP5A

center ZN601

zoom (ZN601 or AP5A or ZN601_waters or ZN601_nearby_res), 6

turn y, 8
turn x, -6

viewport 2200,1800

# ============================================================
# OUTPUT — CORRECT EXISTING FOLDER
# ============================================================

ray 3000,2400

png ..\results\AK1_panelH_Zn\AK1_panelH1_ZN601_clean.png, dpi=600