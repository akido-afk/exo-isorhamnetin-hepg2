reinitialize

# ============================================================
# LOAD ORIGINAL 1Z83
# Includes AP5, Zn ions and crystallographic waters
# ============================================================

load AK1_1Z83.pdb, complex

hide everything


# ============================================================
# BASIC SELECTIONS
# ============================================================

select AK1A, complex and chain A and polymer.protein

select AP5A, complex and chain A and resn AP5 and resi 631

select ZN601, complex and chain A and resn ZN and resi 601
select ZN641, complex and chain A and resn ZN and resi 641


# ============================================================
# GENERAL STYLE
# ============================================================

show cartoon, AK1A
color gray80, AK1A
set cartoon_transparency, 0.68, AK1A

show sticks, AP5A
color cyan, AP5A and elem C
color red, AP5A and elem O
color blue, AP5A and elem N
color orange, AP5A and elem P
set stick_radius, 0.18, AP5A

bg_color white
set opaque_background, on
set ray_opaque_background, on
set orthoscopic, on
set antialias, 2
set depth_cue, 0
set ambient, 0.45
set specular, 0.20


# ============================================================
# H1 — ZN601 / AP5 COORDINATION ENVIRONMENT
# ============================================================

show spheres, ZN601
color purple, ZN601
set sphere_scale, 0.52, ZN601

# Key crystallographic waters coordinating ZN601
select ZN601_waters, complex and chain A and resn HOH and resi 845+838+703+704

show spheres, ZN601_waters
color marine, ZN601_waters
set sphere_scale, 0.22, ZN601_waters


# Nearby protein residues for context
select ZN601_nearby_res, AK1A and resi 22+38+93+138

show sticks, ZN601_nearby_res
color palegreen, ZN601_nearby_res and elem C
color red, ZN601_nearby_res and elem O
color blue, ZN601_nearby_res and elem N
set stick_radius, 0.13, ZN601_nearby_res


# Direct Zn601 -> AP5 coordination
distance ZN601_AP5_O2D, ZN601, AP5A and name O2D
distance ZN601_AP5_O1G, ZN601, AP5A and name O1G

set dash_color, black, ZN601_AP5_O2D
set dash_color, black, ZN601_AP5_O1G
set dash_width, 2.5
set dash_gap, 0.25


# Zn601 -> coordinating waters
distance ZN601_W845, ZN601, complex and chain A and resn HOH and resi 845 and name O
distance ZN601_W838, ZN601, complex and chain A and resn HOH and resi 838 and name O
distance ZN601_W703, ZN601, complex and chain A and resn HOH and resi 703 and name O
distance ZN601_W704, ZN601, complex and chain A and resn HOH and resi 704 and name O

set dash_color, gray40, ZN601_W845
set dash_color, gray40, ZN601_W838
set dash_color, gray40, ZN601_W703
set dash_color, gray40, ZN601_W704


# Distance labels
set label_color, black
set label_font_id, 7
set label_size, 16
set label_outline_color, white


# Zn label
label ZN601, "ZN601"

# AP5 label
pseudoatom AP5_label, AP5A
hide everything, AP5_label
label AP5_label, "Native AP5"
set label_color, black, AP5_label
set label_size, 18, AP5_label


# Camera for H1
orient ZN601 or AP5A
center ZN601
zoom (ZN601 or AP5A or ZN601_waters or ZN601_nearby_res), 6

turn y, 8
turn x, -6

viewport 2200,1800

ray 3000,2400

png ..\results\AK1_panelH_Zn\AK1_panelH1_ZN601_AP5_environment.png, dpi=600


# ============================================================
# PREPARE H2
# Hide H1-specific objects
# ============================================================

hide spheres, ZN601
hide spheres, ZN601_waters
hide sticks, ZN601_nearby_res
hide labels, ZN601

disable ZN601_AP5_O2D
disable ZN601_AP5_O1G
disable ZN601_W845
disable ZN601_W838
disable ZN601_W703
disable ZN601_W704
disable AP5_label


# ============================================================
# H2 — ZN641 / PROTEIN COORDINATION ENVIRONMENT
# ============================================================

show spheres, ZN641
color purple, ZN641
set sphere_scale, 0.52, ZN641

select ZN641_res, AK1A and resi 26+36

show sticks, ZN641_res
color palegreen, ZN641_res and elem C
color red, ZN641_res and elem O
color blue, ZN641_res and elem N
set stick_radius, 0.16, ZN641_res


# Waters around ZN641
select ZN641_waters, complex and chain A and resn HOH and resi 840+846

show spheres, ZN641_waters
color marine, ZN641_waters
set sphere_scale, 0.22, ZN641_waters


# Direct coordination
distance ZN641_GLU26, ZN641, AK1A and resi 26 and name OE2
distance ZN641_HIS36, ZN641, AK1A and resi 36 and name ND1

distance ZN641_W840, ZN641, complex and chain A and resn HOH and resi 840 and name O
distance ZN641_W846, ZN641, complex and chain A and resn HOH and resi 846 and name O

set dash_color, black, ZN641_GLU26
set dash_color, black, ZN641_HIS36

set dash_color, gray40, ZN641_W840
set dash_color, gray40, ZN641_W846


label ZN641, "ZN641"
label AK1A and resi 26 and name CA, "GLU26"
label AK1A and resi 36 and name CA, "HIS36"

set label_color, black
set label_size, 18
set label_outline_color, white


# Fade AP5 further because it is contextual here
set stick_transparency, 0.70, AP5A


# Camera for H2
orient ZN641
center ZN641

zoom (ZN641 or ZN641_res or ZN641_waters), 6

turn y, -10
turn x, 5

ray 3000,2400

png ..\results\AK1_panelH_Zn\AK1_panelH2_ZN641_protein_environment.png, dpi=600