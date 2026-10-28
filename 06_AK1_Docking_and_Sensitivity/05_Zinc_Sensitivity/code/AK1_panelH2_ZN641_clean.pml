reinitialize

# ============================================================
# LOAD ORIGINAL 1Z83
# ============================================================

load AK1_1Z83.pdb, complex

hide everything

# ============================================================
# SELECTIONS
# ============================================================

select AK1A, complex and chain A and polymer.protein
select AP5A, complex and chain A and resn AP5 and resi 631
select ZN641, complex and chain A and resn ZN and resi 641

# ============================================================
# PROTEIN
# ============================================================

show cartoon, AK1A
color gray80, AK1A
set cartoon_transparency, 0.68, AK1A

# ============================================================
# ZN641
# ============================================================

show spheres, ZN641
color purple, ZN641
set sphere_scale, 0.52, ZN641

label ZN641, "ZN641"

# ============================================================
# GLU26 / HIS36
# ============================================================

select ZN641_res, AK1A and resi 26+36

show sticks, ZN641_res

color palegreen, ZN641_res and elem C
color red, ZN641_res and elem O
color blue, ZN641_res and elem N

set stick_radius, 0.16, ZN641_res

label AK1A and resi 26 and name CA, "GLU26"
label AK1A and resi 36 and name CA, "HIS36"

# ============================================================
# COORDINATING WATERS
# Show, but no water-distance labels
# ============================================================

select ZN641_waters, complex and chain A and resn HOH and resi 840+846

show spheres, ZN641_waters
color marine, ZN641_waters
set sphere_scale, 0.22, ZN641_waters

# ============================================================
# ONLY TWO IMPORTANT PROTEIN–ZN641 DISTANCES
# ============================================================

distance ZN641_GLU26, ZN641, AK1A and resi 26 and name OE2
distance ZN641_HIS36, ZN641, AK1A and resi 36 and name ND1

set dash_color, black, ZN641_GLU26
set dash_color, black, ZN641_HIS36

set dash_width, 2.5
set dash_gap, 0.25

# ============================================================
# AP5 ONLY AS FAINT CONTEXT
# ============================================================

show sticks, AP5A

color cyan, AP5A and elem C
color red, AP5A and elem O
color blue, AP5A and elem N
color orange, AP5A and elem P

set stick_radius, 0.13, AP5A
set stick_transparency, 0.72, AP5A

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

orient ZN641

center ZN641

zoom (ZN641 or ZN641_res or ZN641_waters), 6

turn y, -10
turn x, 5

viewport 2200,1800

# ============================================================
# OUTPUT
# ============================================================

ray 3000,2400

png ..\results\AK1_panelH_Zn\AK1_panelH2_ZN641_clean.png, dpi=600