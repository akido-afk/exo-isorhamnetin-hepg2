reinitialize

# ============================================================
# BACKGROUND / RENDER SETTINGS
# ============================================================

bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on
set depth_cue, 0

# ============================================================
# LOAD AK1 + ZN601 RECEPTOR
# ============================================================

load ..\receptors\AK1_chainA_FIXED_ZN601.pdb, AK1_ZN

hide everything

select AK1, AK1_ZN and polymer.protein
select ZN601, AK1_ZN and resn ZN and resi 601

show cartoon, AK1
color gray80, AK1
set cartoon_transparency, 0.72, AK1

show spheres, ZN601
color purple, ZN601
set sphere_scale, 0.52, ZN601

# ============================================================
# LOAD ORIGINAL Zn-FREE DOMINANT REPRESENTATIVE
# seed 9009
# ============================================================

load ..\results\AK1_repr9009\AK1_ISO_seed9009_best.sdf, ISO_Znfree

show sticks, ISO_Znfree

color tv_orange, ISO_Znfree and elem C
color red, ISO_Znfree and elem O
color blue, ISO_Znfree and elem N

set stick_radius, 0.22, ISO_Znfree

# ============================================================
# LOAD ZN601-RETAINED REPRESENTATIVE
# seed 5005
# ============================================================

load ..\results\AK1_ZN_sensitivity\ZN601_10seed\poses\AK1_ZN601_ISO_seed5005_best.sdf, ISO_ZN601

show sticks, ISO_ZN601

color marine, ISO_ZN601 and elem C
color red, ISO_ZN601 and elem O
color blue, ISO_ZN601 and elem N

set stick_radius, 0.22, ISO_ZN601

# ============================================================
# NATIVE AP5 AS FAINT REFERENCE
# ============================================================

load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native

show sticks, AP5_native

color cyan, AP5_native and elem C
color red, AP5_native and elem O
color blue, AP5_native and elem N
color orange, AP5_native and elem P

set stick_radius, 0.12, AP5_native
set stick_transparency, 0.72, AP5_native

# ============================================================
# LOCAL BINDING-SITE RESIDUES
# ============================================================

select pocket, byres ((ISO_Znfree around 4.0) or (ISO_ZN601 around 4.0))

select pocket_protein, pocket and polymer.protein

show sticks, pocket_protein

color palegreen, pocket_protein and elem C
color red, pocket_protein and elem O
color blue, pocket_protein and elem N

set stick_radius, 0.11, pocket_protein
set stick_transparency, 0.25, pocket_protein

# ============================================================
# ZN601 -> ISO CONTACT
# Find oxygen atoms of Zn-retained ISO within 2.8 A
# ============================================================

select ISO_ZN_O, ISO_ZN601 and elem O within 2.8 of ZN601

distance ZN601_ISO_contact, ZN601, ISO_ZN_O

set dash_color, black, ZN601_ISO_contact
set dash_width, 3
set dash_gap, 0.25

set label_color, black
set label_size, 18
set label_outline_color, white
set label_font_id, 7

# ============================================================
# LABEL ZINC ONLY
# ============================================================

label ZN601, "ZN601"

# ============================================================
# CAMERA
# ============================================================

orient ISO_Znfree or ISO_ZN601 or ZN601
center ISO_Znfree or ISO_ZN601 or ZN601
zoom ISO_Znfree or ISO_ZN601 or ZN601 or pocket_protein, 7

turn y, 8
turn x, -5

viewport 2200,1800

# ============================================================
# FINAL RENDER
# ============================================================

ray 3000,2400

png ..\results\AK1_panelJ\AK1_panelJ_Znfree_vs_ZN601_overlay.png, dpi=600