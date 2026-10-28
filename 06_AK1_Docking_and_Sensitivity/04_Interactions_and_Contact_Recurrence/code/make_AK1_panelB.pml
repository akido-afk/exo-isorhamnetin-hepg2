reinitialize


# ============================================================
# LOAD STRUCTURES
# ============================================================

load ..\receptors\AK1_chainA_FIXED.pdb, AK1

load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native

load ..\results\AK1_repr9009\AK1_ISO_seed9009_best.sdf, ISO


# ============================================================
# CLEAN START
# ============================================================

remove solvent
hide everything


# ============================================================
# AK1 PROTEIN
# ============================================================

show cartoon, AK1

color gray85, AK1

set cartoon_transparency, 0.65, AK1


# ============================================================
# NATIVE AP5
# Reference ligand only — deliberately subdued
# ============================================================

show sticks, AP5_native

color cyan, AP5_native and elem C
color red, AP5_native and elem O
color blue, AP5_native and elem N

set stick_radius, 0.10, AP5_native
set stick_transparency, 0.55, AP5_native


# ============================================================
# ISORHAMNETIN — REPRESENTATIVE SEED 9009
# Main ligand of interest
# ============================================================

show sticks, ISO

# Make ISO very obvious
color yellow, ISO and elem C
color red, ISO and elem O
color blue, ISO and elem N

set stick_radius, 0.32, ISO

# Small spheres help distinguish ISO from AP5
show spheres, ISO
set sphere_scale, 0.18, ISO


# ============================================================
# AK1 POCKET ONLY
# Important: restrict selection to protein
# ============================================================

select pocket, byres (AK1 within 4.0 of ISO)

show sticks, pocket

color palegreen, pocket and elem C
color red, pocket and elem O
color blue, pocket and elem N

set stick_radius, 0.12, pocket
set stick_transparency, 0.20, pocket


# ============================================================
# POLAR CONTACTS BETWEEN ISO AND AK1
# ============================================================

distance ISO_contacts, \
    (ISO and (elem O+N)), \
    (pocket and (elem O+N)), \
    3.6, mode=2

set dash_color, black, ISO_contacts
set dash_width, 2.5, ISO_contacts
set dash_gap, 0.25, ISO_contacts
set dash_length, 0.18, ISO_contacts

# Do not show PyMOL's automatic distance numbers
hide labels, ISO_contacts


# ============================================================
# RESIDUE LABELS
# Only residues physically close to ISO
# ============================================================

hide labels

select label_residues, \
    (AK1 and name CA and byres (AK1 within 3.8 of ISO))

label label_residues, "%s%s" % (resn,resi)

set label_size, 18
set label_color, black
set label_font_id, 7
set label_outline_color, white


# ============================================================
# WHITE BACKGROUND
# ============================================================

bg_color white

set opaque_background, on
set ray_opaque_background, on

set antialias, 2
set orthoscopic, on
set depth_cue, 0
set specular, 0.20
set ambient, 0.45


# ============================================================
# VIEW
# ============================================================

orient ISO

center ISO

zoom (ISO or pocket), 7

turn y, 12
turn x, -8


# ============================================================
# IMAGE SIZE
# ============================================================

viewport 2200, 1800


# ============================================================
# RENDER
# ============================================================

ray 3000, 2400

png ..\results\AK1_panelB\AK1_panelB_seed9009_FINAL.png, dpi=600