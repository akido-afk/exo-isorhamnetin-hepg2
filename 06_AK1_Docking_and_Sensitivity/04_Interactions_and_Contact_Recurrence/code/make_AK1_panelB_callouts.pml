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
# AP5 REFERENCE
# ============================================================

show sticks, AP5_native

color cyan, AP5_native and elem C
color red, AP5_native and elem O
color blue, AP5_native and elem N

set stick_radius, 0.10, AP5_native
set stick_transparency, 0.60, AP5_native


# ============================================================
# REPRESENTATIVE ISORHAMNETIN — SEED 9009
# ============================================================

show sticks, ISO

color yellow, ISO and elem C
color red, ISO and elem O
color blue, ISO and elem N

set stick_radius, 0.32, ISO

show spheres, ISO
set sphere_scale, 0.17, ISO


# ============================================================
# GENERAL POCKET
# ============================================================

select pocket, byres (AK1 within 4.0 of ISO)

show sticks, pocket

color palegreen, pocket and elem C
color red, pocket and elem O
color blue, pocket and elem N

set stick_radius, 0.11, pocket
set stick_transparency, 0.25, pocket


# ============================================================
# DIRECT-CONTACT RESIDUES
# ============================================================

select direct_contacts, byres (AK1 within 3.6 of ISO)

# make the actual contacting residues slightly stronger
show sticks, direct_contacts
set stick_radius, 0.16, direct_contacts


# ============================================================
# RED CALLOUT MARKERS
# ============================================================

select direct_CA, direct_contacts and name CA

show spheres, direct_CA

color red, direct_CA

set sphere_scale, 0.34, direct_CA
set sphere_transparency, 0.45, direct_CA


# ============================================================
# LABEL ONLY DIRECT-CONTACT RESIDUES
# ============================================================

hide labels

label direct_CA, "%s%s" % (resn,resi)

set label_size, 19
set label_color, red
set label_font_id, 7
set label_outline_color, white


# ============================================================
# ISO LABEL
# ============================================================

pseudoatom ISO_label, ISO

hide everything, ISO_label

label ISO_label, "Isorhamnetin"

set label_size, 22
set label_color, red
set label_font_id, 7
set label_outline_color, white


# ============================================================
# PREDICTED POLAR CONTACTS
# ============================================================

distance ISO_contacts, \
    (ISO and (elem O+N)), \
    (direct_contacts and (elem O+N)), \
    3.6, mode=2

set dash_color, black, ISO_contacts
set dash_width, 2.5, ISO_contacts
set dash_gap, 0.25, ISO_contacts
set dash_length, 0.18, ISO_contacts

hide labels, ISO_contacts


# ============================================================
# BACKGROUND / APPEARANCE
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

viewport 2200,1800


# ============================================================
# OUTPUT — NEW FILE
# ============================================================

ray 3000,2400

png ..\results\AK1_panelB\AK1_panelB_seed9009_CALLOUTS.png, dpi=600