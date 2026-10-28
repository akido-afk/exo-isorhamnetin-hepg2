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
set cartoon_transparency, 0.70, AK1


# ============================================================
# AP5 — FAINT REFERENCE
# ============================================================

show sticks, AP5_native

color cyan, AP5_native and elem C
color red, AP5_native and elem O
color blue, AP5_native and elem N

set stick_radius, 0.08, AP5_native
set stick_transparency, 0.65, AP5_native


# ============================================================
# REPRESENTATIVE ISORHAMNETIN — SEED 9009
# ============================================================

show sticks, ISO

color yellow, ISO and elem C
color red, ISO and elem O
color blue, ISO and elem N

set stick_radius, 0.30, ISO

show spheres, ISO
set sphere_scale, 0.15, ISO


# ============================================================
# GENERAL AK1 POCKET
# ============================================================

select pocket, byres (AK1 within 4.0 of ISO)

show sticks, pocket

color palegreen, pocket and elem C
color red, pocket and elem O
color blue, pocket and elem N

set stick_radius, 0.10, pocket
set stick_transparency, 0.35, pocket


# ============================================================
# CANDIDATE POLAR-CONTACT RESIDUES
#
# Much more selective than "anything within 3.6 A".
# Protein O/N atoms must be near ISO O/N atoms.
# ============================================================

select polar_res, byres ( \
    (AK1 and (elem O+N)) within 3.6 of \
    (ISO and (elem O+N)) \
)

show sticks, polar_res

# Make contacting residues stronger than surrounding pocket
color salmon, polar_res and elem C
color red, polar_res and elem O
color blue, polar_res and elem N

set stick_radius, 0.18, polar_res


# ============================================================
# SMALL RED CALLOUT MARKERS
# ============================================================

select polar_CA, polar_res and name CA

show spheres, polar_CA
color red, polar_CA

# deliberately small
set sphere_scale, 0.16, polar_CA
set sphere_transparency, 0.25, polar_CA


# ============================================================
# LABEL ONLY POLAR-CONTACT CANDIDATES
# ============================================================

hide labels

label polar_CA, "%s%s" % (resn,resi)

set label_size, 20
set label_color, black
set label_font_id, 7
set label_outline_color, white


# ============================================================
# PREDICTED POLAR CONTACTS
# ============================================================

distance ISO_polar_contacts, \
    (ISO and (elem O+N)), \
    (polar_res and (elem O+N)), \
    3.6, mode=2

set dash_color, black, ISO_polar_contacts
set dash_width, 2.5, ISO_polar_contacts
set dash_gap, 0.25, ISO_polar_contacts
set dash_length, 0.18, ISO_polar_contacts

# keep diagram clean: no automatic distance-number labels
hide labels, ISO_polar_contacts


# ============================================================
# BACKGROUND / LIGHTING
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

zoom (ISO or pocket), 8

turn y, 12
turn x, -8


# ============================================================
# IMAGE SIZE
# ============================================================

viewport 2200,1800


# ============================================================
# NEW OUTPUT — DOES NOT OVERWRITE OLD FIGURES
# ============================================================

ray 3000,2400

png ..\results\AK1_panelB\AK1_panelB_seed9009_POLAR_CALLOUTS.png, dpi=600