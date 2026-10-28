reinitialize

# ============================================================
# LOAD STRUCTURES
# ============================================================

load ..\receptors\AK1_chainA_FIXED.pdb, AK1

load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native

load ..\docking\AK1_AP5_redock_bestpose.sdf, AP5_redocked


# ============================================================
# CLEAN START
# ============================================================

remove solvent
hide everything


# ============================================================
# AK1 PROTEIN
# ============================================================

show cartoon, AK1

color gray80, AK1

set cartoon_transparency, 0.68, AK1


# ============================================================
# BINDING POCKET
# ============================================================

select pocket, byres (AK1 within 4.0 of AP5_native)

show sticks, pocket

color palegreen, pocket and elem C
color red, pocket and elem O
color blue, pocket and elem N

set stick_radius, 0.10, pocket
set stick_transparency, 0.35, pocket


# ============================================================
# CRYSTALLOGRAPHIC AP5
# CYAN
# ============================================================

show sticks, AP5_native

color cyan, AP5_native and elem C
color red, AP5_native and elem O
color blue, AP5_native and elem N

set stick_radius, 0.23, AP5_native


# ============================================================
# REDOCKED AP5
# MAGENTA
# ============================================================

show sticks, AP5_redocked

color magenta, AP5_redocked and elem C
color red, AP5_redocked and elem O
color blue, AP5_redocked and elem N

set stick_radius, 0.18, AP5_redocked


# ============================================================
# APPEARANCE
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
# CLOSE-UP VIEW
# ============================================================

orient AP5_native

center AP5_native

zoom (AP5_native or AP5_redocked or pocket), 7

turn y, 10
turn x, -5


# ============================================================
# IMAGE SIZE
# ============================================================

viewport 2200,1800


# ============================================================
# OUTPUT
# ============================================================

ray 3000,2400

png ..\results\AK1_panelC\AK1_panelC_AP5_redocking_validation.png, dpi=600