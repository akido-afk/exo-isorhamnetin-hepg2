reinitialize

# =========================
# Load protein
# =========================
load ..\receptors\AK1_chainA_FIXED.pdb, AK1
remove solvent
hide everything

show cartoon, AK1
color gray70, AK1
set cartoon_transparency, 0.35, AK1

# =========================
# Load native AP5
# =========================
load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native
show sticks, AP5_native
color cyan, AP5_native
set stick_radius, 0.16, AP5_native

# =========================
# Load representative ISO pose
# seed 9009
# =========================
load ..\results\AK1_repr9009\AK1_ISO_seed9009_best.sdf, ISO_9009
show sticks, ISO_9009

color yellow, ISO_9009 and elem C
color red, ISO_9009 and elem O
color blue, ISO_9009 and elem N

set stick_radius, 0.24, ISO_9009

# =========================
# Binding pocket
# =========================
select pocket, byres (AK1 within 4.0 of ISO_9009)
show sticks, pocket
color palegreen, pocket
set stick_radius, 0.12, pocket
set stick_transparency, 0.20, pocket

# =========================
# General appearance
# =========================
bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on

viewport 1800,1800

# =========================
# Camera / framing
# =========================
orient AK1
center AK1
zoom AK1, 18

# =========================
# Render final PNG
# =========================
ray 2800,2800
png ..\results\AK1_repr9009\AK1_ISO_representative_seed9009_static.png, dpi=600