reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 2
set orthoscopic, on

# -----------------------------
# Load structures
# -----------------------------
load ..\receptors\AK1_chainA_FIXED.pdb, AK1
load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native

# IMPORTANT:
# Replace this filename if your extracted representative pose has another name.
# I assume you extracted seed 9009 best pose as an SDF.
load ..\docking\AK1_ISO_seed9009_bestpose.sdf, ISO

remove solvent
hide everything

# -----------------------------
# Protein style
# -----------------------------
show cartoon, AK1
color gray70, AK1
set cartoon_transparency, 0.35, AK1

# -----------------------------
# Native ligand (reference) = cyan
# -----------------------------
show sticks, AP5_native
color cyan, AP5_native
set stick_radius, 0.16, AP5_native

# -----------------------------
# Docked ISO = bright orange/yellow and thicker
# -----------------------------
show sticks, ISO
color tv_orange, ISO and elem C
color red, ISO and elem O
color blue, ISO and elem N
set stick_radius, 0.24, ISO

# -----------------------------
# Pocket residues
# -----------------------------
select pocket, byres ((AP5_native or ISO) around 4.0)
show sticks, pocket
color palegreen, pocket
set stick_radius, 0.14, pocket

# Optional: make ISO slightly "pop" more
set stick_transparency, 0.00, ISO
set stick_transparency, 0.15, AP5_native
set stick_transparency, 0.10, pocket

# -----------------------------
# View/orientation
# -----------------------------
orient AP5_native or ISO
zoom AP5_native or ISO or pocket, 11

# -----------------------------
# Save one static PNG first
# -----------------------------
ray 2200,2200
png ..\results\AK1_rotate_ISO_orange\AK1_ISO_orange_static.png, dpi=300

# -----------------------------
# Movie frames
# -----------------------------
mset 1 x120
python
for i in range(1, 121):
    cmd.mdo(i, "turn y, 3")
python end

set ray_trace_frames, 1
mpng ..\results\AK1_rotate_ISO_orange\frames\frame_