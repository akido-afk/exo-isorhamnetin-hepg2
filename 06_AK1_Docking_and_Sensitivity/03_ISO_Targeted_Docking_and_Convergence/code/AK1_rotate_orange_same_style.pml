load ..\receptors\AK1_chainA_FIXED.pdb, AK1
remove solvent
hide everything

# -----------------------------
# AK1 protein
# -----------------------------
show cartoon, AK1
color gray70, AK1
set cartoon_transparency, 0.35, AK1

# -----------------------------
# Native AP5
# -----------------------------
load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native
show sticks, AP5_native
color cyan, AP5_native

# -----------------------------
# Redocked AP5
# -----------------------------
load ..\docking\AK1_AP5_redock_bestpose.sdf, AP5_redocked
show sticks, AP5_redocked
color magenta, AP5_redocked

# -----------------------------
# Docked isorhamnetin
# -----------------------------
load ..\docking\AK1_ISO_bestpose.sdf, ISO
show sticks, ISO
color tv_orange, ISO

# -----------------------------
# Binding-pocket residues
# -----------------------------
select pocket, byres (AP5_native around 4.0)
show sticks, pocket
color palegreen, pocket

# recolor ISO after pocket coloring so it remains visible
color tv_orange, ISO

set stick_radius, 0.18
set stick_radius, 0.24, ISO

# -----------------------------
# Appearance
# -----------------------------
bg_color white
set ray_opaque_background, on
set antialias, 2

# -----------------------------
# IMPORTANT: wider framing
# -----------------------------
viewport 1600, 1600

orient AK1
center AK1

# 15 A buffer around whole AK1
zoom AK1, 15

# -----------------------------
# Static figure
# -----------------------------
ray 2400, 2400
png ..\results\AK1_rotate_ISO_orange\AK1_ISO_AP5_overlay_orange.png, dpi=300

# -----------------------------
# Rotating movie
# -----------------------------
mset 1 x120

python
for i in range(1, 121):
    cmd.mdo(i, "turn y, 3")
python end

set ray_trace_frames, 1

mpng ..\results\AK1_rotate_ISO_orange\frames\frame_