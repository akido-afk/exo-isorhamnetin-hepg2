# ============================================================
# AK1 + representative ISO pose
# Representative dominant pose = seed 9009
# ============================================================

reinitialize

# -----------------------------
# AK1 protein
# -----------------------------

load ..\receptors\AK1_chainA_FIXED.pdb, AK1

remove solvent
hide everything

show cartoon, AK1
color gray70, AK1

set cartoon_transparency, 0.30, AK1


# -----------------------------
# Native AP5
# -----------------------------

load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native

show sticks, AP5_native

color cyan, AP5_native

set stick_radius, 0.16, AP5_native


# -----------------------------
# Representative ISO
# seed 9009
# -----------------------------

load ..\results\AK1_final_visuals\poses\AK1_ISO_seed9009_best.sdf, ISO_rep

show sticks, ISO_rep

# ligand carbons
color yellow, ISO_rep and elem C

# oxygens
color red, ISO_rep and elem O

# nitrogen if present
color blue, ISO_rep and elem N

set stick_radius, 0.28, ISO_rep


# -----------------------------
# Binding pocket
# residues within 4 A of ISO
# -----------------------------

select ISO_pocket, byres (AK1 within 4.0 of ISO_rep)

show sticks, ISO_pocket

color palegreen, ISO_pocket

set stick_radius, 0.13, ISO_pocket

set stick_transparency, 0.20, ISO_pocket


# -----------------------------
# General appearance
# -----------------------------

bg_color white

set antialias, 2
set ray_opaque_background, on

set orthoscopic, on

viewport 1600, 1600


# -----------------------------
# Frame WHOLE AK1
# -----------------------------

orient AK1
center AK1

# Extra margin prevents clipping during rotation
zoom AK1, 18


# ============================================================
# STATIC HIGH-RESOLUTION PNG
# ============================================================

ray 2800, 2800

png ..\results\AK1_final_visuals\AK1_ISO_representative_static.png, dpi=600


# ============================================================
# ROTATION
# ============================================================

mset 1 x120

python
for i in range(1, 121):
    cmd.mdo(i, "turn y, 3")
python end

set ray_trace_frames, 1

mpng ..\results\AK1_final_visuals\rotation_frames\AK1_frame_