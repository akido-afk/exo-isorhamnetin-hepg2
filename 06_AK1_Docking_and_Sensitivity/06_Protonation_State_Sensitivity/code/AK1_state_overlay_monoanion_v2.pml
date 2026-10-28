reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on
set depth_cue, 0
set specular, 0.15
set stick_quality, 16
set sphere_quality, 3
set cartoon_fancy_helices, 1
set cartoon_sampling, 14

# =========================
# LOAD RECEPTOR
# =========================
load AK1_1Z83.pdb, prot

hide everything
show cartoon, prot and chain A and polymer.protein
color gray85, prot
set cartoon_transparency, 0.45, prot

# =========================
# LOAD LIGANDS
# Neutral alternative representative = seed 6006
# Monoanion representative = seed 10101
# =========================
load ..\docking\AK1_ISO_seed6006.pdbqt, neutral_alt_all
create neutral_alt, neutral_alt_all, 1, 1
delete neutral_alt_all

load ..\docking\AK1_ISO_7Ominus_seed10101.pdbqt, mono_all
create mono_7O, mono_all, 1, 1
delete mono_all

# =========================
# DEFINE POCKET RESIDUES
# =========================
select pocket, prot and chain A and resi 21+61+97+132+138+149
show sticks, pocket
color palegreen, pocket
set stick_radius, 0.11, pocket

# =========================
# STYLE LIGANDS
# =========================
show sticks, neutral_alt
show sticks, mono_7O

color cyan, neutral_alt
color magenta, mono_7O

set stick_radius, 0.22, neutral_alt
set stick_radius, 0.16, mono_7O

# make the overlap clearer
set stick_transparency, 0.10, neutral_alt
set stick_transparency, 0.00, mono_7O

# =========================
# CLEAN UP
# =========================
hide everything, hydrogens
hide labels
remove solvent

# =========================
# VIEW
# Use the same pocket-centered orientation in both panels
# =========================
orient pocket
zoom pocket or neutral_alt or mono_7O, 11
turn y, 20
turn x, -12

# =========================
# RENDER
# =========================
ray 2400, 1800
png ..\results\AK1_ligand_state_sensitivity\figures_v2\AK1_state_overlay_monoanion_v2.png, dpi=300
quit