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
# Neutral dominant representative = seed 9009
# Dianion representative = seed 9009
# =========================
load ..\docking\AK1_ISO_seed9009.pdbqt, neutral_dom_all
create neutral_dom, neutral_dom_all, 1, 1
delete neutral_dom_all

load ..\docking\AK1_ISO_dianion_seed9009.pdbqt, dianion_all
create dianion_7O_4O, dianion_all, 1, 1
delete dianion_all

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
show sticks, neutral_dom
show sticks, dianion_7O_4O

color cyan, neutral_dom
color orange, dianion_7O_4O

set stick_radius, 0.22, neutral_dom
set stick_radius, 0.16, dianion_7O_4O

# make both visible in overlap
set stick_transparency, 0.10, neutral_dom
set stick_transparency, 0.00, dianion_7O_4O

# =========================
# CLEAN UP
# =========================
hide everything, hydrogens
hide labels
remove solvent

# =========================
# VIEW
# same orientation as monoanion panel
# =========================
orient pocket
zoom pocket or neutral_dom or dianion_7O_4O, 11
turn y, 20
turn x, -12

# =========================
# RENDER
# =========================
ray 2400, 1800
png ..\results\AK1_ligand_state_sensitivity\figures_v2\AK1_state_overlay_dianion_v2.png, dpi=300
quit