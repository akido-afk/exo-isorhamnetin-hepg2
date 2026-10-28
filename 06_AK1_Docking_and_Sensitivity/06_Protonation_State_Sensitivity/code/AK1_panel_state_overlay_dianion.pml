reinitialize

bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on
set depth_cue, off
set specular, 0.2
set stick_radius, 0.16

# =========================
# LOAD RECEPTOR
# =========================
load ..\structures\AK1_1Z83.pdb, AK1
remove not chain A
remove solvent
hide everything

show cartoon, polymer.protein
color gray85, polymer.protein
set cartoon_transparency, 0.55

# =========================
# LOAD DOCKED LIGANDS
# Neutral dominant representative: seed 9009
# Dianion representative: seed 9009
# =========================
load ..\docking\AK1_ISO_seed9009.pdbqt, neutral_dom_all
load ..\docking\AK1_ISO_dianion_seed9009.pdbqt, dianion_all

create neutral_dom, neutral_dom_all, 1, 1
create dianion, dianion_all, 1, 1

delete neutral_dom_all
delete dianion_all

# =========================
# SHOW LIGANDS
# =========================
show sticks, neutral_dom
show sticks, dianion

color cyan, neutral_dom
color orange, dianion

set stick_radius, 0.18, neutral_dom
set stick_radius, 0.18, dianion

# =========================
# DEFINE LOCAL POCKET
# =========================
select lig_pair, neutral_dom or dianion
select pocket, byres (polymer.protein within 4.5 of lig_pair)

show sticks, pocket
color palegreen, pocket
set stick_radius, 0.15, pocket

# Optional: label key recurrent residues
label (pocket and resi 21 and name CA), "LYS21"
label (pocket and resi 61 and name CA), "MET61"
label (pocket and resi 97 and name CA), "ARG97"
label (pocket and resi 132 and name CA), "ARG132"
label (pocket and resi 138 and name CA), "ARG138"
label (pocket and resi 149 and name CA), "ARG149"

set label_size, 18
set label_color, black

# =========================
# VIEW
# =========================
orient lig_pair
zoom lig_pair, 10
turn y, 18
turn x, -12

# Clean up
hide lines
hide nonbonded
set cartoon_fancy_helices, 1

# =========================
# OUTPUT
# =========================
ray 2400, 1800
png ..\results\AK1_ligand_state_sensitivity\pymol_panels\AK1_panel_state_overlay_dianion_vs_neutralDom.png, dpi=300