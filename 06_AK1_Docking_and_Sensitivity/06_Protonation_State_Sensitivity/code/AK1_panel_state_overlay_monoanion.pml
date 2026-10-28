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
# We take MODEL 1 = best pose
# Neutral alternative-family representative: seed 8008
# Monoanion representative: seed 10101
# =========================
load ..\docking\AK1_ISO_seed8008.pdbqt, neutral_alt_all
load ..\docking\AK1_ISO_7Ominus_seed10101.pdbqt, monoanion_all

create neutral_alt, neutral_alt_all, 1, 1
create monoanion, monoanion_all, 1, 1

delete neutral_alt_all
delete monoanion_all

# =========================
# SHOW LIGANDS
# =========================
show sticks, neutral_alt
show sticks, monoanion

color cyan, neutral_alt
color hotpink, monoanion

set stick_radius, 0.18, neutral_alt
set stick_radius, 0.18, monoanion

# =========================
# DEFINE LOCAL POCKET
# =========================
select lig_pair, neutral_alt or monoanion
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
png ..\results\AK1_ligand_state_sensitivity\pymol_panels\AK1_panel_state_overlay_monoanion_vs_neutralAlt.png, dpi=300