reinitialize

# =========================================
# Load protein
# =========================================
load ..\receptors\AK1_chainA_FIXED.pdb, AK1
remove solvent
hide everything

# =========================================
# Load representative ISO pose (seed 9009)
# =========================================
load ..\results\AK1_repr9009\AK1_ISO_seed9009_best.sdf, ISO_9009

# =========================================
# Optional: load native AP5 as faint reference
# =========================================
load ..\ligand\AK1_native_AP5_chainA.pdb, AP5_native

# =========================================
# Define pocket around ISO
# =========================================
select pocket, byres (AK1 within 4.0 of ISO_9009)
select pocket_bb, pocket and backbone
select pocket_sc, pocket and not backbone

# =========================================
# Display settings
# =========================================
hide everything, AK1

show cartoon, AK1
color gray85, AK1
set cartoon_transparency, 0.65, AK1

show sticks, pocket_sc
color palegreen, pocket_sc
set stick_radius, 0.16, pocket_sc

show sticks, ISO_9009
color yellow, ISO_9009 and elem C
color red, ISO_9009 and elem O
color blue, ISO_9009 and elem N
set stick_radius, 0.28, ISO_9009

# AP5 shown lightly for reference only
show lines, AP5_native
color cyan, AP5_native
set line_width, 2.0, AP5_native
set transparency, 0.55, AP5_native

# =========================================
# Residue labels
# =========================================
label pocket and name CA, "%s%s" % (resn, resi)
set label_size, 18
set label_color, black
set label_font_id, 7
set label_outline_color, white

# =========================================
# Polar contacts / H-bond-like distances
# =========================================
# between ISO hetero atoms and nearby protein hetero atoms
distance iso_contacts, (ISO_9009 and (elem O+N)), (pocket_sc and (elem O+N)), 3.6, mode=2
set dash_color, black, iso_contacts
set dash_width, 2.0, iso_contacts
set dash_gap, 0.25, iso_contacts
set dash_length, 0.18, iso_contacts

hide labels, iso_contacts

# =========================================
# Nice appearance
# =========================================
bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on
set depth_cue, 0
set specular, 0.25
set ambient, 0.45

# =========================================
# Camera / focus
# =========================================
orient pocket or ISO_9009
zoom pocket or ISO_9009, 6

# small manual rotation for prettier view
turn x, -10
turn y, 20
turn z, -8

viewport 2200,1800

# =========================================
# Render output
# =========================================
ray 2600,2000
png ..\results\AK1_panelB\AK1_panelB_closeup_seed9009.png, dpi=600