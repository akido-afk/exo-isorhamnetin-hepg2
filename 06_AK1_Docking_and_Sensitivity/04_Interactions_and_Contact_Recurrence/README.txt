04 Interactions and Contact Recurrence
Which AK1 residues contact the representative pose and recur across the dominant pose family?

PAPER PANELS
Figure 8b,d; Figure S14b,c; Figure S15a–c. These panels show the representative interaction map, recurrent polar distances, broad/polar contact recurrence and ISO–AP5 footprint overlap.
Figure 8b, Figure 8d, Figure S14b, Figure S14c, Figure S15a, Figure S15b, Figure S15c

CONTENTS
Polar_Contacts_Recurrence_and_AP5_Overlap [Fig8bd_FigS14bc_FigS15abc].txt — source code and execution instructions
code/ — Python, R, batch and PyMOL scripts
inputs/ — molecular structures and input datasets where applicable
results/ — saved docking results and numerical analyses
figures/ — figure outputs
supporting_files/ — intermediate files needed by the sources
INPUTS_AND_OUTPUTS.txt — script and file correspondence
REPRODUCTION_COMMANDS.txt — commands for running the supplied scripts

METHODS AND PARAMETERS
AK1_extract_ISO_contacts.py uses the repaired receptor and results/AK1_repr9009/AK1_ISO_seed9009_best.sdf. Candidate polar screening uses 3.6 Å; candidate hydrogen-bond distance 3.5 Å is not an angular hydrogen-bond validation. AK1_dominant8_contact_frequency.py uses broad heavy-atom contacts ≤4.0 Å and polar proximity ≤3.6 Å across its explicit eight-seed list. AP5 comparison and recurrent-distance exports have separate scripts.

REPRODUCING THIS CHAPTER
Prepare a working project as described in the docking README.txt. The setup creates the structures/, receptors/, ligand/, docking/ and results/ directories and sets the project and Vina paths. Run each script from structures/.
Use the saved docking outputs for downstream analysis. For new searches, prepare a project with --fresh-docking and run the appropriate docking launcher.
Python scripts: python SCRIPT_NAME.py
R scripts: Rscript SCRIPT_NAME.R
PyMOL scripts: pymol -cq SCRIPT_NAME.pml
Batch scripts: run the chosen .bat from the working project's structures directory.
Each source section in the TXT is a separate script. Use the preparation helper to select the execution scripts and place them in the working project.

FILES
  results/AK1_AP5_ISO_fingerprint/AK1_AP5_vs_ISO_residue_fingerprint.csv
  ../Workbooks/04_AP5_vs_ISO_residue_fingerprint.xlsx
  results/AK1_dominant8_contact_frequency/AK1_dominant8_contact_frequency.csv
  ../Workbooks/04_dominant8_contact_frequency.xlsx
  ../Workbooks/04_dominant8_recurrent_polar_distances.xlsx
  results/AK1_dominant8_contact_frequency/AK1_dominant8_recurrent_polar_distances_long.csv
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed1001_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed1001_best.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed10101_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed10101_best.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed2002_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed2002_best.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed3003_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed3003_best.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed4004_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed4004_best.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed5005_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed5005_best.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed7007_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed7007_best.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed9009_all.sdf
  supporting_files/AK1_dominant8_contact_frequency/poses/AK1_ISO_seed9009_best.sdf
  ../Workbooks/04_ISO_seed9009_contact_analysis.xlsx
  results/AK1_panelB/AK1_ISO_seed9009_polar_contacts.csv
  figures/AK1_panelB/AK1_panelB_closeup_seed9009.png
  figures/AK1_panelB/AK1_panelB_seed9009_CALLOUTS.png
  figures/AK1_panelB/AK1_panelB_seed9009_FINAL.png
  figures/AK1_panelB/AK1_panelB_seed9009_POLAR_CALLOUTS.png
  figures/AK1_panelB/AK1_panelB_seed9009_REFINED_V2.png
  figures/AK1_PanelB_ISO_AP5_interactions.png
  figures/AK1_panelF/AK1_panelF_recurrent_polar_contact_stability_V2.pdf
  figures/AK1_panelF/AK1_panelF_recurrent_polar_contact_stability_V2.png
  figures/AK1_panelG/AK1_panelG_AP5_ISO_polar_fingerprint.pdf
  figures/AK1_panelG/AK1_panelG_AP5_ISO_polar_fingerprint.png
  figures/AK1_panelK_contact_recurrence/AK1_panelK_dominant8_contact_recurrence.pdf
  figures/AK1_panelK_contact_recurrence/AK1_panelK_dominant8_contact_recurrence.png
  results/AK1_panelK_contact_recurrence/AK1_panelK_dominant8_contact_recurrence_data.csv
  ../Workbooks/04_panelK_dominant8_contact_recurrence.xlsx
  figures/AK1_panelL_native_ligand_fingerprint/AK1_panelL_A_polar_contact_fingerprint.pdf
  figures/AK1_panelL_native_ligand_fingerprint/AK1_panelL_A_polar_contact_fingerprint.png
  results/AK1_panelL_native_ligand_fingerprint/AK1_panelL_any_contact_fingerprint.csv
  figures/AK1_panelL_native_ligand_fingerprint/AK1_panelL_B_fingerprint_overlap.pdf
  figures/AK1_panelL_native_ligand_fingerprint/AK1_panelL_B_fingerprint_overlap.png
  results/AK1_panelL_native_ligand_fingerprint/AK1_panelL_fingerprint_overlap_summary.csv
  figures/AK1_panelL_native_ligand_fingerprint/AK1_panelL_native_ligand_fingerprint_comparison.pdf
  figures/AK1_panelL_native_ligand_fingerprint/AK1_panelL_native_ligand_fingerprint_comparison.png
  ../Workbooks/04_panelL_native_ligand_fingerprint.xlsx
  results/AK1_panelL_native_ligand_fingerprint/AK1_panelL_polar_contact_fingerprint.csv
  code/AK1_AP5_ISO_fingerprint_compare.py
  code/AK1_dominant8_contact_frequency.py
  code/AK1_export_recurrent_contact_distances.py
  code/AK1_extract_ISO_contacts.py
  code/AK1_panelB_closeup.pml
  code/AK1_panelB_refined_v2.pml
  code/AK1_panelF_contact_stability.R
  code/AK1_panelF_contact_stability_v2.R
  code/AK1_panelG_AP5_ISO_fingerprint.R
  code/AK1_panelG_contact_recurrence.R
  code/AK1_panelK_contact_recurrence.R
  code/AK1_panelL_native_ligand_fingerprint.R
  code/make_AK1_panelB.pml
  code/make_AK1_panelB_callouts.pml
  code/make_AK1_panelB_callouts_v2.pml

EXECUTION SCRIPTS
The preparation helper automatically selects these files from code/portable/:
AK1_panelF_contact_stability.R — Plot mean polar-contact distances with standard-deviation error bars.

RUNNING WITH THE CONFIGURED ENVIRONMENT
After preparation, work from the new structures/ directory and use:
python ..\run_analysis.py SCRIPT_NAME.py
python ..\run_analysis.py SCRIPT_NAME.R
python ..\run_analysis.py SCRIPT_NAME.pml
Docking launchers additionally require --allow-docking. Use --preview for a rapid PyMOL still-image check; full-resolution settings are the default.

ADDITIONAL EXECUTION SCRIPTS
AK1_panelG_contact_recurrence.R — Plot broad and polar contact recurrence, ordering residues by their shortest polar-contact distance.
