05 Zinc Sensitivity
How does retaining crystallographic ZN601 change docking geometry and scores?

PAPER PANELS
Figure S15d–g: zinc environments, scores, retained-zinc pose RMSD and comparison with Zn-free families. Figure S16: representative Zn-free versus ZN601-retained pose overlay.
Figure S15d, Figure S15e, Figure S15f, Figure S15g, Figure S16

CONTENTS
Zinc_Environment_and_ZN601_Docking_Sensitivity [FigS15defg_FigS16].txt — source code and execution instructions
code/ — Python, R, batch and PyMOL scripts
inputs/ — molecular structures and input datasets where applicable
results/ — saved docking results and numerical analyses
figures/ — figure outputs
supporting_files/ — intermediate files needed by the sources
INPUTS_AND_OUTPUTS.txt — script and file correspondence
REPRODUCTION_COMMANDS.txt — commands for running the supplied scripts

METHODS AND PARAMETERS
build_AK1_FIXED_ZN601.py and AK1_inspect_Zn_environment.py establish the retained-ion model. run_AK1_ZN601_10seeds.bat performs the sensitivity search; AK1_ZN601_10seed_comprehensive.py calculates multi-seed pose and score comparisons. ZN601 and ZN641 environment illustrations are preserved separately. These geometries are Vina sensitivity outputs, not a validated zinc coordination force field.

REPRODUCING THIS CHAPTER
Prepare a working project as described in the docking README.txt. The setup creates the structures/, receptors/, ligand/, docking/ and results/ directories and sets the project and Vina paths. Run each script from structures/.
Use the saved docking outputs for downstream analysis. For new searches, prepare a project with --fresh-docking and run the appropriate docking launcher.
Python scripts: python SCRIPT_NAME.py
R scripts: Rscript SCRIPT_NAME.R
PyMOL scripts: pymol -cq SCRIPT_NAME.pml
Batch scripts: run the chosen .bat from the working project's structures directory.
Each source section in the TXT is a separate script. Use the preparation helper to select the execution scripts and place them in the working project.

FILES
  results/docking/AK1_ZN601_ISO_seed1001.pdbqt
  results/docking/AK1_ZN601_ISO_seed10101.pdbqt
  results/docking/AK1_ZN601_ISO_seed2002.pdbqt
  results/docking/AK1_ZN601_ISO_seed3003.pdbqt
  results/docking/AK1_ZN601_ISO_seed4004.pdbqt
  results/docking/AK1_ZN601_ISO_seed5005.pdbqt
  results/docking/AK1_ZN601_ISO_seed6006.pdbqt
  results/docking/AK1_ZN601_ISO_seed7007.pdbqt
  results/docking/AK1_ZN601_ISO_seed8008.pdbqt
  results/docking/AK1_ZN601_ISO_seed9009.pdbqt
  inputs/receptors/AK1_chainA_FIXED_ZN601.pdb
  inputs/receptors/AK1_ZN601.pdbqt
  figures/AK1_panelH_Zn/AK1_panelH1_ZN601_AP5_environment.png
  figures/AK1_panelH_Zn/AK1_panelH1_ZN601_clean.png
  figures/AK1_panelH_Zn/AK1_panelH2_ZN641_clean.png
  figures/AK1_panelH_Zn/AK1_panelH2_ZN641_protein_environment.png
  figures/AK1_panelH_Zn/AK1_panelH_Zn_environment_FINAL.png
  figures/AK1_panelI/AK1_panelI_ZN601_receptor_state_sensitivity.pdf
  figures/AK1_panelI/AK1_panelI_ZN601_receptor_state_sensitivity.png
  figures/AK1_panelJ/AK1_panelJ_Znfree_vs_ZN601_overlay.png
  ../Workbooks/05_1Z83_Zn_environment_analysis.xlsx
  results/AK1_1Z83_Zn_neighbors_within5A.csv
  supporting_files/AK1_ISO_seed9009_ZN601_allposes.sdf
  supporting_files/AK1_ISO_seed9009_ZN601_best.sdf
  supporting_files/AK1_ISO_seed9009_ZNfree_allposes.sdf
  supporting_files/AK1_ISO_seed9009_ZNfree_best.sdf
  ../Workbooks/05_ZN601_pilot_pose_classification.xlsx
  ../Workbooks/05_ZN601_seed9009_pilot_comparison.xlsx
  supporting_files/original_10_reference_poses/AK1_ISO_seed1001_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed1001_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed10101_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed10101_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed2002_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed2002_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed3003_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed3003_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed4004_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed4004_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed5005_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed5005_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed6006_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed6006_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed7007_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed7007_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed8008_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed8008_best.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed9009_all.sdf
  supporting_files/original_10_reference_poses/AK1_ISO_seed9009_best.sdf
  ../Workbooks/05_ZN601_10seed_comprehensive.xlsx
  results/ZN601_10seed/AK1_ZN601_pose_RMSD_matrix.csv
  results/ZN601_10seed/AK1_ZN601_scores.csv
  results/ZN601_10seed/AK1_ZN601_vs_Znfree_families.csv
  results/ZN601_10seed/AK1_ZN601_Zn_distances.csv
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed1001_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed1001_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed10101_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed10101_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed2002_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed2002_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed3003_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed3003_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed4004_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed4004_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed5005_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed5005_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed6006_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed6006_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed7007_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed7007_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed8008_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed8008_best.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed9009_allposes.sdf
  supporting_files/ZN601_10seed/original_Znfree_reference_poses/AK1_ISO_seed9009_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed1001_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed1001_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed10101_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed10101_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed2002_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed2002_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed3003_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed3003_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed4004_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed4004_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed5005_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed5005_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed6006_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed6006_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed7007_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed7007_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed8008_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed8008_best.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed9009_allposes.sdf
  supporting_files/ZN601_10seed/poses/AK1_ZN601_ISO_seed9009_best.sdf
  code/AK1_combine_panelH.py
  code/AK1_inspect_Zn_environment.py
  code/AK1_panelH1_ZN601_clean.pml
  code/AK1_panelH2_ZN641_clean.pml
  code/AK1_panelH_Zn_environment.pml
  code/AK1_panelI_ZN601_sensitivity.R
  code/AK1_panelJ_Znfree_vs_ZN601_overlay.pml
  code/AK1_ZN601_10seed_comprehensive.py
  code/AK1_ZN601_classify_pose.py
  code/AK1_ZN601_pilot_compare.py
  code/AK1_Zn_environment_to_excel.py
  code/build_AK1_FIXED_ZN601.py
  code/run_AK1_ZN601_10seeds.bat

EXECUTION SCRIPTS
The preparation helper automatically selects these files from code/portable/:
run_AK1_ZN601_10seeds.bat — Run the specified Vina searches, skip completed outputs and stop if a search fails.

RUNNING WITH THE CONFIGURED ENVIRONMENT
After preparation, work from the new structures/ directory and use:
python ..\run_analysis.py SCRIPT_NAME.py
python ..\run_analysis.py SCRIPT_NAME.R
python ..\run_analysis.py SCRIPT_NAME.pml
Docking launchers additionally require --allow-docking. Use --preview for a rapid PyMOL still-image check; full-resolution settings are the default.
