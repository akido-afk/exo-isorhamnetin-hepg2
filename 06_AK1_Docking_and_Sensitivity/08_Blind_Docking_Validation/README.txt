08 Blind Docking Validation
Does a whole-protein search recover the targeted nucleotide-region pose families?

PAPER PANELS
Figure 8h–j; Figure S19c–e. These panels show blind mode-1 RMSD, targeted-versus-blind scores and targeted-family recovery across mode-1 and all retained poses.
Figure 8h, Figure 8i, Figure 8j, Figure S19c, Figure S19d, Figure S19e

CONTENTS
Whole_Protein_Blind_Docking_and_Targeted_Pose_Recovery [Fig8hij_FigS19cde].txt — source code and execution instructions
code/ — Python, R, batch and PyMOL scripts
inputs/ — molecular structures and input datasets where applicable
results/ — saved docking results and numerical analyses
figures/ — figure outputs
supporting_files/ — intermediate files needed by the sources
INPUTS_AND_OUTPUTS.txt — script and file correspondence
REPRODUCTION_COMMANDS.txt — commands for running the supplied scripts

METHODS AND PARAMETERS
AK1_calculate_blind_docking_box.py derives a whole-protein box with 6 Å padding per side. Blind-search settings are center (19.739,45.925,1.072) Å and size (57.316,53.181,53.673) Å; inspect the batch command for exhaustiveness 128, 20 modes and energy_range 10. Pilot seed1001 and full ten-seed analyses are distinct. AK1_analyze_blind_10seed_full.py compares saved modes with both targeted reference families.

REPRODUCING THIS CHAPTER
Prepare a working project as described in the docking README.txt. The setup creates the structures/, receptors/, ligand/, docking/ and results/ directories and sets the project and Vina paths. Run each script from structures/.
Use the saved docking outputs for downstream analysis. For new searches, prepare a project with --fresh-docking and run the appropriate docking launcher.
Python scripts: python SCRIPT_NAME.py
R scripts: Rscript SCRIPT_NAME.R
PyMOL scripts: pymol -cq SCRIPT_NAME.pml
Batch scripts: run the chosen .bat from the working project's structures directory.
Each source section in the TXT is a separate script. Use the preparation helper to select the execution scripts and place them in the working project.

FILES
  results/docking/AK1_ISO_BLIND_seed1001.pdbqt
  results/docking/AK1_ISO_BLIND_seed10101.pdbqt
  results/docking/AK1_ISO_BLIND_seed2002.pdbqt
  results/docking/AK1_ISO_BLIND_seed3003.pdbqt
  results/docking/AK1_ISO_BLIND_seed4004.pdbqt
  results/docking/AK1_ISO_BLIND_seed5005.pdbqt
  results/docking/AK1_ISO_BLIND_seed6006.pdbqt
  results/docking/AK1_ISO_BLIND_seed7007.pdbqt
  results/docking/AK1_ISO_BLIND_seed8008.pdbqt
  results/docking/AK1_ISO_BLIND_seed9009.pdbqt
  supporting_files/AK1_blind_box.txt
  ../Workbooks/08_BLIND_10SEED_Full_Analysis.xlsx
  results/full_10seed_analysis_v1/AK1_BLIND_all_poses.csv
  results/full_10seed_analysis_v1/AK1_BLIND_overall_summary.csv
  results/full_10seed_analysis_v1/AK1_BLIND_seed_summary.csv
  figures/full_10seed_analysis_v1/figures_v1/AK1_blind_panelA_mode1_RMSD.pdf
  figures/full_10seed_analysis_v1/figures_v1/AK1_blind_panelA_mode1_RMSD.png
  figures/full_10seed_analysis_v1/figures_v1/AK1_blind_panelB_targeted_vs_blind_scores.pdf
  figures/full_10seed_analysis_v1/figures_v1/AK1_blind_panelB_targeted_vs_blind_scores.png
  figures/full_10seed_analysis_v1/figures_v1/AK1_blind_panelC_recovery_frequency.pdf
  figures/full_10seed_analysis_v1/figures_v1/AK1_blind_panelC_recovery_frequency.png
  figures/full_10seed_analysis_v1/figures_v1/AK1_BLIND_VALIDATION_ABC.pdf
  figures/full_10seed_analysis_v1/figures_v1/AK1_BLIND_VALIDATION_ABC.png
  ../Workbooks/08_BLIND_VALIDATION_FIGURE_DATA.xlsx
  results/pilot_seed1001/AK1_BLIND_seed1001_all20_modes_analysis.csv
  ../Workbooks/08_BLIND_seed1001_PILOT_Analysis.xlsx
  code/AK1_analyze_blind_10seed_full.py
  code/AK1_analyze_blind_pilot.py
  code/AK1_calculate_blind_docking_box.py
  code/AK1_plot_blind_validation.R
  code/AK1_run_blind_remaining9.bat

RUNNING WITH THE CONFIGURED ENVIRONMENT
After preparation, work from the new structures/ directory and use:
python ..\run_analysis.py SCRIPT_NAME.py
python ..\run_analysis.py SCRIPT_NAME.R
python ..\run_analysis.py SCRIPT_NAME.pml
Docking launchers additionally require --allow-docking. Use --preview for a rapid PyMOL still-image check; full-resolution settings are the default.

ADDITIONAL EXECUTION SCRIPTS
AK1_run_blind_remaining9.bat — Run the nine remaining blind-search seeds with exhaustiveness 128, 20 modes and energy range 10; skip completed outputs.
