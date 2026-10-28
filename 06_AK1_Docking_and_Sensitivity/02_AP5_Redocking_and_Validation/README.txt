02 AP5 Redocking and Validation
Can native AP5 orientations be recovered, and how do AP5 and ISO docking scores compare?

PAPER PANELS
Figure 8a: AP5 crystallographic redocking. Figure S14a: ISO versus AP5 ten-search docking-score distributions.
Figure 8a, Figure S14a

CONTENTS
AP5_Redocking_RMSD_and_Docking_Score_Comparison [Fig8a_FigS14a].txt — source code and execution instructions
code/ — Python, R, batch and PyMOL scripts
inputs/ — molecular structures and input datasets where applicable
results/ — saved docking results and numerical analyses
figures/ — figure outputs
supporting_files/ — intermediate files needed by the sources
INPUTS_AND_OUTPUTS.txt — script and file correspondence
REPRODUCTION_COMMANDS.txt — commands for running the supplied scripts

METHODS AND PARAMETERS
calc_AP5_RMSD.py compares all saved AP5 poses with crystallographic AP5 using AP5_ideal.sdf as a connectivity template. It calculates direct-coordinate heavy-atom RMSD over chemically equivalent atom mappings, without ligand superposition. The AP5 redocking and multi-seed scoring files are retained separately. The PyMOL scripts generate the AP5 overlays.

REPRODUCING THIS CHAPTER
Prepare a working project as described in the docking README.txt. The setup creates the structures/, receptors/, ligand/, docking/ and results/ directories and sets the project and Vina paths. Run each script from structures/.
Use the saved docking outputs for downstream analysis. For new searches, prepare a project with --fresh-docking and run the appropriate docking launcher.
Python scripts: python SCRIPT_NAME.py
R scripts: Rscript SCRIPT_NAME.R
PyMOL scripts: pymol -cq SCRIPT_NAME.pml
Batch scripts: run the chosen .bat from the working project's structures directory.
Each source section in the TXT is a separate script. Use the preparation helper to select the execution scripts and place them in the working project.

FILES
  results/docking/AK1_AP5_redock_bestpose.sdf
  results/docking/AK1_AP5_redock_seed1001.pdbqt
  results/docking/AK1_AP5_redock_seed1001.sdf
  figures/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_10seed_score_comparison.pdf
  figures/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_10seed_score_comparison.png
  ../Workbooks/02_AP5_vs_ISO_10seed_score_comparison.xlsx
  figures/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_10seed_score_comparison_V2.pdf
  figures/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_10seed_score_comparison_V2.png
  figures/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_10seed_score_comparison_V3.pdf
  figures/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_10seed_score_comparison_V3.png
  results/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_all_seed_scores.csv
  results/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_seedwise_comparison.csv
  results/AK1_AP5_vs_ISO_scores/AK1_AP5_vs_ISO_summary_statistics.csv
  figures/AK1_panelC/AK1_panelC_AP5_redocking_validation.png
  figures/AK1_panelC/AK1_panelC_AP5_redocking_validation_LEGEND.png
  code/AK1_AP5_vs_ISO_score_comparison.R
  code/AK1_panelC_AP5_redocking.pml
  code/annotate_AK1_panelC.py
  code/calc_AP5_RMSD.py
  code/run_AK1_AP5_10seeds.bat

EXECUTION SCRIPTS
The preparation helper automatically selects these files from code/portable/:
run_AK1_AP5_10seeds.bat — Run the specified Vina searches, skip completed outputs and stop if a search fails.

RUNNING WITH THE CONFIGURED ENVIRONMENT
After preparation, work from the new structures/ directory and use:
python ..\run_analysis.py SCRIPT_NAME.py
python ..\run_analysis.py SCRIPT_NAME.R
python ..\run_analysis.py SCRIPT_NAME.pml
Docking launchers additionally require --allow-docking. Use --preview for a rapid PyMOL still-image check; full-resolution settings are the default.
