01 Structure and Ligand Preparation
How were the 1Z83 chain-A receptor, crystallographic AP5 and isorhamnetin inputs prepared?

PAPER PANELS
Methods 4.14; preparation supports Figure 8 and Figures S13–S19. Representative-pose selection supports Figure 8b; the representative structural view corresponds to Figure S13b.
Figure S13b

CONTENTS
Structure_Ligand_Preparation_and_Representative_Pose [Fig8_FigsS13-S19].txt — source code and execution instructions
code/ — Python, R, batch and PyMOL scripts
inputs/ — molecular structures and input datasets where applicable
results/ — saved docking results and numerical analyses
figures/ — figure outputs
supporting_files/ — intermediate files needed by the sources
INPUTS_AND_OUTPUTS.txt — script and file correspondence
REPRODUCTION_COMMANDS.txt — commands for running the supplied scripts

METHODS AND PARAMETERS
Key inputs include structures/AK1_1Z83.pdb, ligand/AK1_native_AP5_chainA.pdb, ligand/isorhamnetin_3D.sdf and receptors/AK1_chainA_FIXED.pdb. fix_AK1.py uses PDBFixer; make_AK1_box.py derives the targeted search box. Chain extraction, receptor preparation and ligand preparation commands are supplied in REPRODUCTION_COMMANDS.txt.

REPRODUCING THIS CHAPTER
Prepare a working project as described in the docking README.txt. The setup creates the structures/, receptors/, ligand/, docking/ and results/ directories and sets the project and Vina paths. Run each script from structures/.
Use the saved docking outputs for downstream analysis. For new searches, prepare a project with --fresh-docking and run the appropriate docking launcher.
Python scripts: python SCRIPT_NAME.py
R scripts: Rscript SCRIPT_NAME.R
PyMOL scripts: pymol -cq SCRIPT_NAME.pml
Batch scripts: run the chosen .bat from the working project's structures directory.
Each source section in the TXT is a separate script. Use the preparation helper to select the execution scripts and place them in the working project.

FILES
  inputs/ligand/AK1_native_AP5_chainA.pdb
  inputs/ligand/AP5.pdbqt
  inputs/ligand/AP5_ideal.sdf
  inputs/ligand/isorhamnetin.pdbqt
  inputs/ligand/isorhamnetin_3D.sdf
  inputs/receptors/AK1.pdbqt
  inputs/receptors/AK1_chainA.pdb
  inputs/receptors/AK1_chainA_FIXED.pdb
  inputs/structures/AK1_1Z83.pdb
  code/AK1_dominant_representative.py
  code/AK1_final_representative.pml
  code/AK1_repr9009_static.pml
  code/fix_AK1.py
  code/make_AK1_box.py

RUNNING WITH THE CONFIGURED ENVIRONMENT
After preparation, work from the new structures/ directory and use:
python ..\run_analysis.py SCRIPT_NAME.py
python ..\run_analysis.py SCRIPT_NAME.R
python ..\run_analysis.py SCRIPT_NAME.pml
Docking launchers additionally require --allow-docking. Use --preview for a rapid PyMOL still-image check; full-resolution settings are the default.
