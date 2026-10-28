AK1–ISORHAMNETIN DOCKING: METHODS AND INSTRUCTIONS

The 8 chapters contain the analysis code, molecular inputs, saved results,
figures and reproduction instructions for Figure 8 and Figures S13–S19.


Each chapter has a README and a method-titled code TXT containing its complete
source code. The filename ends with the relevant main/SI figure
references in square brackets, for example [Fig8a_FigS14a]. Consecutive panel letters identify each
panel: Fig8bd means Figure 8b and 8d; FigS15abc means Figure S15a–c.
The executable scripts are in code/.
The chapter TXT documents contain the input/output references, reproduction
commands and source listings. They are not a single combined executable script.

PREPARE A WORKING COPY
From this folder run:
python "01_Structure_and_Ligand_Preparations/supporting_files/reproduction/prepare_project.py" --output "C:/path/to/new_working_project"

Choose a new output directory. The helper reconstructs the original folder
layout and adjusts project/executable paths in that independent copy only.
This step prepares the files, run the analysis scripts separately.
Run the selected script from the new project's structures/ directory.
Use saved docking files to repeat downstream analyses. To rerun docking,
use the docking launcher and a separate new output folder, preserving saved
results. Complete commands and scientific settings are in the chapter files.

SOFTWARE
Docking engine: AutoDock Vina 1.2.7, supplied under chapter 01's
supporting_files/software/. Preparation uses PDBFixer and Meeko 0.8.0.
The analyses use Python 3.10, NumPy, pandas, SciPy, RDKit, xlsxwriter,
gemmi and source-specific imports. Workbook-related sources may require
openpyxl. R sources declare their plotting packages, including ggplot2,
dplyr, readr, tidyr and patchwork where used. Structural figures use PyMOL;
image composition uses Pillow. See each source's import/library declarations.

EXECUTION SCRIPTS
The preparation helper selects scripts from code/portable/ where available
and uses code/ for the remaining scripts. The selected files are copied into
the working project. Their source listings are included in the chapter TXT.

ENVIRONMENT SETUP
Use an activated conda environment with Python 3.10 and install packages using
python -m pip, so package installation and execution use the same interpreter.
Meeko 0.8.0, RDKit, gemmi, SciPy, pandas and xlsxwriter are required for the
relevant docking analyses; install PDBFixer/OpenMM through conda-forge for
receptor preparation. Imageio and Pillow are needed only for image/GIF scripts.
ProDy is not imported by the supplied docking Python scripts and is not a
required installation step. Prepared receptors are supplied; when rebuilding,
repair missing atoms with fix_AK1.py before running Meeko receptor preparation.
Do not remove incomplete residues to bypass preparation errors.

You may add --check-environment to the preparation command to check the selected
Python. To make a new docking search, add --fresh-docking: the supplied poses
are retained in reference_docking/ and new searches write into an empty docking/.

Run R code through Rscript, not directly as CMD commands. Install missing R
packages from an R session, or use:
Rscript -e "install.packages(c('ggplot2','dplyr','readr','tidyr','patchwork','scales'), repos='https://cloud.r-project.org')"
If packages are installed but not found, check the library paths of the chosen
R executable with Rscript -e ".libPaths()" before reinstalling. Run PyMOL source
files through pymol -cq FILE.pml, rather than entering a .pml filename in CMD.

The 7X7S preparation uses the supplied sequence-aware alignment source and
residue mapping tables. Keep these mappings and original atom correspondence.
The blind box intentionally covers the protein and exceeds Vina's default
large-box warning threshold; retain the supplied box and exhaustiveness 128.

RUNNING THE PREPARED PROJECT
The preparation helper also places run_analysis.py in the new project root.
From that project's structures/ directory use:
python ..\run_analysis.py calc_AP5_RMSD.py
python ..\run_analysis.py AK1_plot_blind_validation.R
python ..\run_analysis.py AK1_panelC_AP5_redocking.pml

The launcher uses the active Python interpreter, configures conda DLL paths,
locates the selected R version's user library and reads R sources as UTF-8.
CSV text is exported as UTF-8. PDF plots use Cairo when available, preserving
Unicode labels. The launcher checks that PyMOL creates the expected PNG files.

To select another scientific interpreter or executable, add --python PATH,
--rscript PATH, or --pymol PATH before the source filename. An explicit R
package library may be selected with --r-library PATH. No packages are installed
automatically.

For a quick PyMOL still-image check use:
python ..\run_analysis.py --preview AK1_7X7S_selected_overlay.pml
Preview images are written separately under preview_figures/ at up to 640 pixels.
Normal execution retains the original publication-resolution render settings.
Complex transparent receptor overlays can take several minutes at full size.

Docking launchers require an explicit flag, for example:
python ..\run_analysis.py --allow-docking run_AK1_AP5_10seeds.bat
To generate new searches, first prepare a new --fresh-docking project. Without
that option the provided completed files are retained and skipped by the
portable batch sources.
