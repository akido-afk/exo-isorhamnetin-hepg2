07 Receptor Ensemble Sensitivity
How does docking change across selected open 7X7S NMR conformers aligned to closed 1Z83?

PAPER PANELS
Figure 8e–g; Figure S18c–g; Figure S19a,b. These panels cover receptor alignment/opening, docking scores, AP5-region localization, mapped contacts, closed-state footprint retention and Jaccard similarity.
Figure 8e, Figure 8f, Figure 8g, Figure S18c, Figure S18d, Figure S18e, Figure S18f, Figure S18g, Figure S19a, Figure S19b

CONTENTS
Receptor_Ensemble_Docking_and_Contacts [Fig8efg_FigS18cdefg_FigS19ab].txt — source code and execution instructions
code/ — Python, R, batch and PyMOL scripts
inputs/ — molecular structures and input datasets where applicable
results/ — saved docking results and numerical analyses
figures/ — figure outputs
supporting_files/ — intermediate files needed by the sources
INPUTS_AND_OUTPUTS.txt — script and file correspondence
REPRODUCTION_COMMANDS.txt — commands for running the supplied scripts

METHODS AND PARAMETERS
The input files include the full 20-model AK1_7X7S_NMR.pdb and selected models 20 (least-open), 8 (ensemble medoid) and 1 (most-open). AK1_inspect_7X7S_ensemble.py selects models using C-alpha distances and pairwise RMSD. AK1_prepare_7X7S_selected_models.py aligns sequence-mapped core residues to 1Z83. AK1_fix_7X7S_selected_receptors.py, the batch docking launcher, ensemble analysis and sequence-mapped contact extraction are separate stages. Preserve residue mapping tables.

REPRODUCING THIS CHAPTER
Prepare a working project as described in the docking README.txt. The setup creates the structures/, receptors/, ligand/, docking/ and results/ directories and sets the project and Vina paths. Run each script from structures/.
Use the saved docking outputs for downstream analysis. For new searches, prepare a project with --fresh-docking and run the appropriate docking launcher.
Python scripts: python SCRIPT_NAME.py
R scripts: Rscript SCRIPT_NAME.R
PyMOL scripts: pymol -cq SCRIPT_NAME.pml
Batch scripts: run the chosen .bat from the working project's structures directory.
Each source section in the TXT is a separate script. Use the preparation helper to select the execution scripts and place them in the working project.

FILES
  results/docking/AK1_7X7S_model1_most_open_ISO_seed1001.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed10101.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed2002.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed3003.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed4004.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed5005.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed6006.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed7007.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed8008.pdbqt
  results/docking/AK1_7X7S_model1_most_open_ISO_seed9009.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed1001.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed10101.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed2002.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed3003.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed4004.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed5005.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed6006.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed7007.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed8008.pdbqt
  results/docking/AK1_7X7S_model20_least_open_ISO_seed9009.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed1001.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed10101.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed2002.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed3003.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed4004.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed5005.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed6006.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed7007.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed8008.pdbqt
  results/docking/AK1_7X7S_model8_ensemble_medoid_ISO_seed9009.pdbqt
  inputs/receptors/AK1_7X7S_model1_most_open.pdbqt
  inputs/receptors/AK1_7X7S_model1_most_open_FIXED.pdb
  inputs/receptors/AK1_7X7S_model20_least_open.pdbqt
  inputs/receptors/AK1_7X7S_model20_least_open_FIXED.pdb
  inputs/receptors/AK1_7X7S_model8_ensemble_medoid.pdbqt
  inputs/receptors/AK1_7X7S_model8_ensemble_medoid_FIXED.pdb
  results/7X7S_inspection/AK1_7X7S_candidate_conformers.csv
  ../Workbooks/07_7X7S_ENSEMBLE_INSPECTION.xlsx
  results/7X7S_inspection/AK1_7X7S_model_summary.csv
  results/7X7S_inspection/AK1_7X7S_pairwise_CA_RMSD.csv
  figures/7X7S_selected_models/AK1_1Z83_vs_7X7S_selected_overlay.png
  results/7X7S_selected_models/AK1_7X7S_pocket_residue_mapping.csv
  results/7X7S_selected_models/AK1_7X7S_selected_model_alignment_summary.csv
  ../Workbooks/07_7X7S_SELECTED_MODEL_Prep_V2.xlsx
  results/7X7S_selected_models/AK1_7X7S_to_1Z83_sequence_mapping.csv
  figures/7X7S_selected_models/AK1_7X7S_V2_closed_vs_least_open.png
  figures/7X7S_selected_models/AK1_7X7S_V2_closed_vs_medoid.png
  figures/7X7S_selected_models/AK1_7X7S_V2_closed_vs_most_open.png
  results/docking_analysis_v1/AK1_1Z83_closed_pairwise_RMSD.csv
  results/docking_analysis_v1/AK1_7X7S_AP5_region_metrics.csv
  results/docking_analysis_v1/AK1_7X7S_cross_family_RMSD.csv
  ../Workbooks/07_7X7S_ENSEMBLE_DOCKING_Analysis.xlsx
  results/docking_analysis_v1/AK1_7X7S_ensemble_seed_results.csv
  results/docking_analysis_v1/AK1_7X7S_ensemble_summary.csv
  results/docking_analysis_v1/AK1_7X7S_least_open_pairwise_RMSD.csv
  results/docking_analysis_v1/AK1_7X7S_medoid_pairwise_RMSD.csv
  results/docking_analysis_v1/AK1_7X7S_most_open_pairwise_RMSD.csv
  results/mapped_contact_fingerprint_v1/AK1_ensemble_all_pose_contacts_mapped_to_1Z83.csv
  results/mapped_contact_fingerprint_v1/AK1_ensemble_contact_comparison_matrix.csv
  results/mapped_contact_fingerprint_v1/AK1_ensemble_contact_overlap_vs_closed.csv
  results/mapped_contact_fingerprint_v1/AK1_ensemble_contact_recurrence_mapped_to_1Z83.csv
  ../Workbooks/07_ENSEMBLE_MAPPED_CONTACT_Contacts.xlsx
  results/mapped_contact_fingerprint_v1/AK1_ensemble_recurrent_ANY_contacts_80pct.csv
  results/mapped_contact_fingerprint_v1/AK1_ensemble_recurrent_POLAR_contacts_80pct.csv
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_ANY_contact_heatmap.pdf
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_ANY_contact_heatmap.png
  results/mapped_contact_fingerprint_v1/figures/AK1_ensemble_ANY_contact_heatmap_plot_data.csv
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_closed_fingerprint_retention.pdf
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_closed_fingerprint_retention.png
  results/mapped_contact_fingerprint_v1/figures/AK1_ensemble_closed_fingerprint_retention_plot_data.csv
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_Jaccard_similarity.pdf
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_Jaccard_similarity.png
  results/mapped_contact_fingerprint_v1/figures/AK1_ensemble_Jaccard_similarity_plot_data.csv
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_POLAR_contact_heatmap.pdf
  figures/mapped_contact_fingerprint_v1/figures/AK1_ensemble_POLAR_contact_heatmap.png
  results/mapped_contact_fingerprint_v1/figures/AK1_ensemble_POLAR_contact_heatmap_plot_data.csv
  results/quantitative_figure_v2/AK1_ensemble_quantitative_figure_summary.csv
  figures/quantitative_figure_v2/figures/AK1_ensemble_panelA_opening_metric.pdf
  figures/quantitative_figure_v2/figures/AK1_ensemble_panelA_opening_metric.png
  figures/quantitative_figure_v2/figures/AK1_ensemble_panelB_Vina_scores.pdf
  figures/quantitative_figure_v2/figures/AK1_ensemble_panelB_Vina_scores.png
  figures/quantitative_figure_v2/figures/AK1_ensemble_panelC_AP5_localization.pdf
  figures/quantitative_figure_v2/figures/AK1_ensemble_panelC_AP5_localization.png
  figures/quantitative_figure_v2/figures/AK1_ensemble_quantitative_ABC_combined.pdf
  figures/quantitative_figure_v2/figures/AK1_ensemble_quantitative_ABC_combined.png
  ../Workbooks/07_7X7S_RECEPTOR_Prep.xlsx
  results/receptor_preparation/AK1_7X7S_receptor_preparation_summary.csv
  inputs/structures/AK1_7X7S_model1_most_open_ALIGNED_1Z83.pdb
  inputs/structures/AK1_7X7S_model1_most_open_RAW.pdb
  inputs/structures/AK1_7X7S_model20_least_open_ALIGNED_1Z83.pdb
  inputs/structures/AK1_7X7S_model20_least_open_RAW.pdb
  inputs/structures/AK1_7X7S_model8_ensemble_medoid_ALIGNED_1Z83.pdb
  inputs/structures/AK1_7X7S_model8_ensemble_medoid_RAW.pdb
  inputs/structures/AK1_7X7S_NMR.pdb
  code/AK1_7X7S_pairwise_overlays_v2.pml
  code/AK1_7X7S_selected_overlay.pml
  code/AK1_analyze_7X7S_ensemble_docking.py
  code/AK1_extract_ensemble_mapped_contacts.py
  code/AK1_fix_7X7S_selected_receptors.py
  code/AK1_inspect_7X7S_ensemble.py
  code/AK1_plot_ensemble_contact_fingerprint.R
  code/AK1_plot_ensemble_Jaccard_similarity.R
  code/AK1_plot_ensemble_quantitative_ABC.R
  code/AK1_prepare_7X7S_selected_models.py
  code/AK1_run_7X7S_ensemble_docking.bat

EXECUTION SCRIPTS
The preparation helper automatically selects these files from code/portable/:
AK1_run_7X7S_ensemble_docking.bat — Run the specified Vina searches, skip completed outputs and stop if a search fails.

RUNNING WITH THE CONFIGURED ENVIRONMENT
After preparation, work from the new structures/ directory and use:
python ..\run_analysis.py SCRIPT_NAME.py
python ..\run_analysis.py SCRIPT_NAME.R
python ..\run_analysis.py SCRIPT_NAME.pml
Docking launchers additionally require --allow-docking. Use --preview for a rapid PyMOL still-image check; full-resolution settings are the default.
