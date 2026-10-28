03 ISO Targeted Docking and Convergence
How reproducible are the neutral-isorhamnetin docking poses across ten independent seeds?

PAPER PANELS
Figure 8c: ten-pose pairwise RMSD and dominant/alternative families. Figure S13a: independently seeded docking scores. Figure S13c: ISO placement relative to the nucleotide-binding region.
Figure 8c, Figure S13a, Figure S13c

CONTENTS
Targeted_Docking_Scores_Clustering_and_Pose_Convergence [Fig8c_FigS13ac].txt — source code and execution instructions
code/ — Python, R, batch and PyMOL scripts
inputs/ — molecular structures and input datasets where applicable
results/ — saved docking results and numerical analyses
figures/ — figure outputs
supporting_files/ — intermediate files needed by the sources
INPUTS_AND_OUTPUTS.txt — script and file correspondence
REPRODUCTION_COMMANDS.txt — commands for running the supplied scripts

METHODS AND PARAMETERS
run_AK1_multiseed.py uses receptors/AK1.pdbqt and ligand/isorhamnetin.pdbqt; center (25.044,44.571,-3.708) Å, size (23.687,23.364,30.062) Å, exhaustiveness 16, ten modes, seeds 1001,2002,3003,4004,5005,6006,7007,8008,9009,10101. AK1_pose_convergence.py calculates direct atom-order heavy-atom RMSDs of the ten rank-1 poses, without fitting or symmetry correction. AK1_cluster_poses.py uses average linkage at 2 and 3 Å. AK1_dominant_representative.py fixes the eight dominant seeds and calculates their medoid. build_AK1_comprehensive_excel.py calculates scores/clusters but embeds AP5_REDocking_RMSD=1.776 as metadata.

REPRODUCING THIS CHAPTER
Prepare a working project as described in the docking README.txt. The setup creates the structures/, receptors/, ligand/, docking/ and results/ directories and sets the project and Vina paths. Run each script from structures/.
Use the saved docking outputs for downstream analysis. For new searches, prepare a project with --fresh-docking and run the appropriate docking launcher.
Python scripts: python SCRIPT_NAME.py
R scripts: Rscript SCRIPT_NAME.R
PyMOL scripts: pymol -cq SCRIPT_NAME.pml
Batch scripts: run the chosen .bat from the working project's structures directory.
Each source section in the TXT is a separate script. Use the preparation helper to select the execution scripts and place them in the working project.

FILES
  results/docking/AK1_AP5_seed1001.pdbqt
  results/docking/AK1_AP5_seed10101.pdbqt
  results/docking/AK1_AP5_seed2002.pdbqt
  results/docking/AK1_AP5_seed3003.pdbqt
  results/docking/AK1_AP5_seed4004.pdbqt
  results/docking/AK1_AP5_seed5005.pdbqt
  results/docking/AK1_AP5_seed6006.pdbqt
  results/docking/AK1_AP5_seed7007.pdbqt
  results/docking/AK1_AP5_seed8008.pdbqt
  results/docking/AK1_AP5_seed9009.pdbqt
  results/AK1_ISO_bestpose.sdf
  results/docking/AK1_ISO_seed1001.pdbqt
  results/docking/AK1_ISO_seed1001.sdf
  results/docking/AK1_ISO_seed10101.pdbqt
  results/docking/AK1_ISO_seed2002.pdbqt
  results/docking/AK1_ISO_seed3003.pdbqt
  results/docking/AK1_ISO_seed4004.pdbqt
  results/docking/AK1_ISO_seed5005.pdbqt
  results/docking/AK1_ISO_seed6006.pdbqt
  results/docking/AK1_ISO_seed7007.pdbqt
  results/docking/AK1_ISO_seed8008.pdbqt
  results/docking/AK1_ISO_seed9009.pdbqt
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed1001_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed1001_best.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed10101_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed10101_best.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed2002_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed2002_best.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed3003_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed3003_best.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed4004_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed4004_best.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed5005_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed5005_best.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed7007_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed7007_best.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed9009_ALL.sdf
  supporting_files/AK1_final_visuals/poses/AK1_ISO_seed9009_best.sdf
  figures/AK1_ISO_AP5_overlay.png
  ../Workbooks/03_ISO_Full_DOCKING_RESULTS.xlsx
  results/AK1_ISO_pose_RMSD_matrix.csv
  figures/AK1_ISO_rotation.gif
  figures/AK1_ISO_rotation_clean.gif
  figures/AK1_panelD/AK1_panelD_multiseed_reproducibility.pdf
  figures/AK1_panelD/AK1_panelD_multiseed_reproducibility.png
  figures/AK1_panelD/AK1_panelD_multiseed_reproducibility_V2.pdf
  figures/AK1_panelD/AK1_panelD_multiseed_reproducibility_V2.png
  figures/AK1_panelE/AK1_panelE_pose_RMSD_heatmap.pdf
  figures/AK1_panelE/AK1_panelE_pose_RMSD_heatmap.png
  figures/AK1_panelE/AK1_panelE_pose_RMSD_heatmap_CLEAN.pdf
  figures/AK1_panelE/AK1_panelE_pose_RMSD_heatmap_CLEAN.png
  figures/AK1_panelE/AK1_panelE_pose_RMSD_heatmap_EXPLAINED.pdf
  figures/AK1_panelE/AK1_panelE_pose_RMSD_heatmap_EXPLAINED.png
  figures/AK1_panelE_clean_v2.png
  figures/AK1_panelE_pose_RMSD_heatmap_clean.png
  figures/AK1_repr9009/AK1_ISO_representative_seed9009_static.png
  supporting_files/AK1_repr9009/AK1_ISO_seed9009_allposes.sdf
  supporting_files/AK1_repr9009/AK1_ISO_seed9009_best.sdf
  supporting_files/AK1_rotate/frame_0001.png
  supporting_files/AK1_rotate/frame_0002.png
  supporting_files/AK1_rotate/frame_0003.png
  supporting_files/AK1_rotate/frame_0004.png
  supporting_files/AK1_rotate/frame_0005.png
  supporting_files/AK1_rotate/frame_0006.png
  supporting_files/AK1_rotate/frame_0007.png
  supporting_files/AK1_rotate/frame_0008.png
  supporting_files/AK1_rotate/frame_0009.png
  supporting_files/AK1_rotate/frame_0010.png
  supporting_files/AK1_rotate/frame_0011.png
  supporting_files/AK1_rotate/frame_0012.png
  supporting_files/AK1_rotate/frame_0013.png
  supporting_files/AK1_rotate/frame_0014.png
  supporting_files/AK1_rotate/frame_0015.png
  supporting_files/AK1_rotate/frame_0016.png
  supporting_files/AK1_rotate/frame_0017.png
  supporting_files/AK1_rotate/frame_0018.png
  supporting_files/AK1_rotate/frame_0019.png
  supporting_files/AK1_rotate/frame_0020.png
  supporting_files/AK1_rotate/frame_0021.png
  supporting_files/AK1_rotate/frame_0022.png
  supporting_files/AK1_rotate/frame_0023.png
  supporting_files/AK1_rotate/frame_0024.png
  supporting_files/AK1_rotate/frame_0025.png
  supporting_files/AK1_rotate/frame_0026.png
  supporting_files/AK1_rotate/frame_0027.png
  supporting_files/AK1_rotate/frame_0028.png
  supporting_files/AK1_rotate/frame_0029.png
  supporting_files/AK1_rotate/frame_0030.png
  supporting_files/AK1_rotate/frame_0031.png
  supporting_files/AK1_rotate/frame_0032.png
  supporting_files/AK1_rotate/frame_0033.png
  supporting_files/AK1_rotate/frame_0034.png
  supporting_files/AK1_rotate/frame_0035.png
  supporting_files/AK1_rotate/frame_0036.png
  supporting_files/AK1_rotate/frame_0037.png
  supporting_files/AK1_rotate/frame_0038.png
  supporting_files/AK1_rotate/frame_0039.png
  supporting_files/AK1_rotate/frame_0040.png
  supporting_files/AK1_rotate/frame_0041.png
  supporting_files/AK1_rotate/frame_0042.png
  supporting_files/AK1_rotate/frame_0043.png
  supporting_files/AK1_rotate/frame_0044.png
  supporting_files/AK1_rotate/frame_0045.png
  supporting_files/AK1_rotate/frame_0046.png
  supporting_files/AK1_rotate/frame_0047.png
  supporting_files/AK1_rotate/frame_0048.png
  supporting_files/AK1_rotate/frame_0049.png
  supporting_files/AK1_rotate/frame_0050.png
  supporting_files/AK1_rotate/frame_0051.png
  supporting_files/AK1_rotate/frame_0052.png
  supporting_files/AK1_rotate/frame_0053.png
  supporting_files/AK1_rotate/frame_0054.png
  supporting_files/AK1_rotate/frame_0055.png
  supporting_files/AK1_rotate/frame_0056.png
  supporting_files/AK1_rotate/frame_0057.png
  supporting_files/AK1_rotate/frame_0058.png
  supporting_files/AK1_rotate/frame_0059.png
  supporting_files/AK1_rotate/frame_0060.png
  supporting_files/AK1_rotate/frame_0061.png
  supporting_files/AK1_rotate/frame_0062.png
  supporting_files/AK1_rotate/frame_0063.png
  supporting_files/AK1_rotate/frame_0064.png
  supporting_files/AK1_rotate/frame_0065.png
  supporting_files/AK1_rotate/frame_0066.png
  supporting_files/AK1_rotate/frame_0067.png
  supporting_files/AK1_rotate/frame_0068.png
  supporting_files/AK1_rotate/frame_0069.png
  supporting_files/AK1_rotate/frame_0070.png
  supporting_files/AK1_rotate/frame_0071.png
  supporting_files/AK1_rotate/frame_0072.png
  supporting_files/AK1_rotate/frame_0073.png
  supporting_files/AK1_rotate/frame_0074.png
  supporting_files/AK1_rotate/frame_0075.png
  supporting_files/AK1_rotate/frame_0076.png
  supporting_files/AK1_rotate/frame_0077.png
  supporting_files/AK1_rotate/frame_0078.png
  supporting_files/AK1_rotate/frame_0079.png
  supporting_files/AK1_rotate/frame_0080.png
  supporting_files/AK1_rotate/frame_0081.png
  supporting_files/AK1_rotate/frame_0082.png
  supporting_files/AK1_rotate/frame_0083.png
  supporting_files/AK1_rotate/frame_0084.png
  supporting_files/AK1_rotate/frame_0085.png
  supporting_files/AK1_rotate/frame_0086.png
  supporting_files/AK1_rotate/frame_0087.png
  supporting_files/AK1_rotate/frame_0088.png
  supporting_files/AK1_rotate/frame_0089.png
  supporting_files/AK1_rotate/frame_0090.png
  supporting_files/AK1_rotate/frame_0091.png
  supporting_files/AK1_rotate/frame_0092.png
  supporting_files/AK1_rotate/frame_0093.png
  supporting_files/AK1_rotate/frame_0094.png
  supporting_files/AK1_rotate/frame_0095.png
  supporting_files/AK1_rotate/frame_0096.png
  supporting_files/AK1_rotate/frame_0097.png
  supporting_files/AK1_rotate/frame_0098.png
  supporting_files/AK1_rotate/frame_0099.png
  supporting_files/AK1_rotate/frame_0100.png
  supporting_files/AK1_rotate/frame_0101.png
  supporting_files/AK1_rotate/frame_0102.png
  supporting_files/AK1_rotate/frame_0103.png
  supporting_files/AK1_rotate/frame_0104.png
  supporting_files/AK1_rotate/frame_0105.png
  supporting_files/AK1_rotate/frame_0106.png
  supporting_files/AK1_rotate/frame_0107.png
  supporting_files/AK1_rotate/frame_0108.png
  supporting_files/AK1_rotate/frame_0109.png
  supporting_files/AK1_rotate/frame_0110.png
  supporting_files/AK1_rotate/frame_0111.png
  supporting_files/AK1_rotate/frame_0112.png
  supporting_files/AK1_rotate/frame_0113.png
  supporting_files/AK1_rotate/frame_0114.png
  supporting_files/AK1_rotate/frame_0115.png
  supporting_files/AK1_rotate/frame_0116.png
  supporting_files/AK1_rotate/frame_0117.png
  supporting_files/AK1_rotate/frame_0118.png
  supporting_files/AK1_rotate/frame_0119.png
  supporting_files/AK1_rotate/frame_0120.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0001.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0002.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0003.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0004.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0005.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0006.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0007.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0008.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0009.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0010.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0011.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0012.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0013.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0014.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0015.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0016.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0017.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0018.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0019.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0020.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0021.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0022.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0023.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0024.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0025.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0026.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0027.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0028.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0029.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0030.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0031.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0032.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0033.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0034.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0035.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0036.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0037.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0038.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0039.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0040.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0041.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0042.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0043.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0044.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0045.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0046.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0047.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0048.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0049.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0050.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0051.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0052.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0053.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0054.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0055.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0056.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0057.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0058.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0059.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0060.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0061.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0062.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0063.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0064.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0065.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0066.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0067.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0068.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0069.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0070.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0071.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0072.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0073.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0074.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0075.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0076.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0077.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0078.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0079.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0080.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0081.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0082.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0083.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0084.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0085.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0086.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0087.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0088.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0089.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0090.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0091.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0092.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0093.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0094.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0095.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0096.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0097.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0098.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0099.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0100.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0101.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0102.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0103.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0104.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0105.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0106.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0107.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0108.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0109.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0110.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0111.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0112.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0113.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0114.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0115.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0116.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0117.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0118.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0119.png
  supporting_files/AK1_rotate seed 1001 frames/frame_0120.png
  figures/AK1_rotate_ISO_orange/AK1_ISO_AP5_overlay_orange.png
  figures/AK1_rotate_ISO_orange/AK1_ISO_orange_rotation_clean.gif
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0001.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0002.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0003.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0004.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0005.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0006.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0007.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0008.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0009.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0010.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0011.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0012.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0013.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0014.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0015.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0016.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0017.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0018.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0019.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0020.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0021.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0022.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0023.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0024.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0025.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0026.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0027.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0028.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0029.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0030.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0031.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0032.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0033.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0034.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0035.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0036.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0037.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0038.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0039.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0040.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0041.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0042.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0043.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0044.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0045.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0046.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0047.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0048.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0049.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0050.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0051.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0052.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0053.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0054.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0055.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0056.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0057.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0058.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0059.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0060.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0061.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0062.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0063.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0064.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0065.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0066.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0067.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0068.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0069.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0070.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0071.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0072.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0073.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0074.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0075.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0076.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0077.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0078.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0079.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0080.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0081.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0082.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0083.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0084.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0085.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0086.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0087.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0088.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0089.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0090.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0091.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0092.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0093.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0094.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0095.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0096.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0097.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0098.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0099.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0100.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0101.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0102.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0103.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0104.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0105.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0106.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0107.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0108.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0109.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0110.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0111.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0112.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0113.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0114.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0115.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0116.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0117.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0118.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0119.png
  supporting_files/AK1_rotate_ISO_orange/frames/frame_0120.png
  figures/AK1_rotate_v2.gif
  ../Workbooks/03_ISO_AK_Docking_Master_Results.xlsx
  code/AK1_cluster_poses.py
  code/AK1_cluster_scores.py
  code/AK1_panelD_multiseed.R
  code/AK1_panelD_multiseed_v2.R
  code/AK1_panelE_clean_for_paper.R
  code/AK1_panelE_clean_v2.R
  code/AK1_panelE_heatmap_clean.R
  code/AK1_panelE_paper_final.R
  code/AK1_panelE_RMSD_heatmap.R
  code/AK1_panelE_RMSD_heatmap_explained.R
  code/AK1_pose_convergence.py
  code/AK1_rotate_ISO_orange.pml
  code/AK1_rotate_orange_same_style.pml
  code/build_AK1_comprehensive_excel.py
  code/build_docking_results.py
  code/extract_AK1_seed9009_best.py
  code/extract_best_pose.py
  code/make_AK1_gif.py
  code/make_AK1_gif_clean.py
  code/make_AK1_gif_v2.py
  code/make_AK1_orange_gif_clean.py
  code/make_AK1_visuals.pml
  code/prepare_AK1_dominant_visuals.py
  code/run_AK1_multiseed.py

EXECUTION SCRIPTS
The preparation helper automatically selects these files from code/portable/:
AK1_rotate_ISO_orange.pml — Load the representative seed9009 pose and render the structural view.

RUNNING WITH THE CONFIGURED ENVIRONMENT
After preparation, work from the new structures/ directory and use:
python ..\run_analysis.py SCRIPT_NAME.py
python ..\run_analysis.py SCRIPT_NAME.R
python ..\run_analysis.py SCRIPT_NAME.pml
Docking launchers additionally require --allow-docking. Use --preview for a rapid PyMOL still-image check; full-resolution settings are the default.

ADDITIONAL EXECUTION SCRIPTS
build_docking_results.py — Build the primary docking workbook from the neutral 1Z83 AK1_ISO_seed files. Zinc and receptor-ensemble results are analyzed in their respective chapters.
