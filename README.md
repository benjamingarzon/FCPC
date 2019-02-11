#Functional connectome prediction of cognition (FCPC)

Files:

process/compute_modules.m: compute community detection based on connectivity matrices

process/organize_data.R: merges data from different files into one matrix

process/name_ROIS.sh: name the ROIs (maps) from the ICA analysis based on HarvardOxford atlas

stats/mediation_report.R: create a report from output of mediation.R

stats/mediation.R: calls do_mediation and do_prediction

stats/prediction_funcs.R: functions for prediction of cognitive scores (do_mediation and do_prediction)

stats/aux_pls.R: wrappers for PLS functions

stats/slices_summary_2: allows to overlay atlas on ICA group maps for better identification of edges. Requires overlay.py and needs to be run before plot_weights

stats/plot_weights.m: plot pairs of ROIs with highest weights

stats/plot_results.Rmd: plot accuracy and create mediation matrices

mri_vol2vol --mov ../../ICA_ME/nspn_ME/ica200.gica/modules_1.3_smooth.nii.gz --targ /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz --regheader --o modules.nii.gz
#slices_summary modules.nii.gz 0.01 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz modules.sum

~/Software/FCPC/stats/slices_summary_2 modules.nii.gz 0.01 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz modules.sum 1 0.5 -all
# run slices dir for small ones with more contrast