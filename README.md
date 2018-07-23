#Functional connectome prediction of cognition (FCPC)

Files:

process/compute_modules.m: compute community detection based on connectivity matrices

process/organize_data.R: merges data from different files into one matrix

process/name_ROIS.sh: name the ROIs (maps) from the ICA analysis based on HarvardOxford atlas


stats/mediation_report.R: create a report from output of mediation.R

stats/mediation.R: calls do_mediation and do_prediction

stats/prediction_funcs.R: functions for prediction of cognitive scores (do_mediation and do_prediction)

stats/aux_pls.R: wrappers for PLS functions

stats/plot_weights.m: plot pairs of ROIs with highest weights


