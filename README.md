
Scripts for connectivity analysis for "A generic decision-making ability predicts psychopathology in adolescents and young adults and is reflected in brain connectivity patterns that are distinct from general intelligence" (Michael Moutoussis and Benjamín Garzón, Sharon Neufeld, Dominik R. Bach, Francesco Rigoli, NSPN Consortium, Marc Guitart-Masip and Raymond J. Dolan)


## Data processing
Resting-state data were initially processed with MEICA (https://afni.nimh.nih.gov/pub/dist/doc/program_help/meica.py.html), Melodic (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MELODIC), and FSLnets (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets).

process/connectivity/run_meica.sh: run meica analysis on all subjects

process/connectivity/run_melodic.sh: run Melodic

process/connectivity/run_melodicpost.sh: runs Melodic and dual regression

process/compute_modules.m: compute community detection based on connectivity matrices

process/organize_data.R: merges data from different files into one matrix

process/name_ROIS.sh: name the ROIs (maps) from the ICA analysis based on HarvardOxford atlas


## Statistical analyses
stats/run_analyses.R: main analysis script

stats/mediation.R: calls do_mediation and do_prediction

stats/prediction_funcs.R: functions for prediction of cognitive scores (do_mediation and do_prediction)

stats/aux_pls.R: wrappers for PLS functions

stats/plot_results.Rmd: plot accuracy and create mediation matrices
