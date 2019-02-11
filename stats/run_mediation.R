
# do
#

# add xvalidation noise to charts (do it 5 times)
# try with intelligence, does it make more sense?
# run with iter 20 same results?
##########################################
# xxxxxxadd the IQ missing people before running it
###############################################
# try with depressed
##################################


# when I remove visual, motor or basal ganglia, the rediction is not mediated by age or sex
# removing visuomotor makes it not significant anymore
# visual is the only one which predicts alone, and that is not mediated by age


rm(list=ls())
options(warn = 0)
OUTPUT_DIR = "~/Data/NSPN/module_stats/"

fig_count = 0
setwd("~/Software/FCPC/")
source("stats/aux_pls.R")
source("stats/prediction_funcs.R")
source("process/organize_data.R")


#########################################################
# define parameters

DIR_ICA_ME = '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica'
demo.pars = demo.pars.2k
NPROCS = 30
NITER = 200 # 50
NITER_FULL = 200 # 200
N_FOLDS = 30
MED_ITERS = 2000
DIR_PREFFIX = "mediation_data_2k_7mods_200_"

#modules_file = file.path(DIR_ICA_ME, 'BrainMapPartition.txt')
#module_names = read.table('~/Data/NSPN/BrainMapNames.txt', sep = ';', header = T)
modules_file = file.path(DIR_ICA_ME, 'modules_1.3.txt')
modules = read.table(modules_file)
module_names = seq(max(modules))

#1 motor, 2 ventral frontal, 3 BG + salience, 4 dorsolateral frontal, 5 dorsomedial frontal, 6 visual, 7 DMN

module_labels = c("SMT", #somatosensory + motor",
                  "IMOFC", # ventral frontal
                  "BGTMP", # basal ganglia + salience
                  "DLPFC", # dorsolateral prefrontal
                  "DMPFC", # dorsomedial prefrontal
                  "VIS",  # dorsomedial visual
                  "PDMN")  # default mode


mycovars = c("brain_vol", "fd", "UCL", "CBU", "meica_dof")#, "TIV", "processing")
to_select_list = NULL

# enough to do it for the single networks
for (j in 1:length(module_names)) to_select_list[[j]] = c(j, j);

#########################################################

# decision acuity
var = "decAc"
source("stats/mediation.R")

# IQ matrix
var = "IQmatrix"
source("stats/mediation.R")


# IQ vocabulary
var = "IQvocab"
source("stats/mediation.R")
