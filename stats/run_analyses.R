

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
NPROCS = 30
NITER = 200 
NITER_FULL = 200 
N_FOLDS = 30
MED_ITERS = 2000

modules_file = file.path(DIR_ICA_ME, 'modules_1.3.txt')
modules = read.table(modules_file)
module_names = seq(max(modules))


module_labels = c("SMT", #somatosensory + motor",
                  "IMOFC", # ventral frontal
                  "BGTMP", # basal ganglia + salience
                  "DLPFC", # dorsolateral prefrontal
                  "DMPFC", # dorsomedial prefrontal
                  "VIS",  # dorsomedial visual
                  "PDMN")  # default mode


to_select_list = NULL

# enough to do it for the single networks
for (j in 1:length(module_names)) to_select_list[[j]] = c(j, j);

#########################################################


# ONLY NON-DEPRESSED SAMPLE
demo.pars = demo.pars.nondep
select_medical = "epression"
DIR_PREFFIX = "mediation_data_nondep_7mods_200_1"
mycovars = c("brain_vol", "fd", "UCL", "CBU", "meica_dof", "IQcomp")

# decision acuity
var = "decAc"
controlIQ = T
source("stats/mediation.R")

mycovars = c("brain_vol", "fd", "UCL", "CBU", "meica_dof")
var = "IQcomp"
controlIQ = F
source("stats/mediation.R")

