
###############################################################################
# Predicting decision acuity based on functional connectivity on NSPN dataset
# benjamin.garzon@gmail.com
###############################################################################

rm(list=ls())
options(warn = 0)
OUTPUT_DIR = "~/Data/NSPN/results/"

fig_count = 0
setwd("~/Software/FCPC/")
source("stats/aux_pls.R")
source("stats/prediction_funcs.R")
#source("process/organize_data.R")

###############################################################################
# define parameters
###############################################################################
NPROCS = 5# 35
N_FOLDS = 20
NITER = 3 #200
NPERM = 5 #101
DEFAULTMAXCOMP = 10

NETMATS_FILE = 'netmats_ridge.csv'
SEPARATOR = ","

DIR_ICA_ME.BSL = '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica'
DIR_ICA_ME.FU = '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica'

DIR_PREFFIX0 = paste0("prediction_data_nondep_ICA_ridge_", N_FOLDS, "folds_", 
                      NITER, "iters_")

# load behaviour
PARAMS_DATA_FILE_NSPN = '~/Data/NSPN/params_NSPN_nondep.csv'
data.pars = read.table(PARAMS_DATA_FILE_NSPN, sep = ';')


# define modules  
modules_file = file.path('~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica', 
                         'modules_optimal.txt')
module_labels = c("SMT", # 1 sensory/somatomotorsomatosensory + motor",
                  "VIS", # 2 visual
                  "MPC", # 3 medial prefrontal cortex
                  "FPN", # 4 frontoparietal control network
                  "FPL", # 5 frontal pole
                  "OFC", # 6 orbitofrontal cortex,medial and lateral
                  "RDC", # 7 right DLPFC
                  "ATC", # 8 anterior temporal cortex (including MTL)
                  "OPC", # 9 opercular cortex
                  "PTC", # 10 posterior temporal cortex
                  "PCC", # 11 posterior cingulate cortex
                  "LDC", # 12 left DLPFC
                  "SUB", # 13 subcortical
                  "SAN")  # 14 salience network 

modules = read.table(modules_file)
module_names = seq(max(modules))

# list of modules
to_select_list = NULL

for (j in 1:max(module_names)) to_select_list[[j]] = c(j, j);
to_select_list[[max(module_names)+1]] = c(0, 0);

###############################################################################
# Correcting for demographic factors, on complete data
###############################################################################

DIR_PREFFIX = paste0(DIR_PREFFIX0, "democorrected_")

print("------------------------------")
print("Demo corrected")
print("------------------------------")

##################
# decision acuity
##################
if (T){
# baseline
testmodels = NULL
DIR_ICA_ME = DIR_ICA_ME.BSL
var = "decAc.bsl"
age_var = 'age_scan1'
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", 
             "age.2", "sex.num", "age.1.sex", "age.2.sex")

var.network.results = 
  do_predictions_loop(OUTPUT_DIR, 
                      DIR_PREFFIX, DIR_ICA_ME, 
                      data.pars, var, mycovars, 
                      modules_file, to_select_list,
                      maxcomp = DEFAULTMAXCOMP, 
                      N_FOLDS = N_FOLDS,
                      NPROCS = NPROCS, # processors to use
                      NITER = NITER, # spls bagging iteration
                      invert = F, 
                      netmats_file = NETMATS_FILE,
                      age_var = age_var,
                      NPERM = NPERM,
                      mymodels = testmodels, 
                      savecoefs = T)

# project to followup
testmodels = var.network.results
rm(var.network.results)
DIR_ICA_ME = DIR_ICA_ME.FU
var = "decAc.fu"
age_var = 'age_scan2'
mycovars = c("brain_vol", "fd", "UCL.2", "CBU.2", "meica_dof", "age.1", 
             "age.2", "sex.num", "age.1.sex", "age.2.sex")

var.network.results = 
  do_predictions_loop(OUTPUT_DIR, 
                      DIR_PREFFIX, DIR_ICA_ME, 
                      data.pars, var, mycovars, 
                      modules_file, to_select_list,
                      maxcomp = DEFAULTMAXCOMP, 
                      N_FOLDS = N_FOLDS,
                      NPROCS = NPROCS, # processors to use
                      NITER = NITER, # spls bagging iteration
                      invert = F, 
                      netmats_file = NETMATS_FILE,
                      age_var = age_var,
                      NPERM = NPERM,
                      mymodels = testmodels, 
                      savecoefs = F)

rm(testmodels, var.network.results)
}
##################
# IQ composite
##################
if (F){
# only baseline
testmodels = NULL
DIR_ICA_ME = DIR_ICA_ME.BSL
var = "IQcomp.bsl"
age_var = 'age_scan1'
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", 
             "age.2", "sex.num", "age.1.sex", "age.2.sex")
var.network.results = 
  do_predictions_loop(OUTPUT_DIR, 
                      DIR_PREFFIX, DIR_ICA_ME, 
                      data.pars, var, mycovars, 
                      modules_file, to_select_list,
                      maxcomp = DEFAULTMAXCOMP, 
                      N_FOLDS = N_FOLDS,
                      NPROCS = NPROCS, # processors to use
                      NITER = NITER, # spls bagging iteration
                      invert = F, 
                      netmats_file = NETMATS_FILE,
                      age_var = age_var,
                      NPERM = NPERM,
                      mymodels = testmodels, 
                      savecoefs = T)

# project to followup
testmodels = var.network.results
rm(var.network.results)
DIR_ICA_ME = DIR_ICA_ME.FU
var = "IQcomp.fu"
age_var = 'age_scan2'
mycovars = c("brain_vol", "fd", "UCL.2", "CBU.2", "meica_dof", "age.1", 
             "age.2", "sex.num", "age.1.sex", "age.2.sex")
var.network.results = 
  do_predictions_loop(OUTPUT_DIR, 
                      DIR_PREFFIX, DIR_ICA_ME, 
                      data.pars, var, mycovars, 
                      modules_file, to_select_list,
                      maxcomp = DEFAULTMAXCOMP, 
                      N_FOLDS = N_FOLDS,
                      NPROCS = NPROCS, # processors to use
                      NITER = NITER, # spls bagging iteration
                      invert = F, 
                      netmats_file = NETMATS_FILE,
                      age_var = age_var,
                      NPERM = NPERM,
                      mymodels = testmodels, 
                      savecoefs = F)

rm(testmodels, var.network.results)
}

#############################################################
# Repeat analyses, correcting for cognition (decAc / IQ)
#############################################################
DIR_PREFFIX = paste0(DIR_PREFFIX0, "cogcorrected_")

print("------------------------------")
print("Cognition corrected")
print("------------------------------")

##################
# decision acuity
##################

# only baseline
testmodels = NULL
if (T){
DIR_ICA_ME = DIR_ICA_ME.BSL
var = "decAc.bsl"
age_var = 'age_scan1'
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", 
             "age.2", "sex.num", "age.1.sex", "age.2.sex", "IQcomp.bsl")
var.network.results = 
  do_predictions_loop(OUTPUT_DIR, 
                      DIR_PREFFIX, DIR_ICA_ME, 
                      data.pars, var, mycovars, 
                      modules_file, to_select_list,
                      maxcomp = DEFAULTMAXCOMP, 
                      N_FOLDS = N_FOLDS,
                      NPROCS = NPROCS, # processors to use
                      NITER = NITER, # spls bagging iteration
                      invert = F, 
                      netmats_file = NETMATS_FILE,
                      age_var = age_var,
                      NPERM = NPERM,
                      mymodels = testmodels, 
                      savecoefs = F)
rm(var.network.results)
}
##################
# IQ composite
##################

# only baseline
DIR_ICA_ME = DIR_ICA_ME.BSL
var = "IQcomp.bsl"
age_var = 'age_scan1'
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", 
             "age.2", "sex.num", "age.1.sex", "age.2.sex", "decAc.bsl")
var.network.results = 
  do_predictions_loop(OUTPUT_DIR, 
                      DIR_PREFFIX, DIR_ICA_ME, 
                      data.pars, var, mycovars, 
                      modules_file, to_select_list,
                      maxcomp = DEFAULTMAXCOMP, 
                      N_FOLDS = N_FOLDS,
                      NPROCS = NPROCS, # processors to use
                      NITER = NITER, # spls bagging iteration
                      invert = F, 
                      netmats_file = NETMATS_FILE,
                      age_var = age_var,
                      NPERM = NPERM,
                      mymodels = testmodels, 
                      savecoefs = F)
rm(var.network.results)


#############################################################
# Not correcting for age
#############################################################

DIR_PREFFIX = paste0(DIR_PREFFIX0, "notagecorrected_")
data.pars = data.pars.nondep

print("------------------------------")
print("Not age corrected")
print("------------------------------")

##################
# decision acuity
##################

if (F){
  # baseline
  testmodels = NULL
  DIR_ICA_ME = DIR_ICA_ME.BSL
  var = "decAc.bsl"
  age_var = 'age_scan1'
  mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "sex.num")
  var.network.results = 
    do_predictions_loop(OUTPUT_DIR, 
                        DIR_PREFFIX, DIR_ICA_ME, 
                        data.pars, var, mycovars, 
                        modules_file, to_select_list,
                        maxcomp = DEFAULTMAXCOMP, 
                        N_FOLDS = N_FOLDS,
                        NPROCS = NPROCS, # processors to use
                        NITER = NITER, # spls bagging iteration
                        invert = F, 
                        netmats_file = NETMATS_FILE,
                        age_var = age_var,
                        NPERM = NPERM,
                        mymodels = testmodels, 
                        savecoefs = F)
  

  
}

