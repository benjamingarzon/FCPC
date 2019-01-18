
# do
#

# add xvalidation noise to charts (do it 5 times)
# try with intelligence, does it make more sense?
# run with iter 20 same results?
##########################################
# add the IQ missing people before running it
###############################################
# try with depressed
##################################


# when I remove visual, motor or basal ganglia, the rediction is not mediated by age or sex
# removing visuomotor makes it not significant anymore
# visual is the only one which predicts alone, and that is not mediated by age

#1 motor, 2 ventral frontal, 3 BG + salience, 4 dorsolateral frontal, 5 dorsomedial frontal, 6 visual, 7 DMN
# best combination to predict 4 + 6

# removing DMN also important effect
# check results with depressed

# plot results
# repeat with only diag terms....

rm(list=ls())
options(warn = 0)
OUTPUT_DIR = "~/Data/NSPN/module_stats/"

fig_count = 0
setwd("~/Software/FCPC/")
#source("~/Software/DAD/VariabilityAnalysis/aux_funcs.R")
#source("aux_analyze.R")
source("stats/aux_pls.R")
source("stats/prediction_funcs.R")
source("process/organize_data.R")


# define parameters

DIR_ICA_ME = '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica'
demo.pars = demo.pars.2k

RESULTS_DIR = file.path(OUTPUT_DIR, "mediation_data_2k_7mods_single_decAc")
dir.create(RESULTS_DIR)
RESULTS_FILE = file.path(RESULTS_DIR, "data.RData")
#load(RESULTS_FILE)
#modules_file = file.path(DIR_ICA_ME, 'BrainMapPartition.txt')
#module_names = read.table('~/Data/NSPN/BrainMapNames.txt', sep = ';', header = T)
modules_file = file.path(DIR_ICA_ME, 'modules_1.3.txt')
modules = read.table(modules_file)
module_names = seq(max(modules))
#1 motor, 2 ventral frontal, 3 BG + salience, 4 dorsolateral frontal, 5 dorsomedial frontal, 6 visual, 7 DMN
module_labels = c("SMT", #somatosensory + motor",
                  "VPFC", # ventral frontal
                  "BGSAL", # basal ganglia + salience
                  "DLPFC", # dorsolateral prefrontal
                  "DMPFC", # dorsomedial prefrontal
                  "VIS",  # dorsomedial visual
                  "DMN")  # default mode

NPROCS = 30
NITER = 50 # 200
N_FOLDS = 50
MED_ITERS = 2000



mycovars = c("brain_vol", "fd", "UCL", "CBU", "meica_dof")#, "TIV", "processing")
to_select_list = NULL

# enough to do it for the single networks
for (j in 1:length(module_names)) to_select_list[[j]] = c(j, j);
     
var = "decAc"
#var = "wasi_zl_matrix_raw_score"
#var = "wasi_za_vocab_raw_score"

# do prediction for different networks
var.network.results = do_predictions_loop(DIR_ICA_ME, 
                                         demo.pars,
                                         var,
                                         mycovars,
                                         modules_file,
                                         to_select_list,
                                         NPROCS = NPROCS, 
                                         N_FOLDS = N_FOLDS, 
                                         NITER = NITER
                                         )

# do mediation for each combination tested
var.mediation.network.results = NULL

for (i in seq(length(var.network.results))){
  print(i)
  var.mediation.network.results[[i]] = do_mediation(DIR_ICA_ME, var.network.results[[i]], NSIMS = MED_ITERS)
}


NITER = 200

save.image(RESULTS_FILE)

# do whole-brain prediction
var.prediction.results = do_prediction(DIR_ICA_ME, 
                                         demo.pars, 
                                         var, 
                                         mycovars, 
                                         N_FOLDS = N_FOLDS, 
                                         NPROCS = NPROCS, 
                                         NITER = NITER)

# calculate mediation
var.mediation.results = do_mediation(DIR_ICA_ME, var.prediction.results, write_coefs = T)

save.image(RESULTS_FILE)


# check the effect of increasing the number of iterations
if(F){
  iter.results = list()
  NITERS = c(rep(1, 5), rep(20, 5), rep(50, 5), rep(100, 5))
  for (NIT in NITERS){
    print(NIT)
    iter.result =  do_prediction(DIR_ICA_ME, 
                                 demo.pars, 
                                 var, 
                                 mycovars, 
                                 N_FOLDS = N_FOLDS, 
                                 NPROCS = NPROCS, 
                                 NITER = NIT, 
                                 modules_file = modules_file,
                                 to_select = to_select_list[[21]])
    iter.results = c(iter.results, list(iter.result))
    
    iter.rho = sapply(iter.results, function(x) ifelse(!x$empty, x$cor.test$estimate, NA))
    plot(NITERS[1:length(iter.rho)], iter.rho, xlab = "NITER", ylab = "Correlation between true and predicted", ylim = c(-0.1, 0.3), pch = 20)
    
  } 
}
