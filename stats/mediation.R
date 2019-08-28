

RESULTS_DIR = file.path(OUTPUT_DIR, paste0(DIR_PREFFIX, var))
dir.create(RESULTS_DIR)
RESULTS_FILE = file.path(RESULTS_DIR, "data.RData")


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

# do whole-brain prediction
var.prediction.results = do_prediction(DIR_ICA_ME, 
                                       demo.pars, 
                                       var, 
                                       mycovars, 
                                       N_FOLDS = N_FOLDS, 
                                       NPROCS = NPROCS, 
                                       NITER = NITER_FULL)

save.image(RESULTS_FILE)

# calculate mediations

var.mediation.network.results = NULL

for (i in seq(length(var.network.results))){
  print(i)
  var.mediation.network.results[[i]] = do_mediation(DIR_ICA_ME, 
                                                    var.network.results[[i]], 
                                                    NSIMS = MED_ITERS,
                                                    RESULTS_DIR = RESULTS_DIR, 
                                                    controlIQ = F)
}

var.mediation.results = do_mediation(DIR_ICA_ME, 
                                     NSIMS = MED_ITERS,
                                     var.prediction.results, 
                                     RESULTS_DIR = RESULTS_DIR, 
                                     controlIQ = F)

all.mediation.results = var.mediation.network.results
all.mediation.results[[length(var.network.results) + 1]] = var.mediation.results


all_module_labels = c(module_labels, "Full")
save.image(RESULTS_FILE)

