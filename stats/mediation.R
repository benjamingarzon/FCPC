

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

# do mediation for each combination tested
var.mediation.network.results = NULL

for (i in seq(length(var.network.results))){
  print(i)
  var.mediation.network.results[[i]] = do_mediation(DIR_ICA_ME, 
                                                    var.network.results[[i]], 
                                                    NSIMS = MED_ITERS,
                                                    RESULTS_DIR = RESULTS_DIR)
}

save.image(RESULTS_FILE)

# do whole-brain prediction
var.prediction.results = do_prediction(DIR_ICA_ME, 
                                       demo.pars, 
                                       var, 
                                       mycovars, 
                                       N_FOLDS = N_FOLDS, 
                                       NPROCS = NPROCS, 
                                       NITER = NITER_FULL)

# calculate mediation
var.mediation.results = do_mediation(DIR_ICA_ME, 
                                     NSIMS = MED_ITERS,
                                     var.prediction.results, 
                                     RESULTS_DIR = RESULTS_DIR)

all.mediation.results = var.mediation.network.results
all.mediation.results[[i + 1]] = var.mediation.results

all_module_labels = c(module_labels, "Full")
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
