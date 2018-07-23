rm(list=ls())
options(warn = 0)
FIGS_DIR = "~/Data/NSPN/module_stats/"
OUTPUT_DIR = "~/Data/NSPN/module_stats/"
fig_count = 0
setwd("~/Software/FCPC/")
#source("~/Software/DAD/VariabilityAnalysis/aux_funcs.R")
#source("aux_analyze.R")
source("stats/aux_pls.R")
source("stats/prediction_funcs.R")
source("process/organize_data.R")


DIR_ICA_ME = '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica'

#modules_file = file.path(DIR_ICA_ME, 'BrainMapPartition.txt')
#module_names = read.table('~/Data/NSPN/BrainMapNames.txt', sep = ';', header = T)

modules_file = file.path(DIR_ICA_ME, 'modules.txt')
modules = read.table(modules_file)
module_names = seq(max(modules))
#read.table('~/Data/NSPN/module_names.txt', sep = ';', header = T)


NPROCS = 10 #30
NITER = 5 # 200 #
N_FOLDS = 3# 50
mycovars = c("brain_vol", "fd", "UCL", "CBU", "meica_dof")#, "TIV", "processing")
#to_select_list = as.list(c(seq(4), seq(6, 18))) # check one module per iteration
#to_select_list = as.list(c(seq(11))) # check one module per iteration
to_select_list = NULL
k = 1
for (i in 1:length(module_names))
  for (j in 1:length(module_names))
    if (i >= j)   {
      to_select_list[[k]] = c(i, j);
      k = k + 1
    } 

RESULTS_FILE = './mediation_data.RData'

var = "decAc"


# do prediction for different networks
decAc.network.results = do_predictions_loop(DIR_ICA_ME, 
                                         demo.pars.NSPN,
                                         var,
                                         mycovars,
                                         modules_file,
                                         to_select_list,
                                         NPROCS = NPROCS, 
                                         N_FOLDS = N_FOLDS, 
                                         NITER = NITER
                                         )

#for (i in seq(length(decAc.network.results))) decAc.network.results[[i]]$to_select = to_select_list[[i]] 
rho = data.frame(
                 mean = sapply(decAc.network.results, function(x) ifelse(!x$empty, x$cor.test$estimate, NA)), 
                 low = sapply(decAc.network.results, function(x) ifelse(!x$empty, x$cor.test$conf.int[1], NA)),
                 high = sapply(decAc.network.results, function(x) ifelse(!x$empty, x$cor.test$conf.int[2], NA)),
                 p = sapply(decAc.network.results, function(x) ifelse(!x$empty, x$cor.test$p.value, NA)),
                 nfeatures = sapply(decAc.network.results, function(x) ifelse(!x$empty, x$nfeatures, NA)),
                 from = as.factor(sapply(decAc.network.results, function(x) x$to_select[1])),
                 to = as.factor(sapply(decAc.network.results, function(x) x$to_select[2]))                                      
                  )        

# check the prediction is not driven by the number of features
print(cor.test(rho$nfeatures, abs(rho$mean)))  
print(cor.test(rho$nfeatures, rho$mean))  

ggplot(rho, aes( x = nfeatures, y = mean, ymin = low, ymax = high )) + 
  geom_pointrange() + 
  xlab("Number of features") +
  ylab("Correlation between true and predicted") + 
  theme_classic() 


rho$mean[rho$p > 0.05] = NA

ggplot(rho, aes(x = from, y = to, fill = mean)) +
  #  ggtitle('Correlation matrix  between true and predicted') +
  geom_tile() +
  labs(title = 'Correlation matrix  between true and predicted',
       x = 'Network',
       y = 'Network',
       fill = 'Correlation'
  ) + 
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  scale_fill_gradient2(midpoint = .1, mid ="orange", high = "yellow", low = "red", limits = c(0, +.2)) +
  geom_text(aes(x = from, y = to, label = round(mean, 2)), color = "black", fontface = "bold", size = 3) +
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.5, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5) )


save.image(RESULTS_FILE)

# do whole-brain prediction
decAc.prediction.results = do_prediction(DIR_ICA_ME, 
                                         demo.pars.NSPN, 
                                         var, 
                                         mycovars, 
                                         N_FOLDS = N_FOLDS, 
                                         NPROCS = NPROCS, 
                                         NITER = NITER)

# calculate mediation
decAc.mediation.results = do_mediation(DIR_ICA_ME, decAc.prediction.results)

save.image(RESULTS_FILE)


save.image(RESULTS_FILE)

# do check the effect of increasing the number of 

NITERS = c(50, 100, 150, 200, 250)
iter.results = list()
for (NITER in seq(NITERS)) iter.results = c(iter_results,
                                       do_prediction(DIR_ICA_ME, 
                                         demo.pars.NSPN, 
                                         var, 
                                         mycovars, 
                                         N_FOLDS = N_FOLDS, 
                                         NPROCS = NPROCS, 
                                         NITER = NITER)
)

iter.rho = sapply(iter.results, function(x) ifelse(!x$empty, x$cor.test$estimate, NA))
plot(NITERS, iter.rho, xlab = "NITER", ylab = "Correlation between true and predicted")


stophere
# rho = t(sapply(decAc.network.results, function(x) x$rho )

# rho.mean.order = order(rho$mean)
# mod_names = module_names$Name[unlist(to_select_list)]
# rho$name = factor(mod_names, levels = mod_names[rho.mean.order]) 
# 
# ggplot(rho, aes( x = name, y = mean, ymin = low, ymax = high )) + 
#   geom_pointrange() + 
#   geom_hline(yintercept = 0, lwd = 0.2, lty = 2) +
#   xlab("Network") +
#   ylab("Correlation between true and predicted") + 
#   theme_classic() + 
#   coord_flip()
# 

# rho.mean = rowMeans(rho)
# rho.mean.order = order(rho.mean)
# rownames(rho) = module_names$Name[unlist(to_select_list)] 
# rho = rho[rho.mean.order, ]
# rho.melt = melt(rho)[c(1, 3)]
# colnames(rho.melt) = c('net', 'value')
# 
# ggplot(rho.melt, aes(y = value, x = as.factor(net))) + geom_boxplot() + 
#   #  geom_jitter(height = 0, width = 0.1, size = 0.5) + 
#   xlab("Network") +
#   ylab("Correlation between true and predicted") + 
#   coord_flip() +
#   theme_classic()



stophere

save.image(RESULTS_FILE)
var = "decAc"
decAc.hierarchical.results = do_spls_hierarchical(DIR_ICA_ME, 
                                                  demo.pars.NSPN, 
                                                  var, 
                                                  mycovars, 
                                                  NPROCS = NPROCS, 
                                                  N_FOLDS = N_FOLDS, 
                                                  TEST_PROP = 0.05, 
                                                  NITER = NITER)

save.image(RESULTS_FILE)

load(RESULTS_FILE)
plot_results_hierarchical(decAc.hierarchical.results)


# # SEX AGE INTERACTION
# var = "decAc"
# mylist = list()
# for( i in seq(5)){
# mylist[[i]] = do_spls_hierarchical(DIR_ICA_ME, demo.pars.NSPN, 
#                                        var, mycovars, NPROCS = NPROCS, N_FOLDS = N_FOLDS, 
#                                        TEST_PROP = 0.05, NITER = NITER)
# }
# 
# unlist(lapply(seq(5), function (x) mylist[[x]]$covars.cor$estimate))
# unlist(lapply(seq(5), function (x) mylist[[x]]$full.cor$estimate))
# unlist(lapply(seq(5), function (x) mylist[[x]]$nets.cor$estimate))
# 
#decAc.hierachical.results10.1

# var = "IQ_matrix"
# IQ_matrix.mediation.results = do_mediation(DIR_ICA_ME, demo.pars.NSPN, 
#                                            var, mycovars, NPROCS = NPROCS, N_FOLDS = N_FOLDS, 
#                                            TEST_PROP = 0.05, NITER = NITER)
# 
# var = "IQ_vocab"
# IQ_vocab.mediation.results = do_mediation(DIR_ICA_ME, demo.pars.NSPN, 
#                                           var, mycovars, NPROCS = NPROCS, N_FOLDS = N_FOLDS, 
#                                           TEST_PROP = 0.05, NITER = NITER)



# save.image(RESULTS_FILE)



# var = "IQ_matrix"
# IQ_matrix.hierarchical.results = do_spls_hierarchical(DIR_ICA_ME, demo.pars.NSPN, 
#                                            var, mycovars, NPROCS = NPROCS, N_FOLDS = N_FOLDS, 
#                                             TEST_PROP = 0.05, NITER = NITER)
# var = "IQ_vocab"
# IQ_vocab.hierarchical.results = do_spls_hierarchical(DIR_ICA_ME, demo.pars.NSPN, 
#                                           var, mycovars, NPROCS = NPROCS, N_FOLDS = N_FOLDS, 
#                                           TEST_PROP = 0.05, NITER = NITER)
# 
# 
# save.image(RESULTS_FILE)




# plot_results_hierarchical(IQ_matrix.hierarchical.results)
# plot_results_hierarchical(IQ_vocab.hierarchical.results)

