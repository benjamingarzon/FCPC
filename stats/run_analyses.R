#Of the healthy cohort (n=297), 2 subjects were excluded due to low quality images, 
#1 was excluded due to gross radiological abnormalities, 
#4 were excluded due to missing convergence in ME-ICA pre-processing, 
#and 9 were excluded due to excessive motion during the resting-state functional scan 
#(5 subjects with maximum framewise displacement larger than 1.3 mm and 
#4 subjects with mean framewise displacement #of 0.3 mm using calculation by Power et al 201236).

# 342 imaging
# 671 bhaviour

# 288 with imaging and behaviour

# corrected p-values for mediation

# check meica convergence
# confounds in BLCS model

# larger eta more sparse

# squareform for R, matlab and numpy are equivalent

#1: folds taking subjects into account, aggregate FD; brain, connnectivity

# compute corrected values separately for baseline and followup
# plot coefficients

# check that offset is used properly
# should work without covariates

# try with power parcellation?
######################

# trade-off : maxcomp=15, NITER=100, NPERM=200, NPROCS=35, NFOLDS=30, 
rm(list=ls())
options(warn = 0)
OUTPUT_DIR = "~/Data/NSPN/results/"

fig_count = 0
setwd("~/Software/FCPC/")
source("stats/aux_pls.R")
source("stats/prediction_funcs.R")
source("process/organize_data.R")

#########################################################
# define parameters
#########################################################
NPROCS = 33
N_FOLDS = 20
NITER = 50
NPERM1 = 500
NPERM2 = 500
if (F){
  # using Power parcelation
  NETMATS_FILE = 'netmats_Power_ledoit_Fischer.csv'
  SEPARATOR = " "
  parcel_file = file.path('~/Data/NSPN/templates/Power264.csv')
  parcel = read.table(parcel_file, sep = ";", header = T)
  parcel$Name <- gsub(" ", " ", parcel$Name)
  module_table = unique(parcel[c("Module", "Name", "Color")]) %>% arrange(Module)
  module_names = module_table$Module
  module_labels = module_table$Name
  
  modules_file = file.path('~/Data/NSPN/templates/modules.csv')
  write.table(parcel$Module, file = modules_file)
  write.table(module_table, file = '~/Data/NSPN/templates/module_table.csv')
  
  DIR_ICA_ME.BSL = '~/Data/NSPN/ICA_ME/nspn_ME/Power'
  DIR_ICA_ME.FU = '~/Data/NSPN/ICA_ME/nspn_ME_fu/Power'
  
  DIR_PREFFIX0 = paste0("prediction_data_nondep_Power_ledoit_", N_FOLDS, "folds_", NITER, "iters_")
  
} else {
  NETMATS_FILE = 'netmats_ridge.csv'
  SEPARATOR = ","
  
  # modules_file = file.path('~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica', 'modules_1.3.txt')
  # module_labels = c("SMT", #somatosensory + motor",
  #                   "IMOFC", # ventral frontal
  #                   "BGTMP", # basal ganglia + salience
  #                   "DLPFC", # dorsolateral prefrontal
  #                   "DMPFC", # dorsomedial prefrontal
  #                   "VIS",  # dorsomedial visual
  #                   "PDMN")  # default mode

  modules_file = file.path('~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica', 'modules_optimal.txt')
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
  
  DIR_ICA_ME.BSL = '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica'
  DIR_ICA_ME.FU = '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica'
#  DIR_ICA_ME.LONG = '~/Data/NSPN/ICA_ME/nspn_ME_long'
  
  DIR_PREFFIX0 = paste0("prediction_data_nondep_ICA_ridge_", N_FOLDS, "folds_", NITER, "iters_")
}

#########################################################

# ONLY NON-DEPRESSED SAMPLE
select_medical = "epression"

# list of modules
to_select_list = NULL
# enough to do it for the single networks
for (j in 1:max(module_names)) to_select_list[[j]] = c(j, j);
to_select_list[[max(module_names)+1]] = c(0, 0);
#for (j in 1:2) to_select_list[[j]] = c(j, j);

#############################################################
# Not correcting for age
#############################################################
NPERM = NPERM1
DIR_PREFFIX = paste0(DIR_PREFFIX0, "notagecorrected_")
data.pars = data.pars.nondep

print("------------------------------")
print("Not age corrected")
print("------------------------------")

##################
# decision acuity
##################
savecoefs = T

if (F){
# baseline
testmodels = NULL
DIR_ICA_ME = DIR_ICA_ME.BSL
var = "decAc.bsl"
age_var = 'age_scan1'
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "sex.num")
source("stats/prediction.R")

# project to followup
testmodels = var.network.results
rm(var.network.results)
DIR_ICA_ME = DIR_ICA_ME.FU
var = "decAc.fu"
age_var = 'age_scan2'
mycovars = c("brain_vol", "fd", "UCL.2", "CBU.2", "meica_dof", "sex.num")
source("stats/prediction.R")

rm(testmodels, var.network.results)

}
##############################
# latent change scores model
##############################
if (F){
RESULTS_FILE = file.path(OUTPUT_DIR, paste0(DIR_PREFFIX, "decAc.bsl"),  "data.RData") 
load(RESULTS_FILE)
results.df.bsl = results.df %>% dplyr::select(Subject, y.pred, y.test, label, age) 

RESULTS_FILE = file.path(OUTPUT_DIR, paste0(DIR_PREFFIX, "decAc.fu"),  "data.RData") 
load(RESULTS_FILE)
results.df.fu = results.df %>% dplyr::select(Subject, y.pred, y.test, label, age)

results.all = merge(results.df.bsl, results.df.fu, by = c("Subject", "label"), suffixes = c(".bsl", ".fu")) %>%
  rename(COG_T1 = y.test.bsl,
         COG_T2 = y.test.fu,
         NEU_T1 = y.pred.bsl,
         NEU_T2 = y.pred.fu,
         AGE_T1 = age.bsl,
         AGE_T2 = age.fu) %>% 
  mutate(DELTA_COG = COG_T2 - COG_T1,
         DELTA_NEU = NEU_T2 - NEU_T1)

results.sub = subset(results.all, label == "All")

plot(COG_T1 ~ COG_T2, data = results.sub, pch = 20)
plot(NEU_T1 ~ NEU_T2, data = results.sub, pch = 20)
plot(COG_T1 ~ COG_T2, data = results.sub, pch = 20)
plot(NEU_T1 ~ COG_T1, data = results.sub, pch = 20)
plot(DELTA_NEU ~ AGE_T1, data = results.sub, pch = 20)

points(NEU_T2 ~ COG_T2, data = results.sub,  pch = 20, col = "red")
plot(DELTA_NEU ~ DELTA_COG, data = results.sub)

par(mfrow = c(1, 2))
plot(DELTA_COG ~ AGE_T1, data = results.sub, pch = 20)
plot(DELTA_NEU ~ AGE_T1, data = results.sub, pch = 20)

colMeans(results.sub[-c(1, 2)])
par(mfrow = c(1, 2))
hist(results.sub$NEU_T1, 20, col = "blue")
hist(results.sub$NEU_T2, add=T, 20, col = rgb(1, 0, 0, 0.5))
hist(results.sub$COG_T1, 20, col = "blue")
hist(results.sub$COG_T2, add=T, 20, col = rgb(1, 0, 0, 0.5))


results.BLCS = NULL
for (mylabel in unique(results.all$label))
{
  print(mylabel)
  results.BLCS = rbind(results.BLCS, fitBLCSmodel(subset(results.all, label == mylabel)))
}

results.BLCS.change = subset(results.BLCS, type == "change")
results.BLCS.level = subset(results.BLCS, type == "level")
save(results.BLCS.change, results.BLCS.level, file = file.path(OUTPUT_DIR, paste0("BLCS_", DIR_PREFFIX, "decAc.RData")))
print(results.BLCS.change[c("label", "est", "se", "pvalue")])
# manually
results.corr = NULL
for (l in unique(results.all$label)) {
  X = subset(results.all, label == l)
  mycor = cor.test(X$DELTA_COG, X$DELTA_NEU)
#  print(paste(l, mycor$estimate, mycor$se, mycor$p.value))
  results.corr = rbind(results.corr, list(l, mycor$estimate, mycor$p.value))
  if (mycor$p.value < 0.05) plot(X$DELTA_COG, X$DELTA_NEU, type = "p", main = l, xlab = "d_2 - d_1 (measured)", ylab = "d_2 - d_1 (predicted)")
}
colnames(results.corr) = c("label", "est", "pvalue")
results.corr = as.data.frame(results.corr)

View(
  format(results.corr, digits = 1, nsmall= 2, scientific = F)
)
# print nicely
View(
  format(results.BLCS.change[c("label", "est", "se", "pvalue")], digits = 1, nsmall= 2, scientific = F)
)

View(
  format(results.BLCS.level[c("label", "est", "se", "pvalue")], digits = 1, nsmall= 2, scientific = F)
)

# do cross-sectional and longitudinal mediation

if (F){
results.bsl = results.df.bsl %>% 
  rename(COG = y.test, 
         NEU = y.pred,
         AGE = age)

results.long = merge(results.df.bsl, results.df.fu, by = c("Subject", "label"), suffixes = c(".bsl", ".fu")) %>% 
  mutate(COG = y.test.fu - y.test.bsl, 
         NEU = y.pred.fu - y.pred.bsl,
         AGE = age.fu - age.bsl)


cl <- makeCluster(NPROCS)
registerDoParallel(cl)

res.cross = foreach (mylabel = unique(results.bsl$label), .packages=c('mediation'), 
                     .export=c('do_mediation')) %dopar% do_mediation(subset(results.bsl, label == mylabel), mylabel)

res.long = foreach (mylabel = unique(results.long$label), .packages=c('mediation'), 
                     .export=c('do_mediation')) %do% do_mediation(subset(results.long, label == mylabel), mylabel)

stopCluster(cl)

results.mediation.cross = do.call("rbind", lapply(res.cross, function(x) x$mediation.age.coef))
results.mediation.long = do.call("rbind", lapply(res.long, function(x) x$mediation.age.coef))
results.model.cross = do.call("rbind", lapply(res.cross, function(x) x$model.results))
results.model.long = do.call("rbind", lapply(res.long, function(x) x$model.results))

print(subset(results.mediation.cross, param == "ACME" ))
print(subset(results.mediation.long, param == "ACME" ))
save(results.mediation.cross,
     results.mediation.long,
     results.model.cross,
     results.model.long,
     file = file.path(OUTPUT_DIR, paste0("mediation_", DIR_PREFFIX, "decAc.RData")))
}
}

NPERM = NPERM2
#############################################################
# Correcting for demographic factors, on complete data
#############################################################
data.pars = data.pars.nondep
DIR_PREFFIX = paste0(DIR_PREFFIX0, "democorrected_")

print("------------------------------")
print("Demo corrected")
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
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex")
source("stats/prediction.R")

# project to followup
testmodels = var.network.results
rm(var.network.results)
DIR_ICA_ME = DIR_ICA_ME.FU
var = "decAc.fu"
age_var = 'age_scan2'
mycovars = c("brain_vol", "fd", "UCL.2", "CBU.2", "meica_dof", "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex")
source("stats/prediction.R")

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
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex")
source("stats/prediction.R")

# project to followup
testmodels = var.network.results
rm(var.network.results)
DIR_ICA_ME = DIR_ICA_ME.FU
var = "IQcomp.fu"
age_var = 'age_scan2'
mycovars = c("brain_vol", "fd", "UCL.2", "CBU.2", "meica_dof", "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex")
source("stats/prediction.R")

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
savecoefs = F
# only baseline
testmodels = NULL
DIR_ICA_ME = DIR_ICA_ME.BSL
var = "decAc.bsl"
age_var = 'age_scan1'
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex", "IQcomp.bsl")
source("stats/prediction.R")
rm(var.network.results)

##################
# IQ composite
##################

# only baseline
DIR_ICA_ME = DIR_ICA_ME.BSL
var = "IQcomp.bsl"
age_var = 'age_scan1'
mycovars = c("brain_vol", "fd", "UCL.1", "CBU.1", "meica_dof", "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex", "decAc.bsl")
source("stats/prediction.R")
rm(var.network.results)

