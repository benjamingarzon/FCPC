######## Prediction functions ##########################

library(pryr)
library(spls)
#library(pls)
library(mediation)
library(foreach)
library(doParallel)
library(doMC)
library(caret)
library(pracma)
library(reshape2)
library(dplyr)
library(oro.nifti)
library(RColorBrewer)
library(scales)
###library(ppcor)
library(lavaan)
library(RVAideMemoire)
library(psychometric)
library(boot)
library(rlist)
library(psych)
library(rjson)
library(abind)
FOLDSEED = sample(seq(10000), 1)

MYPALETTE=brewer.pal(11, "Spectral")[seq(11, 1, -1)]# colorspace::diverge_hsv(15, power=2)

FINELINE = 0.5
GGTEXTSIZE2 = 10
GGTEXTSIZE3 = 10
NDIGITS = 3

DEFAULTMAXCOMP = 10
NPERMSAPPROX = 10000
###############################################
# save figure
###############################################
save_fig = function(figname=NULL, width=6.5, height=6.5, res=600, jpg=F){
  
  #size in inches
  if (dev.cur()!=1) dev.off()
  if (is.null(figname)) {
    figname = paste0('plot', fig_count)
    
    fig_count <<- fig_count+1
  }
  
  print(paste("Generating figure: ", figname))
  if (jpg){
    figname = file.path(FIGS_DIR, paste0( figname, '.jpg'))
    jpeg(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  } else {
    figname = file.path(FIGS_DIR, paste0( figname, '.png'))
    png(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  }
  
}

test_cor = function(N, s){
  # r = rep(0, NPERMSAPPROX)
  # for (i in seq(NPERMSAPPROX)){
  #   x = rnorm(N)
  #   y = rnorm(N)
  #   r[i] = cor.test(x, y)$estimate
  # }
  # hist(r, 100)
  # hist(fisherz2r(rnorm(NPERMSAPPROX)/sqrt(N - 3)), 100, add = T, col = rgb(1, 0, 0, 0.5))
  
  X = mvrnorm(n = N, mu = c(0, 0), Sigma = matrix(c(1, s, s, 1), 2, 2))
  mycor = cor.test(X[,1], X[,2])
  perm = fisherz2r(rnorm(NPERMSAPPROX)/sqrt(N - 3))
  perm[1] = mycor$estimate
  pval = sum(perm[1] <= perm)/length(perm)
  #  print(mycor)
  #  print(pval)
  return(c(mycor$p.value, pval))
}
#xx = matrix(0, 200, 2)
#for (i in seq(200)) xx[i, ] = test_cor(250, .24)
#scatterHist(xx)

r2score = function(x, y){
  #  r2 = sum(z^2)/sum((y - mean(y))^2)
  
  y = y - mean(y) # we don't care about the offset
  z = lm(y ~ x)$residuals
  r2 = 1 - var(z)/var(y) 
  return(r2)
}

squareform2 = function (x){
  if (is.vector(x)) {
    n <- length(x)
    m <- floor(sqrt(2 * n))
    if (m * (m + 1) != 2 * n) 
      stop("Argument 'x' does not correspond to a distance matrix.")
    inds <- c()
    k <- m + 1
    for (i in 1:(k - 1)) inds <- c(inds, (1 + i + (i - 1) * 
                                            k):(i * k))
    y <- numeric(k * k)
    y[inds] <- x
    y <- matrix(y, k, k) + t(matrix(y, k, k))
  }
  else if (is.matrix(x)) {
    m <- nrow(x)
    n <- ncol(x)
    if (m != n) 
      stop("Argument 'x' must be a vector or a square matrix.")
    if (n == 1) 
      return(c())
    inds <- c()
    for (i in 1:(n - 1)) inds <- c(inds, (1 + i + (i - 1) * 
                                            n):(i * n))
    y <- x[inds]
  }
  return(y)
}


###############################################
# select connections in a module
###############################################

resultstodf = function(results, number){
  results$demo$sex.num = NULL
  results.df = cbind(results[ c("Subject", "y.pred", "y.test")], results$covars.exp, results$demo)
  results.df$var = results$var
  results.df$number = number
  results.df$label = results$label
  
  if (!is.null(results$group)) results.df$group = results$group
  
  return(as.data.frame(results.df))
}


###############################################
# select connections in a module
###############################################
select_modules = function(modules_file, to_select, invert = F){
  
  modules = unlist(read.table(modules_file))
  #print(modules)
  n.rois = length(modules)
  if(invert) {
    adj = matrix(T, n.rois, n.rois)
    adj[modules == to_select[1], ] = F
    adj[, modules == to_select[1]] = F
    adj[modules == to_select[2], ] = F
    adj[, modules == to_select[2]] = F
    
  } else {
    adj = matrix(F, n.rois, n.rois)
    
    adj[modules == to_select[1], ] = T
    adj[, modules == to_select[1]] = T
    adj[modules == to_select[2], ] = T
    adj[, modules == to_select[2]] = T
    
    if (to_select[1] == 0 ){
      adj[modules == to_select[2], ] = T
      adj[, modules == to_select[2]] = T
    }
    if (to_select[2] == 0){
      adj[modules == to_select[1], ] = T
      adj[, modules == to_select[1]] = T
    }
    if (to_select[1] == 0 & to_select[2] == 0){
      adj[, ] = T
    }
    
  }  
  
  valid = which(squareform2(adj))
  print(paste("Connections between module(s) ", to_select[1], ", " , to_select[2], ": " , length(valid), "Inverted: ", invert))
  return(valid)
}

###############################################
# run prediction for different partitions
###############################################
do_predictions_loop = function(DIR_ICA_ME, # working dir
                               data.pars, # predictors
                               var, # independent variable
                               mycovars, # covariates
                               modules_file,
                               to_select_list,
                               maxcomp = DEFAULTMAXCOMP, # max components for spls
                               N_FOLDS = 10, # number of xvalidation folds
                               NPROCS = 10, # processors to use
                               NITER = 50, # spls bagging iteration
                               invert = F, 
                               netmats_file = '',
                               age_var = '',
                               NPERM = 100,
                               mymodels = NULL, 
                               savecoefs = F)
{
  
  results = list()
  for (i in seq(length(to_select_list)))
  {
    to_select = to_select_list[[i]]
    #print(to_select)
    if(!is.null(mymodels)) mymodel = mymodels[[i]] else mymodel = NULL
    
    results[[i]] = do_prediction(DIR_ICA_ME, 
                                 data.pars, 
                                 var,
                                 mycovars,
                                 maxcomp,
                                 N_FOLDS,
                                 NPROCS,
                                 NITER,
                                 modules_file,
                                 to_select, 
                                 invert, 
                                 netmats_file,
                                 age_var,
                                 NPERM = NPERM,
                                 mymodel = mymodel,
                                 savecoefs = savecoefs)
    #    View(results)
    print(paste("Done", i/length(to_select_list)*100, "%"))
    
  }
  
  #  print(mean(results[[i]]$rho))
  return(results)
}

apply_model = function(X, y, covars, model, mysubjects){
  # for each fold, find model and average across folds
  y.test.mat = y.pred.mat = mysubjects.test = NULL
  for (i in seq(length(model$preprocessing))){
    trainers = match(model$training_subjects[[i]], mysubjects)
    trainers = trainers[!is.na(trainers)]
    X.test = X[-trainers, ]
    mysubjects.test[[i]] = mysubjects[-trainers]   

    if (!is.null(covars)) y.test = y - predict(model$preprocessing[[i]]$mycovarmod, covars)
    else y.test = y

    mu = model$preprocessing[[i]]$mu
    sigma = model$preprocessing[[i]]$sigma
    
    for (j in seq(ncol(X))){
      X.test[, j] = (X.test[, j] - mu[j])/sigma[j]     
    }

    y.pred.perm = matrix(NA, nrow(X), ncol(model$coefs.perm[, , i]))
    y.pred.perm[-trainers, ] = X.test %*% model$coefs.perm[, , i] + model$offset[i] 
    y.test.mat = cbind(y.test.mat, y.test)
    y.pred.mat = abind(y.pred.mat, y.pred.perm, along = 3)
  }  
  return(list(y.test = rowMeans(y.test.mat), 
              y.pred.perm = apply(y.pred.mat, c(1, 2), mean, na.rm = T),
              y.test.orig = y,
              testing_subjects = mysubjects.test)
         )
}

###############################################
# fit and test predictive model
###############################################

do_prediction = function(DIR_ICA_ME, # working dir
                         data.pars, # predictors
                         var, # independent variable
                         mycovars, # covariates
                         maxcomp = DEFAULTMAXCOMP, # max components for spls
                         N_FOLDS = 10, # number of xvalidation folds
                         NPROCS = 10, # processors to use
                         NITER = 50, # spls bagging iteration
                         modules_file = NULL,
                         to_select = NULL, 
                         invert =  F, 
                         netmats_file = 'netmats_ridge.csv',
                         age_var = '',
                         fdfiltering = 0.3,
                         longformat = F, 
                         NPERM = 100,
                         mymodel = NULL, savecoefs = F
                         # instead of estimating, apply a previous model
)
{ 
  print("-----------------------------------------------------------------")
  print(var)
  
  tic("model fitting")
  results = list()
  results$var = var
  results$empty = F
  # needed files
  FILE_ICA_ME = file.path(DIR_ICA_ME, netmats_file) #  FILE_ICA_ME = file.path(DIR_ICA_ME, 'metrics_ridge.csv')
  FILE_BRAIN_VOL = file.path(DIR_ICA_ME, 'brain_volumes.txt')
  FILE_DOF = file.path(DIR_ICA_ME, 'meica_DOF_nom.txt')
  FILE_FD = file.path(DIR_ICA_ME, 'fd.txt')
  FILE_DVARS = file.path(DIR_ICA_ME, 'dvars.txt')
  FILE_SUBJECTS = file.path(DIR_ICA_ME, 'subjects.txt')
  USE_COVARS = ifelse(!is.null(mycovars), 1, 0)
  # load and organize covariates
  brain_vol = read.table(FILE_BRAIN_VOL)
  meica_dof = read.table(FILE_DOF)
  fd = read.table(FILE_FD)
  dvars = read.table(FILE_DVARS)
  ica_data = read.table(FILE_ICA_ME, sep = SEPARATOR, head=F)
  subjects = as.numeric(read.table(FILE_SUBJECTS))
  
#  browser()
  # sync data and covariates
  if (longformat){  
    fd = cbind(fd[c(1, 2)], rowMeans(fd[-c(1, 2, 3)]))
    dvars = cbind(fd[c(1, 2)], rowMeans(dvars[-c(1, 2, 3)]))
    colnames(ica_data) = c("group", paste0("feature.", seq(ncol(ica_data)-1)))
    ica_data$Subject = subjects
    
    mydata = merge(ica_data, data.pars, by=c("Subject", "group"))
    
    covars = merge(dvars, fd, by = c("V1", "V2"))
    covars = merge(covars, brain_vol, by = c("V1", "V2"))
    covars = merge(covars, meica_dof, by = c("V1", "V2"))
    colnames(covars) = c("group", "Subject","dvars", "fd","brain_vol","meica_dof")
    mydata = merge(mydata, covars, by = c("Subject", "group"))
    
  } else {
    fd = cbind(fd[1], rowMeans(fd[-c(1,2)]))
    dvars = cbind(fd[1], rowMeans(dvars[-c(1,2)]))
    colnames(ica_data) = paste0("feature.", seq(ncol(ica_data)))
    ica_data$Subject = subjects
    
    mydata = merge(ica_data, data.pars, by="Subject")
    
    covars = merge(dvars, fd, by = "V1")
    covars = merge(covars, brain_vol, by = "V1")
    covars = merge(covars, meica_dof, by = "V1")
    colnames(covars) = c("Subject","dvars", "fd","brain_vol","meica_dof")
    mydata = merge(mydata, covars, by ="Subject")
  }
  features = grep('feature', colnames(mydata))
  
  # select features for valid modules
  if (!is.null(modules_file)) 
  {
    valid.cols = select_modules(modules_file, to_select, invert)
    features = features[valid.cols]
    results$to_select = to_select
    results$inverted = invert
    
    if (length(valid.cols) < 5){
      print("Not enough features!")
      results$empty = T
      return(results)
    }
  } 
  mydata$age = unlist(mydata[age_var])
  
  # add higher order terms and interactions
  mydata$age.1 = as.numeric(scale(mydata$age))
  mydata$age.2 = mydata$age.1^2
  mydata$age.1.sex = mydata$age.1*mydata$sex.num
  mydata$age.2.sex = mydata$age.2*mydata$sex.num
  
  mydata.exp = cbind(
    mydata[features],
    mydata[c(var, unique(c(mycovars, "age", "sex.num")))]   
  ) 
  mydata.exp$Subject = mydata$Subject

  print(dim(mydata))
  
  if (longformat) mydata.exp$group = mydata$group
  
  if (!is.null(fdfiltering)){
    print(paste("Removing", sum(mydata.exp$fd > fdfiltering), "subjects with high mean FD above", fdfiltering ))
    results$fdfiltered = subset(mydata.exp, fd > fdfiltering)$Subject
    mydata.exp = subset(mydata.exp, fd <= fdfiltering)
  }
  
  mydata.exp = mydata.exp[complete.cases(mydata.exp), ]
  y = unlist(mydata.exp[var])
  
  #remove outliers if any
  valid.rows = (y - mean(y)) < 3*sd(y)
  print(paste("Removed", sum(!valid.rows), "outliers of", var))
  
  covars.exp = mydata.exp[valid.rows, mycovars]
  if ("UCL.1" %in% colnames(covars.exp)) covars.exp = covars.exp %>% rename(UCL = UCL.1, CBU = CBU.1) # ugly hack...
  if ("UCL.2" %in% colnames(covars.exp)) covars.exp = covars.exp %>% rename(UCL = UCL.2, CBU = CBU.2) # ugly hack...
  
  features = grep('feature', colnames(mydata.exp))
  ica_data = as.matrix(mydata.exp[valid.rows, features])
  y = unlist(mydata.exp[var])[valid.rows]
  print(dim(ica_data))
  
  set.seed(FOLDSEED) # so that the folds are in sync for different structures
  
  if (N_FOLDS == 0) N_FOLDS = length(y)
  
  if (longformat) {
    unique_subjects = unique(mydata.exp$Subject)
    folds_sub = createFolds(unique_subjects, k = N_FOLDS)
    folds = NULL
    for (i in seq(length(folds_sub))) {
      folds = c(folds, list(Fold = which(mydata.exp$Subject %in% unique_subjects[folds_sub[[i]]])))
    }
  }
  else folds = createFolds(y, k = N_FOLDS)
  
  mysubjects = mydata.exp$Subject[valid.rows]
  mygroup = mydata.exp$group[valid.rows]
  if (is.null(mymodel)){
    
    cl <- makeCluster(NPROCS)
    registerDoParallel(cl)
    cv = foreach(fold_index = seq(length(folds))) %do% 
      do_crossvalidate_spls_covars_perm_par(fold_index, list(X = ica_data, y = y, covars = covars.exp, 
                                                             folds = folds, 
                                                             subject = mysubjects,
                                                             group = mygroup), 
                                            maxcomp = maxcomp, cluster = cl, 
                                            NITER = NITER, NPERM = NPERM)
    
    stopCluster(cl)
    # extract and reorder so that it's in sync
    results$fold = unlist(lapply(cv, function(x){x$fold}))
    results$training_subjects = lapply(cv, function(x){mysubjects[-x$fold]})
    results$y.pred.perm = list.rbind(lapply(cv, function(x){x$y.pred.perm}))[order(results$fold), ]
    results$y.test.orig = unlist(lapply(cv, function(x){x$y.test.orig}))[order(results$fold)]
    results$y.pred = unlist(lapply(cv, function(x){x$y.pred}))[order(results$fold)]
    results$y.test = unlist(lapply(cv, function(x){x$y.test}))[order(results$fold)]
    results$nfeatures.perm = list.rbind(lapply(cv, function(x){x$nfeatures.perm}))  
    results$mysample.perm = lapply(cv, function(x){as.integer(x$mysample.perm)})  
    results$eta.iter = unlist(lapply(cv, function(x){x$eta.iter}))      
    results$K.iter = unlist(lapply(cv, function(x){x$K.iter}))  
    results$offset = unlist(lapply(cv, function(x){x$offset}))  
    results$coefs.cv = do.call("cbind", lapply(cv, function(x){x$coefs}) )
    
    results$coefs.mean = rowMeans(results$coefs.cv)
    results$coefs.stable = apply(results$coefs.cv, 1,
                                 function(x) ifelse( prod(range(x)) < 0, 0, mean(x)))
    
    # save these to disk to spare memory
    if (savecoefs) {
    coefs.perm = do.call("abind", list(lapply(cv, function(x){x$coefs.perm}), along = 3))
    results$coefs.perm.file = tempfile(pattern = "coefs", fileext = ".rda") 
    save(coefs.perm, file = results$coefs.perm.file) 
    
    results$preprocessing.file = tempfile(pattern = "preproc", fileext = ".rda")
    preprocessing = lapply(cv, function(x) x$preprocessing) 
    save(preprocessing, file = results$preprocessing.file) 
    }
    results$projected = F  
  } else{
    # do not estimate, apply previously computed model
    results$projected = T 
    results$fold = NULL
    results$nfeatures.perm = mymodel$nfeatures.perm  
    results$eta.iter = mymodel$eta.iter  
    results$K.iter = mymodel$K.iter
    
    load(mymodel$coefs.perm.file)
    load(mymodel$preprocessing.file)
    mymodel$coefs.perm = coefs.perm
    mymodel$preprocessing = preprocessing
    
    results$coefs.mean = mymodel$coefs.mean
    results$coefs.stable = mymodel$coefs.stable
#    browser()
    
    results.fit = apply_model(ica_data, y, covars.exp, mymodel, mysubjects)
    results$y.test = results.fit$y.test
    results$mysubjects.test = results.fit$mysubjects.test
    results$y.test.orig = results.fit$y.test.orig
    results$y.pred.perm = results.fit$y.pred.perm
    results$y.pred = results.fit$y.pred.perm[, 1]
  }
  
  results$covars.exp = covars.exp
  results$demo = mydata.exp[valid.rows, c("age", "sex.num")]
  results$nfeatures = ncol(ica_data)
  results$Subject = mydata.exp$Subject[valid.rows]
  
  if (longformat) results$group = mydata.exp$group
  
  results$cor.test = results.R2.test = results$cor.test.approx = NULL
  results$cor.test$perm = apply(results$y.pred.perm, 2, function(x) cor.test(x, results$y.test)$estimate)
  results$cor.test$estimate = results$cor.test$perm[1]
  results$cor.test$p.value = sum(results$cor.test$perm[1] <= results$cor.test$perm)/length(results$cor.test$perm)
  
  # so that one can also show null distributions
  results$cor.test.approx$perm = fisherz2r(rnorm(NPERMSAPPROX)/sqrt(length(results$y.test) - 3))
  results$cor.test.approx$perm[1] = results$cor.test$perm[1]
  results$cor.test.approx$estimate = results$cor.test$perm[1]
  results$cor.test.approx$p.value.perm = sum(results$cor.test$perm[1] <= results$cor.test.approx$perm)/length(results$cor.test.approx$perm)
  results$cor.test.approx$p.value = cor.test(results$y.pred, results$y.test)$p.value
  
  results$R2.test$perm = apply(results$y.pred.perm, 2, function(x) r2score(x, results$y.test))
  results$R2.test$estimate = results$R2.test$perm[1]
  results$R2.test$p.value = sum(results$R2.test$perm[1] <= results$R2.test$perm)/length(results$R2.test$perm)
  
  
  print(paste0("Estimate (exact):", results$cor.test$estimate, ",", results$cor.test$p.value ))
  print(paste0("Estimate (approx):", results$cor.test.approx$estimate, ",", results$cor.test.approx$p.value, ",", results$cor.test.approx$p.value.perm))
  print(paste0("R2:", results$R2.test$estimate, ",", results$R2.test$p.value ))
  
  results$time = toc()

  print(paste0("Memory used:", mem_used()))
  print(paste0("Tempdir:", tempdir()))
  return(results)
}


###############################################
# mediation analysis
###############################################
do_mediation = function(results.df,
                        label,
                        NSIMS = 2000 #,
                        #                            write_coefs = F,
                        #                            RESULTS_DIR = file.path(DIR_ICA_ME, 'mediation/')
)
{ 
  
  mydata = results.df
  #  browser()
  # cross-sectional mediation
  print(paste("Analyzing", nrow(mydata), "observations."))
  # 1) X -> Y
  model.0 = lm(COG ~ AGE, data = mydata)  
  
  # 2) X -> M
  model.M = lm(NEU ~ AGE, data = mydata)
  
  # 3) X + M -> Y
  model.Y = lm(COG ~ AGE + NEU, data = mydata)  
  
  print(summary(model.0))
  print(summary(model.M))
  print(summary(model.Y))
  
  # longitudinal mediation
  
  # mediation 
  mediation.age = mediation::mediate(model.M, model.Y, treat = "AGE", mediator = "NEU", boot = T, sims = NSIMS)
  
  print(summary(mediation.age))
  #  results$mediation.age.p = mediation.age$d0.p
  
  model.0.results = round(summary(model.0)$coefficients[c('AGE'), c(1, 3, 4)], digits = NDIGITS)
  model.M.results = round(summary(model.M)$coefficients[c('AGE'), c(1, 3, 4)], digits = NDIGITS)
  model.Y.results = round(summary(model.Y)$coefficients[c('AGE', 'NEU'), c(1, 3, 4)], digits = NDIGITS)
  names(model.0.results) = names(model.M.results) = names(model.Y.results) = c("Estimate", "t", "p")
  model.results = as.data.frame(rbind(model.0.results, model.M.results, model.Y.results))
  model.results$param = c("AGE.0", "NEU.M", "AGE.Y", "NEU.Y")
  
  mediation.age.sum = summary(mediation.age) # ACME d.avg, d0.ci, d0.p /ADE z... / TOTAL tau /prop n.avg
  
  mediation.age.coef = as.data.frame(rbind(unlist(mediation.age.sum[ c('d.avg', 'd.avg.ci', 'd.avg.p')]),
                                           unlist(mediation.age.sum[ c('z.avg', 'z.avg.ci', 'z.avg.p')]),
                                           unlist(mediation.age.sum[ c('tau.coef', 'tau.ci', 'tau.p')]),
                                           unlist(mediation.age.sum[ c('n.avg', 'n.avg.ci', 'n.avg.p')])))
  
  colnames(mediation.age.coef) = c("Average", "Cil", "Cih", "p")
  mediation.age.coef$param = c("ACME", "ADE", "Total", "Proportion")
  
  mediation.age.coef$label = label
  model.results$label = label
  
  return(list(mediation.age.coef = mediation.age.coef, model.results = model.results))
}



###############################################
# compute results but corrected for covariates
###############################################
do_correction = function(DIR_ICA_ME, 
                         results, 
                         data.pars,
                         controlfor = NULL)
{ 
  #browser()
  data.M = results$covars.exp
  data.M$y.obs = results$y.test.mean
  data.M$y.pred = results$y.pred.mean
  data.M$Subject = results$Subject
  data.all = merge(data.M, data.pars, by = "Subject", suffixes = c("", ".y"), no.dups = T)
  #  data.all$age.1 = scale(data.all$age)
  #  data.all$age.1.sex = data.all$age.1*data.all$sex.num
  
  #  results$cor.test.twotailed = pcor.test(data.M$y.obs, data.M$y.pred, data.all[c("age.1", "sex.num", "age.1.sex", controlfor)]) 
  results$cor.test.twotailed = pcor.test(data.M$y.obs, data.M$y.pred, data.M[c("age.1", "sex.num", "age.1.sex", controlfor)]) 
  #  results$cor.test.twotailed$conf.int = c(0, 0)
  
  return(results)
}


addstar = function(x, fdr_correct) {
  
  fdrp = p.adjust(x[fdr_correct], method = "fdr")
  sig = ifelse(x < alpha, "* ", "  ")
  sig[fdr_correct][fdrp < alpha] = "**"
  return(paste0(x, sig))
}

###############################################
# make a table with results of mediation analyses
###############################################
make_mediation_table = function(results, module_labels, fdr_correct = seq(len(module_names))){
  
  # change col labels
  model.0.results.age = as.data.frame(t(sapply(results, function(x) x$model.0.results['age.1', ]))) %>% mutate(p = addstar(p, fdr_correct))
  model.M.results.age = as.data.frame(t(sapply(results, function(x) x$model.M.results['age.1', ]))) %>% mutate(p = addstar(p, fdr_correct))
  model.Y.results.age = as.data.frame(t(sapply(results, function(x) x$model.Y.results['age.1', ]))) %>% mutate(p = addstar(p, fdr_correct))
  model.0.results.sex = as.data.frame(t(sapply(results, function(x) x$model.0.results['sex1', ]))) %>% mutate(p = addstar(p, fdr_correct))
  model.M.results.sex = as.data.frame(t(sapply(results, function(x) x$model.M.results['sex1', ]))) %>% mutate(p = addstar(p, fdr_correct))
  model.Y.results.sex = as.data.frame(t(sapply(results, function(x) x$model.Y.results['sex1', ]))) %>% mutate(p = addstar(p, fdr_correct))
  
  model.Y.results.M = as.data.frame(t(sapply(results, function(x) x$model.Y.results['M', ]))) %>% mutate(p = addstar(p, fdr_correct))
  
  mediation.age.ACME = round(as.data.frame(t(sapply(results, 
                                                    function(x) x$mediation.age['ACME', ]))), 3) %>% 
    mutate(p = addstar(p, fdr_correct), Cil = paste0('(', Cil, ', ', Cih, ')')) %>% rename(CI = Cil) %>% dplyr::select(-one_of("Cih"))
  
  mediation.age.prop = round(as.data.frame(t(sapply(results, 
                                                    function(x) x$mediation.age['Proportion', ]))), 3) %>% 
    mutate(p = addstar(p, fdr_correct), Average = ifelse(Average >0, Average, 0))
  
  mediation.sex.ACME = round(as.data.frame(t(sapply(results, function(x) x$mediation.sex['ACME', ]))) , 3) %>%
    mutate(p = addstar(p, fdr_correct), Cil = paste0('(', Cil, ', ', Cih, ')')) %>% rename(CI = Cil) %>% dplyr::select(-one_of("Cih"))
  
  mediation.sex.prop = round(as.data.frame(t(sapply(results, function(x) x$mediation.sex['Proportion', ]))), 3) %>% 
    mutate(p = addstar(p, fdr_correct), Average = ifelse(Average >0, Average, 0))
  
  mediation.age.ACME$Proportion = mediation.age.prop$Average
  mediation.sex.ACME$Proportion = mediation.sex.prop$Average
  
  ## model.M, model.0, model.Y
  myorder = order(module_labels)
  dt.Y.M = cbind(module_labels, model.Y.results.M[, -1]) [myorder, ]
  dt.age = cbind(module_labels, model.M.results.age[, -1], model.Y.results.age[, -1], mediation.age.ACME[, -1])[myorder, ]
  dt.sex = cbind(module_labels, model.M.results.sex[, -1], model.Y.results.sex[, -1], mediation.sex.ACME[, -1])[myorder, ]
  
  print("--------------")
  print("M coef model.Y")
  print(kable(dt.Y.M),  format = 'html')
  #  print(dt.Y.M)
  print("--------------")
  print("Age coefs model.0")
  print(model.0.results.age[1, ])
  print("Age coefs model.M, model.Y")
  print(kable(dt.age), format = 'html')
  #  print(dt.age)
  print("--------------")
  print("Sex coefs model.0")
  print(model.0.results.sex[1, ])
  print("Sex coefs model.M, model.Y")
  print(kable(dt.sex),  format = 'html')
  #  print(dt.sex)
  
}


###############################################
# plot results for modules
###############################################

check_parameters = function(myresults){
  
  labels = as.factor(sapply(myresults, function(x) ifelse(!x$empty, x$label, NA)))
  nfeatures.perm = sapply(myresults, function(x) colMeans(x$nfeatures.perm))
  K.iter = sapply(myresults, function(x) x$K.iter)
  eta.iter = sapply(myresults, function(x) x$eta.iter)
  colnames(nfeatures.perm) = colnames(K.iter) = colnames(eta.iter) = labels
  
  # motion parameters
  print("Model parameters")
  
  print(ggplot(melt(K.iter), aes(x = value)) + geom_density() + facet_grid(. ~ Var2))
  print(ggplot(melt(eta.iter), aes(x = value)) + geom_density() + facet_grid(. ~ Var2))
  print(ggplot(melt(nfeatures.perm), aes(x = Var1, y = value)) + geom_line() + geom_point() + facet_grid(. ~ Var2))
  
  # motion parameters / same in all modules
  print("Motion parameters")
  fd = myresults[[1]]$covars.exp$fd  
  meica_dof = myresults[[1]]$covars.exp$meica_dof
  par(mfrow = c(1, 2))
  hist(fd, 30)
  hist(meica_dof, 30)
  print(summary(fd))
  print(summary(meica_dof))
}

statstodf = function(myresults, group_name, approx = F, check_pars = T, groups = NULL){
  print("--------------------------")
  if (check_pars)  check_parameters(myresults)  
  alpha = 0.05
  module_names = sapply(myresults, function(x) x$label)
  if (!approx)
    rho = data.frame(
      mean = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test$estimate, NA)), 
      #      low = sapply(myresults, function(x) ifelse(!x$empty, quantile(x$cor.test$perm, 0.025), NA)),
      #      high = sapply(myresults, function(x) ifelse(!x$empty, quantile(x$cor.test$perm, 0.975), NA)),
      low = sapply(myresults, function(x) ifelse(!x$empty, min(x$cor.test$perm), NA)),
      high = sapply(myresults, function(x) ifelse(!x$empty, quantile(x$cor.test$perm, 0.95), NA)),
      p = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test$p.value, NA)),
      nfeatures = sapply(myresults, function(x) ifelse(!x$empty, x$nfeatures, NA)),
      from = as.factor(sapply(myresults, function(x) ifelse(!x$empty, x$label, NA))),
      time = sapply(myresults, function(x) ifelse(!x$empty, x$time, NA))/3600
    )        
  else
    rho = data.frame(
      mean = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.approx$estimate, NA)), 
      low = sapply(myresults, function(x) ifelse(!x$empty, quantile(x$cor.test.approx$perm, 0.025), NA)),
      high = sapply(myresults, function(x) ifelse(!x$empty, quantile(x$cor.test.approx$perm, 0.975), NA)),
      p = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.approx$p.value, NA)),
      nfeatures = sapply(myresults, function(x) ifelse(!x$empty, x$nfeatures, NA)),
      from = as.factor(sapply(myresults, function(x) ifelse(!x$empty, x$label, NA))),
      time = sapply(myresults, function(x) ifelse(!x$empty, x$time, NA))/3600
    )        
  rho$group_name = group_name
  
  cor.perm = sapply(myresults, function(x) x$cor.test$perm)
  R2.perm = sapply(myresults, function(x) x$R2.test$perm)
  
  # get perm data
  if (!is.null(groups))
  {
    rho.groups = NULL
    cors.perm.group = list()
    for (g in groups){
      print(g)
      rho.group = rho
      rho.group$group_name = g
      for (j in seq(length(myresults))){
        result = myresults[[j]]
        cor.perm.group = apply(result$y.pred.perm[ result$group == g, ], 2, function(x) cor.test(x, result$y.test[ result$group == g])$estimate)
        rho.group$mean[j] = cor.perm.group[1] 
        rho.group$high[j] = quantile(cor.perm.group, 0.975)
        rho.group$low[j] = quantile(cor.perm.group, 0.025)
        rho.group$p[j] = sum(cor.perm.group[1] <= cor.perm.group)/length(cor.perm.group)
        
        cors.perm.group[[g]] = rbind(cors.perm.group[[g]], cor.perm.group)
      }
      rho.groups = rbind(rho.groups, rho.group) 
    }
    #    rho = rho.groups
    r.bsl = subset(rho.groups, group_name == "bsl" & from != "All")$mean
    r.fu = subset(rho.groups, group_name == "fu" & from != "All")$mean
    
    print(cor.test(r.bsl, r.fu))
    
    r.perm = NULL
    for (j in seq(length(cor.perm.group))) r.perm[j] = cor.test(cors.perm.group[['bsl']][-1, j], cors.perm.group[['fu']][-1, j], method = "spearman")$estimate 
    
    pvalue = sum(r.perm[1] <= r.perm)/length(r.perm)
    
    save_fig(figname = "GroupCorrelation", res = BWRES)
    
    par(mfrow=c(1,1))
    plot(r.bsl, r.fu, xlab = "Baseline r", ylab = "Follow-up r", pch = 20, col = "white" )
    text(r.bsl, r.fu,
         labels = subset(rho.groups, group_name == "bsl" & from != "All")$from)
    abline(0, 1, lty = 2)
    legend("topleft", legend = paste("r=", format(r.perm[1], digits = 2), 
                                     "p=", pvalue))
    dev.off()
    
  }
  #browser()
  # get max statistic and corresponding p-values
  nullmax = apply(cor.perm, 1, max)
  rho$pvalues.FWE = sapply(cor.perm[1, ], function(x) sum(x <= nullmax)/length(nullmax))
  
  # check the prediction is not driven by the number of features
  print("Correlation between number of features and mean rho (abs and original)")
  print(cor.test(rho$nfeatures, abs(rho$mean)))  
  print(cor.test(rho$nfeatures, rho$mean), method = "spearman")#, alternative = "greater"))  
  par(mfrow=c(1,1))
  plot(rho$nfeatures, rho$mean)
  
  rho$width = ifelse(rho$p < alpha, 1, 0.5)
  rho$size = ifelse(rho$p < alpha, 5, 4)
  
  rho$network = as.factor(rho$from) 
  ps = rho$p
  ps.adjusted = p.adjust(ps, method = "fdr")
  print("Adjusted correlation all features")
  rho$fdrp = ps.adjusted 
  rho$sig = ifelse(rho$p < alpha, "*",  "")
  #rho$sig[rho$fdrp < alpha] = "**"
  rho$sig[rho$pvalues.FWE < alpha] = "**"
  
  rho = rho %>% arrange(network)
  print(rho)
  
  mytimeplot =  
    ggplot(rho, aes( x = network, y = time )) + 
    geom_col() + 
    theme_bw() + 
    xlab('\nModule') + 
    ylab('Time') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  print(mytimeplot)
  
  print(paste("Total time", sum(rho$time)))
  return(rho)
}
plot_combinations = function(rho, measure, myfigname = 'Accuracy', group = F){

  if (group) {
    myplot =  
      ggplot(rho, aes( x = network, y = mean, ymin = low, ymax = high, label = rho$sig, col = group_name)) +
      geom_errorbar(position = position_dodge(0.4), size = 2, linetype = 1, width = 0.4) + 
      ggtitle(paste('Correlation between observed and predicted', measure, "\n")) + xlab('Network') + 
      geom_point(position = position_dodge(0.4), size = 15, shape = 19) + 
      geom_text(aes( x = network, y = mean + .02), position = position_dodge(0.4), size = CEX_TEXT) +
      ylim(-.3, .3) + 
      theme_bw() + 
      xlab('\nModule') + 
      ylab('r') +
      theme(
        legend.position = "none",
        #legend.text = element_text(size = CEX_AXIS - 10),
        axis.text = element_text(size = CEX_AXIS - 10),
        axis.text.x = element_text(angle = 45, hjust = 1),# vjust = -1.5, 
        axis.text.y = element_text(hjust = -.2),
        axis.title = element_text(size = CEX_AXIS, vjust = 2),
        plot.title = element_text(hjust = 0.5, vjust = -5, size = CEX_TITLE - 10),
        #        plot.margin = margin(1, 0, 2, 0, "cm"),
        axis.line = element_line(size = 2),
        axis.ticks.length = unit(20, "pt"),
        axis.ticks = element_line(size = 2)
      )
  } else {
    myplot = 
      ggplot(rho, aes( x = network, y = mean, ymin = low, ymax = high, label = rho$sig)) + 
      geom_errorbar(position = position_dodge(0.4), size = 2, linetype = 1, width = 0.4) + 
      ggtitle(paste('Correlation between observed and predicted', measure, "\n")) + xlab('Network') + 
      geom_point(position = position_dodge(0.4), size = 15, shape = 19) + 
      geom_text(aes( x = network, y = mean + .02), position = position_dodge(0.4), size = CEX_TEXT) +
      ylim(-.3, .3) + 
      theme_bw() + 
      xlab('\nModule') + 
      ylab('r') +
      theme(
        legend.position = "none",
        #      legend.text = element_text(size = CEX_AXIS - 10),
        axis.text = element_text(size = CEX_AXIS - 10),
        axis.text.x = element_text(angle = 45, hjust = 1),# vjust = -1.5, 
        axis.text.y = element_text(hjust = -.2),
        axis.title = element_text(size = CEX_AXIS, vjust = 2),
        plot.title = element_text(hjust = 0.5, vjust = -5, size = CEX_TITLE - 10),
        #        plot.margin = margin(1, 0, 2, 0, "cm"),
        axis.line = element_line(size = 2),
        axis.ticks.length = unit(20, "pt"),
        axis.ticks = element_line(size = 2)
      )
  }
  save_fig(figname = myfigname, res = BWRES)
  print(myplot)
  dev.off()
  
}

plot_perms = function(RESULTS_FILE, measure, myfigname, selection = NULL){

  load(RESULTS_FILE)
  if (!is.null(selection)){
    mylist = NULL
    j = 1
    for (res in var.network.results) if (res$label %in% selection) {
      mylist[[j]] = res
      j = j + 1
    }
  var.network.results = mylist
  }
  
  rho = statstodf(var.network.results, group_name  = "") %>% arrange(network) 
  
  cor.perm = sapply(var.network.results, function(x) x$cor.test$perm)
  colnames(cor.perm) = sapply(var.network.results, function(x) x$label)
  cor.perm.melt = melt(cor.perm, value.name = "cor", varnames = c("permutation", "network")) %>% 
    mutate(network = factor(network, levels = rho$network)) %>% arrange(network)
  myplot = 
    ggplot(data = cor.perm.melt, aes(x = network, y = cor)) +
    geom_violin(adjust = 0.5, size = 2) +
    geom_point(data = rho, aes( x = network, y = mean), position = position_dodge(0.4), size = 15, shape = 19) +
    geom_text(data = rho, aes( x = network, y = mean + .02, label = rho$sig), position = position_dodge(0.4), size = CEX_TEXT) +
    #  geom_jitter(height = 0, width = 0.1) +
    ggtitle(paste('Correlation between observed and predicted', measure, "\n")) + xlab('Network') + 
    ylim(-.25, .25) + 
    theme_bw() + 
    xlab('\nModule') + 
    ylab('r') +
    theme(
      legend.position = "none",
      #legend.text = element_text(size = CEX_AXIS - 10),
      axis.text = element_text(size = CEX_AXIS - 10),
      axis.text.x = element_text(angle = 45, hjust = 1),# vjust = -1.5, 
      axis.text.y = element_text(hjust = -.2),
      axis.title = element_text(size = CEX_AXIS, vjust = 2),
      plot.title = element_text(hjust = 0.5, vjust = -5, size = CEX_TITLE - 10),
      #        plot.margin = margin(1, 0, 2, 0, "cm"),
      axis.line = element_line(size = 2),
      axis.ticks.length = unit(20, "pt"),
      axis.ticks = element_line(size = 2)
    )
  
  save_fig(figname = myfigname, res = BWRES)
  print(myplot)
  dev.off()
  return(rho)
}



###############################################
# Fit the Bivariate Latent Change Score model 
###############################################

fitBLCSmodel = function(datBLCS){
  
  BLCS<-'

COG_T2 ~ 1*COG_T1     # This parameter regresses COG_T2 perfectly on COG_T1
dCOG1 =~ 1*COG_T2     # This defines the latent change score factor as measured perfectly by scores on COG_T2
dCOG1 ~ 1             # This estimates the intercept of the change score 
COG_T1 ~  1           # This estimates the intercept of COG_T1 
COG_T2 ~ 0*1          # This constrains the intercept of COG_T2 to 0

NEU_T2 ~ 1*NEU_T1     # This parameter regresses NEU_T2 perfectly on NEU_T1
dNEU1 =~ 1*NEU_T2     # This defines the latent change score factor as measured perfectly by scores on NEU_T2
NEU_T2 ~ 0*1          # This line constrains the intercept of NEU_T2 to 0
NEU_T2 ~~ 0*NEU_T2    # This fixes the variance of the NEU_T1 to 0  

dCOG1 ~~  dCOG1       # This estimates the variance of the change scores
COG_T1 ~~   COG_T1    # This estimates the variance of the COG_T1 
COG_T2 ~~ 0*COG_T2    # This fixes the variance of the COG_T2 to 0  

dNEU1 ~ 1             # This estimates the intercept of the change score 
NEU_T1 ~ 1            # This estimates the intercept of NEU_T1 
dNEU1 ~~ dNEU1        # This estimates the variance of the change scores 
NEU_T1 ~~ NEU_T1      # This estimates the variance of NEU_T1 

dNEU1~COG_T1+NEU_T1   # This estimates the COG to NEU coupling parameter and the COG to COG self-feedback
dCOG1~NEU_T1+COG_T1   # This estimates the NEU to COG coupling parameter and the NEU to NEU self-feedback
#dNEU1~NEU_T1   # This estimates the COG to NEU coupling parameter and the COG to COG self-feedback
#dCOG1~COG_T1   # This estimates the NEU to COG coupling parameter and the NEU to NEU self-feedback

COG_T1 ~~  NEU_T1     # This estimates the COG_T1 NEU_T1 covariance
dCOG1~~dNEU1          # This estimates the dCOG and dNEU covariance
'
  #browser()
  
  fitBLCS <- lavaan(BLCS, data=datBLCS[c("COG_T1", "COG_T2", "NEU_T1", "NEU_T2")], estimator='mlr',fixed.x=FALSE,missing='fiml')
  
  #summary(fitBLCS, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
  
  # collect results
  results = rbind(
    subset(parameterEstimates(fitBLCS), lhs == "COG_T1" & rhs ==  "NEU_T1"),
    subset(parameterEstimates(fitBLCS), lhs == "dCOG1" & rhs ==  "dNEU1")
  )
  
  results$type = c("level", "change")
  results$label = datBLCS$label[1]
  return(results)
}