
###############################################################################
# Several utility functions for prediction
###############################################################################

# install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(mclust, pryr, spls, foreach, doParallel, doMC, caret, pracma, 
               dplyr, RColorBrewer, scales, rlist, abind)

#library(RColorBrewer)
#library(scales)
#library(lavaan)
#library(RVAideMemoire)
#library(psychometric)
#library(boot)
#library(psych)
#library(rjson)

FOLDSEED = sample(seq(10000), 1)
MYPALETTE=brewer.pal(11, "Spectral")[seq(11, 1, -1)]
FINELINE = 0.5
GGTEXTSIZE2 = 10
GGTEXTSIZE3 = 10
NDIGITS = 3

NPERMSAPPROX = 10000


###############################################################################
# Do prediction for different networks
###############################################################################

do_predictions_loop = function(OUTPUT_DIR,
                               DIR_PREFFIX, 
                               DIR_ICA_ME, # working dir
                               data.pars, # predictors
                               var, # independent variable
                               mycovars, # covariates
                               modules_file,
                               to_select_list,
                               maxcomp = DEFAULTMAXCOMP, # max spls components 
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

  RESULTS_DIR = file.path(OUTPUT_DIR, paste0(DIR_PREFFIX, var))
  dir.create(RESULTS_DIR)
  RESULTS_FILE = file.path(RESULTS_DIR, "data.RData")
  
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
                                 savecoefs = ifelse(savecoefs, RESULTS_DIR, '')
)
    # View(results)
    print(paste("Done", i/length(to_select_list)*100, "%"))
    
  }
  
  results.df = NULL
  for (j in 1:length(results)) {
    
    if(to_select_list[[j]][1] == 0){
      results[[j]]$label = "All"
    }
    else {
      results[[j]]$label = module_labels[module_names == to_select_list[[j]][1]]
    }
    results.df = rbind(results.df, 
                       resultstodf(results[[j]], to_select_list[[j]][1])
    )
  }
  save(results.df, results, file = RESULTS_FILE)

  return(results)
}


###############################################################################
# Fit and test predictive model, main function
###############################################################################

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
                         mymodel = NULL, savecoefs = '', 
                         N_REPS = 5 # number of repetition of xvalidation
)
{ 
  print("-----------------------------------------------------------------")
  print(var)
  
  tic("model fitting")
  results = list()
  results$var = var
  results$empty = F
  
  # needed files
  FILE_ICA_ME = file.path(DIR_ICA_ME, netmats_file) 
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
    colnames(covars) = c("group", "Subject","dvars", "fd","brain_vol",
                         "meica_dof")
    mydata = merge(mydata, covars, by = c("Subject", "group"))
    
  } else {
    fd = cbind(fd[1], rowMeans(fd[-c(1,2)]))
    dvars = cbind(fd[1], rowMeans(dvars[-c(1,2)]))
    colnames(ica_data) = paste0("feature.", seq(ncol(ica_data)))
    ica_data$Subject = subjects
    
    mydata = merge(ica_data, data.pars, by="Subject", all.x = T)
    
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
  
  mydata.exp = mydata.exp[complete.cases(mydata.exp), ]
  
  if (!is.null(fdfiltering)){
    print(paste("Removing", sum(mydata.exp$fd > fdfiltering), 
                "subjects with high mean FD above", fdfiltering ))
    results$fdfiltered = subset(mydata.exp, fd > fdfiltering)$Subject
    mydata.exp = subset(mydata.exp, fd <= fdfiltering)
  }
  

  y = unlist(mydata.exp[var])
  
  #remove outliers if any
  valid.rows = (y - mean(y)) < 3*sd(y)
  print(paste("Removed", sum(!valid.rows), "outliers of", var))
  
  covars.exp = mydata.exp[valid.rows, mycovars]
  if ("UCL.1" %in% colnames(covars.exp)) covars.exp = covars.exp %>% 
    rename(UCL = UCL.1, CBU = CBU.1) # ugly hack...
  if ("UCL.2" %in% colnames(covars.exp)) covars.exp = covars.exp %>% 
    rename(UCL = UCL.2, CBU = CBU.2) # ugly hack...
  
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
      folds = c(folds, list(Fold = which(mydata.exp$Subject %in% 
                                           unique_subjects[folds_sub[[i]]])))
    }
  }
  else folds = multiFolds(y, k = N_FOLDS, times = N_REPS) 
  #createFolds(y, k = N_FOLDS)

  mysubjects = mydata.exp$Subject[valid.rows]
  mygroup = mydata.exp$group[valid.rows]
  
  if (is.null(mymodel)){
    
    cl <- makeCluster(NPROCS)
    registerDoParallel(cl)
    cv = foreach(fold_index = seq(length(folds))) %do% 
      do_crossvalidate_spls_covars_perm_par(fold_index, 
                                            list(X = ica_data, 
                                                 y = y, 
                                                 covars = covars.exp, 
                                                 folds = folds, 
                                                 subject = mysubjects,
                                                 group = mygroup), 
                                            maxcomp = maxcomp, cluster = cl, 
                                            NITER = NITER, NPERM = NPERM)
    
    stopCluster(cl)

    # collect results
    results$N_REPS = N_REPS
    results$projected = F  
    
    results$fold = unlist(lapply(cv, function(x){x$fold}))
    results$training_subjects = lapply(cv, function(x){mysubjects[-x$fold]})
    results$testing_subjects = lapply(cv, function(x){mysubjects[x$fold]})
    results$y.pred.perm = list.rbind(lapply(cv, function(x){x$y.pred.perm}))  
    results$y.test.orig = unlist(lapply(cv, function(x){x$y.test.orig}))
    results$y.pred = unlist(lapply(cv, function(x){x$y.pred}))
    results$y.test = unlist(lapply(cv, function(x){x$y.test}))
    
    results$nfeatures.perm = list.rbind(lapply(cv, 
                                               function(x){x$nfeatures.perm}))  
    results$mysample.perm = lapply(cv, 
                                   function(x){as.integer(x$mysample.perm)})  
    results$eta.iter = unlist(lapply(cv, function(x){x$eta.iter}))      
    results$K.iter = unlist(lapply(cv, function(x){x$K.iter}))  
    results$offset = unlist(lapply(cv, function(x){x$offset}))  
    results$coefs.cv = do.call("cbind", lapply(cv, function(x){x$coefs}) )
    
    results$coefs.mean = rowMeans(results$coefs.cv)
    results$coefs.stable = apply(results$coefs.cv, 1,
                                 function(x) 
                                   ifelse( prod(range(x)) < 0, 0, mean(x)))
    

    # save these to disk to spare memory
    if ( savecoefs != '') {
    coefs.perm = do.call("abind", 
                         list(lapply(cv, function(x){x$coefs.perm}), along = 3))
    results$coefs.perm.file = tempfile(pattern = "coefs", 
                                       tmpdir = savecoefs, fileext = ".rda") 
    save(coefs.perm, file = results$coefs.perm.file) 
    
    results$preprocessing.file = tempfile(pattern = "preproc", 
                                          tmpdir = savecoefs, fileext = ".rda")
    preprocessing = lapply(cv, function(x) x$preprocessing) 
    save(preprocessing, file = results$preprocessing.file) 
    }
    
  } else{
    
    # do not estimate, apply previously computed model
    results$projected = T 
    results$fold = NULL
    results$nfeatures.perm = mymodel$nfeatures.perm  
    results$eta.iter = mymodel$eta.iter  
    results$K.iter = mymodel$K.iter
    
    load(mymodel$coefs.perm.file)
    mymodel$coefs.perm = coefs.perm

    load(mymodel$preprocessing.file)
    mymodel$preprocessing = preprocessing
    
    results$coefs.mean = mymodel$coefs.mean
    results$coefs.stable = mymodel$coefs.stable

    results.fit = apply_model(ica_data, y, covars.exp, mymodel, 
                              mysubjects, N_REPS, N_FOLDS)
    results$fold = results.fit$fold
    results$y.test = results.fit$y.test
    results$mysubjects.test = results.fit$mysubjects.test
    results$y.test.orig = results.fit$y.test.orig
    results$y.pred.perm = results.fit$y.pred.perm
    results$y.pred = results.fit$y.pred.perm[, 1]

  }

  results$covars.exp = covars.exp
  results$demo.orig = mydata.exp[valid.rows, c("age", "sex.num")]
  results$demo = results$demo.orig[results$fold, ]
  results$nfeatures = ncol(ica_data)
  results$Subject.orig = mysubjects
  results$Subject = mysubjects[results$fold]
  
  if (longformat) results$group = mydata.exp$group

  results$cor.test = results.R2.test = results$cor.test.approx = NULL
  results$nobs = nrow(ica_data)
  
  # average the coefficients across reps
  myrep = unlist(lapply(seq(N_REPS), function(x) rep(x, results$nobs)))
  cor.perm = NULL
  for (k in seq(N_REPS)){
    cor.perm = 
      rbind(cor.perm, apply(results$y.pred.perm[myrep == k, ], 2, 
      function(x) cor.test(x, results$y.test[myrep == k])$estimate))
  }
  
  results$cor.test$perm = fisherz2r(colMeans(fisherz(cor.perm)))
  results$cor.test$estimate = results$cor.test$perm[1]
  results$cor.test$p.value = 
    sum(results$cor.test$perm[1] <= 
          results$cor.test$perm)/length(results$cor.test$perm)
  
  # approximateby fitting a Gaussian
  z.perm = fisherz(results$cor.test$perm)
  mysigma = sqrt(sum(z.perm[-1]^2)/length(z.perm[-1]))
  results$cor.test$p.value.approx = 1 - pnorm(z.perm[1]/mysigma)
  
  # so that one can also show null distributions
  results$cor.test.approx$perm = 
    fisherz2r(rnorm(NPERMSAPPROX)/sqrt(length(results$y.test) - 3))
  results$cor.test.approx$perm[1] = results$cor.test$perm[1]
  results$cor.test.approx$estimate = results$cor.test$perm[1]
  results$cor.test.approx$p.value.perm = 
    sum(results$cor.test$perm[1] <= 
          results$cor.test.approx$perm)/length(results$cor.test.approx$perm)
  results$cor.test.approx$p.value = 
    cor.test(results$y.pred, results$y.test)$p.value
  
  results$R2.test$perm = 
    apply(results$y.pred.perm, 2, function(x) r2score(x, results$y.test))
  results$R2.test$estimate = results$R2.test$perm[1]
  results$R2.test$p.value = 
    sum(results$R2.test$perm[1] <= 
          results$R2.test$perm)/length(results$R2.test$perm)
  
  print(paste0("Estimate:", results$cor.test$estimate, ",  
               exact p: ", results$cor.test$p.value, ", 
               approx p: ",  results$cor.test$p.value.approx))

  results$time = toc()

  print(paste0("Memory used:", mem_used()/10^9, " Gb"))
  print(paste0("Tempdir:", tempdir()))
  return(results)
}

###############################################################################
# Apply previously learned model to new dataset
###############################################################################

apply_model = function(X, y, covars, model, mysubjects, nrep, nfold){
  
  # for each fold, find model and average across folds
  fold = NULL
  y.test.mat = y.pred.mat = mysubjects.test = NULL
  
  for (i in seq(length(model$preprocessing))){
    trainers = as.vector(na.omit(match(model$training_subjects[[i]], 
                                       mysubjects)))
    trainers = trainers[!is.na(trainers)]
    testers = setdiff(mysubjects, trainers)
    fold = c(fold, as.vector(na.omit(match(testers, mysubjects))))
    X.test = X[-trainers, ]
    mysubjects.test[[i]] = mysubjects[-trainers]   
    
    if (!is.null(covars)) 
      y.test = y - predict(model$preprocessing[[i]]$mycovarmod, covars)
    else y.test = y
    y.test[trainers] = NA
    
    mu = model$preprocessing[[i]]$mu
    sigma = model$preprocessing[[i]]$sigma
    
    for (j in seq(ncol(X))){
      X.test[, j] = (X.test[, j] - mu[j])/sigma[j]     
    }
    
    y.pred.perm = matrix(NA, nrow(X), ncol(model$coefs.perm[, , i]))
    y.pred.perm[-trainers, ] = 
      X.test %*% model$coefs.perm[, , i] + model$offset[i] 
    y.test.mat = cbind(y.test.mat, y.test)
    y.pred.mat = abind(y.pred.mat, y.pred.perm, along = 3)
    
  }
  y.test = y.pred.perm = NULL 
  myrep = unlist(lapply(seq(nrep), function(x) rep(x, nfold)))
  
  for (l in unique(myrep)) {
    y.test = c(y.test, rowMeans(y.test.mat[, myrep == l], na.rm = T))
    y.pred.perm = 
      rbind(y.pred.perm, apply(y.pred.mat[, , myrep == l], c(1, 2), mean, na.rm = T))
  }
  
  return(list(y.test = y.test, 
              y.pred.perm = y.pred.perm, 
              y.test.orig = y,
              testing_subjects = mysubjects.test,
              fold = fold)
  )
}


###############################################################################
# Select connections in a module
###############################################################################
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


###############################################################################
# Generate folds for repeated k-fold cross-validation
###############################################################################
multiFolds = function (y, k = 10, times = 5) 
{
  if (class(y)[1] == "Surv") 
    y <- y[, "time"]
  prettyNums <- paste("Rep", gsub(" ", "0", format(1:times)), 
                      sep = "")
  for (i in 1:times) {
    tmp <- createFolds(y, k = k, list = TRUE, returnTrain = F)
    names(tmp) <- paste("Fold", gsub(" ", "0", format(seq(along = tmp))), 
                        ".", prettyNums[i], sep = "")
    out <- if (i == 1) 
      tmp
    else c(out, tmp)
  }
  out
}

###############################################################################
# Compute r2
###############################################################################
r2score = function(x, y){
  
  y = y - mean(y) # we don't care about the offset
  z = lm(y ~ x)$residuals
  r2 = 1 - var(z)/var(y) 
  return(r2)
}

###############################################################################
# Squareform modified version
###############################################################################
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

###############################################################################
# Convert list of results to data frame
###############################################################################

resultstodf = function(results, number){
  
  results$demo$sex.num = NULL
  results.df = cbind(results[ c("Subject", "y.pred", "y.test")], 
                     results$covars.exp, results$demo)
  results.df$var = results$var
  results.df$number = number
  results.df$label = results$label
  
  if (!is.null(results$group)) results.df$group = results$group
  
  return(as.data.frame(results.df))
}
