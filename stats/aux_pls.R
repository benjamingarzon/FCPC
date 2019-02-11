######## PLS methods


do_bagging_RFE = function(X.train, y.train, n_folds, alphafeat, maxcomp = 10, selected = NULL, n_steps = 20){
  
  
  # univeariate feature filter
  #  cs = apply(X.train, 2, function(x) cor.test(y.train, x, method="spearman"))
  cs = apply(X.train, 2, function(x) cor.test(y.train, x))
  ps = sapply(cs, function(x) x$p.value)
  rs = sapply(cs, function(x) x$estimate)
  
  if (is.null(selected)) {
    selected = (ps < alphafeat)
  } else {
    selected = (ps < alphafeat) | selected
  }
  
  # calculate folds and steps
  folds = createFolds(y.train, k = n_folds, list = TRUE, returnTrain = FALSE)
  
  #steps = rep(.9, n_steps) 
  steps = seq(1 - 1/(n_steps + 1), 0, -1/(n_steps + 1))
  n_init = sum(selected)
  
  #load("aux_pls_image.RData")
  params = expand.grid(mystep = seq(length(steps)), comp = seq(maxcomp), RMSEP = 1)
  #myparam = 1
  RMSEP.min = NULL
  comp_data = list()
  for (comp in seq(maxcomp)){
    print(comp)
    selected.comp = selected
    myRMSEP.mean = myselected = mycoefs.mean.mat = Xmeans.mat = Ymeans.mat = NULL
    
    for (mystep in seq(length(steps))){
      Xmeans = Ymeans =  mycoefs.mat = myRMSEP = NULL  
      
      maxcomp = min(maxcomp, sum(selected.comp))
      
      for (n_fold in seq(n_folds)){
        fold = folds[[n_fold]]
        mypls = plsr(y.train[-fold] ~ X.train[-fold, selected.comp], ncomp = comp)
        myRMSEP[n_fold] = sqrt(sum((y.train[fold] - predict(mypls, X.train[fold, selected.comp]))**2))
        
        Ymeans[n_fold] = mypls$Ymean
        mycoefs = Xmean = ps *0
        mycoefs[selected.comp] = coef(mypls, ncomps = comp)
        mycoefs.mat = rbind(mycoefs.mat, mycoefs) 
        Xmean[selected.comp] = mypls$Xmean
        Xmeans = rbind(Xmeans, Xmean)
      }
      
      #mypls = plsr(y.train ~ X.train[, selected], ncomp = comp)
      #Xmeans = mypls$Xmean
      #Ymeans = mypls$Ymean
      #mycoefs = ps * 0
      #mycoefs[selected.comp] = coef(mypls, ncomps = comp)
      #mycoefs.mat = mycoefs 
      
      myRMSEP.mean[mystep] = mean(myRMSEP)
      #    params$RMSEP[myparam] = mean(myRMSEP)
      
      # rank features and select a fraction of them    
      mycoefs.t = apply(mycoefs.mat, 2, function(x) t.test(x)$statistic)
      mycoefs.t[is.na(mycoefs.t)] = 0
      
      mycoefs.mean.mat = rbind(mycoefs.mean.mat, colMeans(mycoefs.mat))
      Xmeans.mat = rbind(Xmeans.mat, colMeans(Xmeans))
      Ymeans.mat = c(Ymeans.mat, mean(Ymeans))
      nfeatures = ceiling(steps[mystep] * n_init)
      valid = order(abs(mycoefs.t), decreasing = T)[1:nfeatures]
      myselected = rbind(myselected, selected.comp)
      selected.comp[selected.comp] = F
      selected.comp[valid] = T
      #    myparam = myparam + 1
    }
    
    RMSEP.min[comp] = min(myRMSEP.mean)
    RMSEP.which.min = which.min(myRMSEP.mean)
    comp_data[[comp]] = list()
    comp_data[[comp]]$RMSEP.which.min = which.min(myRMSEP.mean)
    comp_data[[comp]]$Ymean = Ymeans.mat[RMSEP.which.min] 
    comp_data[[comp]]$mycoefs.mean = mycoefs.mean.mat[RMSEP.which.min, ]
    comp_data[[comp]]$Xmean = Xmeans.mat[RMSEP.which.min, ]
    comp_data[[comp]]$nonzero = myselected[RMSEP.which.min, ]
    comp_data[[comp]]$nfeatures = rowSums(myselected)
  } # components
  
  
  which.comp.min  = which.min(RMSEP.min)
  nonzero =  comp_data[[which.comp.min]]$nonzero
  mycoefs.mean = comp_data[[which.comp.min]]$mycoefs.mean
  Ymean = comp_data[[which.comp.min]]$Ymean
  Xmean = comp_data[[which.comp.min]]$Xmean
  
  #plot(nfeatures, myRMSEP.mean)
  #ggplot(params, aes(comp, mystep)) + geom_tile(aes(fill = RMSEP))
  #ggplot(params, aes(x=mystep, y = RMSEP, col = comp, group = comp)) + geom_line()                                              
  nfeatures = sum(nonzero)
  comps = which.comp.min
  return(list(mycoefs.mean = mycoefs.mean, Ymean = Ymean, Xmean = Xmean, comps = comps, nfeatures = nfeatures))
}


do_bagging_RFE_1 = function(X.train, y.train, n_folds, alphafeat, maxcomp = 10, selected = NULL, n_steps = 20){
  
  # univeariate feature filter
  cs = apply(X.train, 2, function(x) cor.test(y.train, x))
  ps = sapply(cs, function(x) x$p.value)
  rs = sapply(cs, function(x) x$estimate)
  
  if (is.null(selected)) {
    selected = (ps < alphafeat)
  } else {
    selected = (ps < alphafeat) | selected
  }
  
  # calculate folds and steps
  folds = createFolds(y.train, k = n_folds, list = TRUE, returnTrain = FALSE)
  
  steps = seq(1 - 1/(n_steps + 1), 0, -1/(n_steps + 1))
  
  n_init = sum(selected)
  
  myRMSEP.mean = myselected = mycoefs.mean.mat = Xmeans.mat = Ymeans.mat = comps.mat = NULL
  
  for (mystep in seq(length(steps))){
    #print(mystep)
    mycoefs.mat = NULL
    myRMSEP = comps = NULL  
    Xmeans = Ymeans = NULL
    maxcomp = min(maxcomp, sum(selected))
    
    
    
    for (n_fold in seq(n_folds)){
      fold = folds[[n_fold]]
      mypls = plsr(y.train[-fold] ~ X.train[-fold, selected], ncomp = maxcomp, validation = "CV")
      comp = selectNcomp(mypls, method="randomization")    
      #      mypls = plsr(y.train[-fold] ~ X.train[-fold, selected], ncomp = 2, validation = "CV")
      #      comp = 2
      mycoefs = Xmean = ps * 0
      mycoefs[selected] = coef(mypls, ncomps = comp)
      mycoefs.mat = rbind(mycoefs.mat, mycoefs)
      Xmean[selected] = mypls$Xmean
      Xmeans = rbind(Xmeans, Xmean)
      Ymeans[n_fold] = mypls$Ymean
      comps[n_fold] = comp
      # calculate error
      myRMSEP[n_fold] = sqrt(sum((y.train[fold] - predict(mypls, X.train[fold, selected]))**2))
      
    }
    
    myRMSEP.mean[mystep] = mean(myRMSEP)
    
    # rank features and select a fraction of them    
    mycoefs.t = apply(mycoefs.mat, 2, function(x) t.test(x)$statistic)
    mycoefs.t[is.na(mycoefs.t)] = 0
    
    mycoefs.mean.mat = rbind(mycoefs.mean.mat, colMeans(mycoefs.mat))
    Xmeans.mat = rbind(Xmeans.mat, colMeans(Xmeans))
    Ymeans.mat = c(Ymeans.mat, mean(Ymeans))
    comps.mat = c(comps.mat, mean(comps))
    nfeatures = ceiling(steps[mystep] * n_init) 
    valid = order(abs(mycoefs.t), decreasing = T)[1:nfeatures]
    myselected = rbind(myselected, selected)
    selected[selected] = F
    selected[valid] = T
  }
  
  # get the optimal number of features and model coefficients
  nfeatures = rowSums(myselected)
  RMSEP.min = which.min(myRMSEP.mean)
  nonzero = myselected[RMSEP.min, ]
  
  mycoefs.mean = mycoefs.mean.mat[RMSEP.min, ]
  Ymean = Ymeans.mat[RMSEP.min]
  Xmean = Xmeans.mat[RMSEP.min,]
  #plot(nfeatures, myRMSEP.mean)
  
  nfeatures = sum(nonzero)
  comps = comps.mat[RMSEP.min]
  return(list(mycoefs.mean = mycoefs.mean, Ymean = Ymean, Xmean = Xmean, comps = comps, nfeatures = nfeatures))
}

do_crossvalidate_pls_diff = function(fold, input, covars, alphafeat = 0.1, NITER = 20, NSTEPS = 20){
  
  # compare with NULL model (only covariates)
  cv = do_crossvalidate_pls(fold, input, covars, alphafeat = alphafeat, NITER = NITER, NSTEPS = NSTEPS)
  input$X = NULL
  # include all covariates
  cv.0 = do_crossvalidate_pls(fold, input, covars, alphafeat = 1.1, NITER = NITER, NSTEPS = NSTEPS)
  
  cv$rho.diff = cv$rho - cv.0$rho # expected larger
  cv$RMSE.diff = cv$RMSE - cv.0$RMSE  # expected smaller
  cv$rho.0 = cv.0$rho 
  cv$RMSE.0 = cv.0$RMSE  
  
  return(cv)
}

do_crossvalidate_pls = function(fold, input, covars, alphafeat = 0.1, NITER = 20, NSTEPS = 20){
  
  # selecting all features by default
  
  X = input$X
  y = input$y
  
  # join fMRI and covars as features
  if(is.null(X)){
    X.train = as.matrix(covars[-fold, , drop=FALSE])
    X.test = as.matrix(covars[fold, , drop=FALSE])
    selected = rep(T, ncol(covars)) 
  } else {
    if(is.null(covars)){
      X.train = as.matrix(X[-fold, , drop=FALSE])
      X.test = as.matrix(X[fold, , drop=FALSE])
      selected = rep(F, ncol(X)) 
    } else {
      X.train = as.matrix(cbind(X[-fold, , drop=FALSE], covars[-fold, , drop=FALSE]))
      X.test = as.matrix(cbind(X[fold, , drop=FALSE], covars[fold, , drop=FALSE]))
      selected = c(rep(F, ncol(X)), rep(T, ncol(covars)))   
    }
  }
  y.train = y[-fold]  
  y.test = y[fold]
  
  # demean and scale
  for (j in seq(ncol(X.train))){
    
    mu = mean(X.train[, j])
    sigma = sd(X.train[, j])
    X.train[, j] = (X.train[, j] - mu)/sigma     
    X.test[, j] = (X.test[, j] - mu)/sigma    
  }
  maxcomp = 10
  # get model
  if (is.null(X)){
    # if only covariates, directly train a model on them, no RFE
    mypls = plsr(y.train ~ X.train, ncomp = ncol(X.train), validation = "CV")
    comp = selectNcomp(mypls, method="randomization")    
    mycoefs.mean = coef(mypls, ncomps = comp)
    Xmean = mypls$Xmean
    Ymean = mypls$Ymean
    nfeatures = ncol(X.train)
    comps = comp
  } else {
    bagging = do_bagging_RFE_1(X.train, y.train, NITER, alphafeat, maxcomp = maxcomp, selected = selected, n_steps = NSTEPS)
    
    mycoefs.mean = bagging$mycoefs.mean
    Ymean = bagging$Ymean
    Xmean = bagging$Xmean
    nfeatures = bagging$nfeatures
    comps = bagging$comps
  }
  
  
  # get predictions and scores
  y.pred =  t(t(X.test) - Xmean) %*% mycoefs.mean + Ymean 
  y.trainpred =  t(t(X.train) - Xmean) %*% mycoefs.mean + Ymean
  
  RMSE = sqrt(mean((y.pred - y.test)^2 ))
  
  if (var(y.pred)==0) {
    rho = 0
  } else {
    rho = cor.test(y.pred, y.test)$estimate
  }
  
  #plot(y.trainpred, y.train, asp = 1)
  #points(y.pred, y.test, pch=20, col = "red")  
  #abline(0, 1)
  
  # remove covariates, leave only coefficients for RS data
  coefs = mycoefs.mean[!selected]
  
  return(list(fold = fold, y.pred = y.pred, y.test = y.test, rho = rho, RMSE = RMSE, Ymean = Ymean, 
              Xmean = Xmean, coefs = coefs,
              nfeatures = nfeatures, comps = comps))
}


do_crossvalidate_spls = function(fold, input, maxcomp, NITER = 100, steps = seq(0.1, 0.9, 0.2)){
  
  # selecting all features by default  
  X = input$X
  y = input$y
  
  X.train = as.matrix(X[-fold, , drop=FALSE])
  X.test = as.matrix(X[fold, , drop=FALSE])
  
  y.train = y[-fold]  
  y.test = y[fold]
  
  # demean and scale ?
    for (j in seq(ncol(X.train))){
      mu = mean(X.train[, j])
      sigma = sd(X.train[, j])
      X.train[, j] = (X.train[, j] - mu)/sigma     
      X.test[, j] = (X.test[, j] - mu)/sigma    
    }
  
  if (NITER > 1) {
    coefs.iter = y.pred.iter = y.train.pred.iter = NULL    
    for (iter in seq(NITER)){
      print(iter)
      mysample = sample(seq(length(y.train)), replace = TRUE)
      cv <- cv.spls( X.train[mysample, ], y.train[mysample], eta = steps, K = c(1:min(maxcomp, ncol(X.train))), plot.it = F, scale.x = F )
      mypls <- spls( X.train[mysample, ], y.train[mysample], eta = cv$eta.opt, K = cv$K.opt, scale.x = F  )
      #  mypls <- spls( X.train, y.train, eta = 0.3, K = 3  )
      coefs.iter = cbind(coefs.iter, coef.spls(mypls))
      y.pred.iter = cbind(y.pred.iter, predict(mypls, X.test))
      y.train.pred.iter = cbind(y.train.pred.iter, predict(mypls, X.train))
    }
        
    coefs = rowMeans(coefs.iter)
    y.pred = rowMeans(y.pred.iter)
    y.train.pred = rowMeans(y.train.pred.iter)
    
  } else {
    # only one iteration, no resampling
    cv <- cv.spls( X.train, y.train, eta = steps, K = c(1:min(maxcomp, ncol(X.train))), plot.it = F, scale.x = F )
    mypls <- spls( X.train, y.train, eta = cv$eta.opt, K = cv$K.opt, scale.x = F  )
    coefs = coefs.iter = coef.spls(mypls)
    y.pred = predict(mypls, X.test)
    y.train.pred = predict(mypls, X.train)
  }
    
  #print(y.test)
  #print(y.pred)
  #ci = ci.spls(mypls)
  #coefs.correct = correct.spls(ci, plot.it=TRUE)
  #y.pred = as.matrix(ica_data) %*% coefs.correct
  
  RMSE = sqrt(mean((y.pred - y.test)^2 ))
  if (var(y.pred)==0) {
    rho = 0
  } else {
    rho = cor.test(y.pred, y.test)$estimate
  }
  
  return(list(y.pred = y.pred, fold = fold, RMSE=RMSE,
              y.train = y.train, y.train.pred = y.train.pred, 
              y.test = y.test, coefs = coefs, rho = rho,
              eta = cv$eta.opt, K = cv$K.opt, coefs.iter = coefs.iter ))
}

do_crossvalidate_spls_covars = function(fold, input, maxcomp, NITER = 20, covars){
  
  cv = do_crossvalidate_spls(fold, input, maxcomp, 
                             NITER = NITER)
  
  # fit a model with the covariates
  input.covars = list(X = covars, y = input$y) 
  cv.covars = do_crossvalidate_spls(fold, input.covars, maxcomp = ncol(covars), 
                                    NITER = NITER)
  
  # now combine the two models
  y.train = cv$y.train
  y.test = cv$y.test  
  
  X.train = data.frame(nets = cv$y.train.pred, covars = cv.covars$y.train.pred)
  X.test = data.frame(nets = cv$y.pred, covars = cv.covars$y.pred)
  
  
  #mymodel = lm(y.train ~ nets + covars, data = X.train)
  ctrl = trainControl(method = "repeatedcv", number=20, repeats=10, savePredictions = F, allowParallel = T)
  MLmethod = "glmnet" 
  tuneGrid =  expand.grid(.lambda = 10**seq(-2, 2, 0.2), .alpha = c(0, .1, .2, .3, .5, .7) )

  mymodel = train(X.train, y.train,  method = MLmethod, trControl = ctrl, tuneGrid=tuneGrid) 
  y.pred = predict(mymodel, X.test)
  
  RMSE = sqrt(mean((y.pred - y.test)^2 ))
  
  if (var(y.pred)==0) {
    rho = 0
  } else {
    rho = cor.test(y.pred, y.test)$estimate
  }
  
  full = list(
    y.pred = y.pred, 
    y.test = y.test, 
    rho = rho,
    RMSE = RMSE
  )
  return(list(nets = cv, covars = cv.covars, full = full))  
}
