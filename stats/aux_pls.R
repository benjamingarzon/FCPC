
###############################################
# cross-validated spls
###############################################
do_crossvalidate_spls = function(fold, input, maxcomp, NITER = 100, steps = seq(0.1, 0.9, 0.2)){
  
  # selecting all features by default  
  X = input$X
  y = input$y
  print(fold)
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
  if (length(y.test)>1) {
    if (var(y.pred)==0){
      rho = 0
    } else {
      rho = cor.test(y.pred, y.test)$estimate
    }
  } else { 
    rho = 0
  }
  
  return(list(y.pred = y.pred, fold = fold, RMSE=RMSE,
              y.train = y.train, y.train.pred = y.train.pred, 
              y.test = y.test, coefs = coefs, rho = rho,
              eta = cv$eta.opt, K = cv$K.opt, coefs.iter = coefs.iter ))
}

