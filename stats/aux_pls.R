
###############################################################################
# generate permutations taking into account group structure
# for permutation testing
###############################################################################
generate_permutation = function(subject, group){
  
  if (is.null(group)){
    mysample = sample(seq(length(subject)), replace = F)
  } else {
    mysample = rep(0, length(subject))
    
    group_size = table(subject)
    unique_subjects = as.numeric(names(group_size))
    
    for (g in sort(unique(group_size))){
#print("------------")
#print(g)
      # shuffle around groups of size g
      mysubjects = unique_subjects[group_size == g]
      if (length(mysubjects) > 1) myshuffle = sample(mysubjects, replace = F)
      else myshuffle = mysubjects
      
#      print(mysubjects)
#      print(myshuffle)
      if (g == 1) {
        indices = sapply(mysubjects, function(x) which(subject == x))
        new.indices = sapply(myshuffle, function(x) which(subject == x))
        mysample[indices] = new.indices
      } else {
        for (j in seq(length(mysubjects))){
          sub = mysubjects[j]
          newsub = myshuffle[j]
          indices = which(subject == sub)
          new.indices = which(subject == newsub)
# cat("indices", indices, "\n")
          mysample[indices] = sample(new.indices, replace = F)
# cat("new", new.indices, "\n")
        }
      }
    }
  }
  return(mysample)
}

###############################################################################
# fit model for each permutation
###############################################################################
doperm = function(perm, y.train, X.train, X.test, maxcomp, subject, group, 
                  filterthr, nonzero, NITER, eta.steps = seq(0.1, 0.9, 0.1)){
  
  # larger eta more sparse
  K.iter = eta.iter = coefs.iter = offset.iter = NULL
  
  set.seed(perm) # to keep iterations consistent across folds
  
  # permutate training data
  if (perm == 1) mysample = seq(length(y.train))
  else mysample = generate_permutation(subject, group)
  
  y.pred.iter = y.train.pred.iter = nfeatures.iter = NULL    
  
  filtered.p = apply(X.train[, nonzero], 2, 
                     function(x) cor.test(y.train[mysample], 
                                          x)$p.value)
  filtered = filtered.p < filterthr
  if(sum(filtered) < maxcomp) {
    filtered = nonzero[order(filtered.p)[1:maxcomp]]
  }
  #  filtered = apply(X.train[, nonzero], 2, 
  #                   function(x) cor.test(y.train[mysample], 
  #                                        x)$p.value < filterthr)
  #  K.steps = seq(1, min(maxcomp, sum(filtered)))
  K.steps = seq(maxcomp)
  
  mycv <- cv.spls( X.train[, nonzero][, filtered], y.train[mysample], 
                   eta = eta.steps, K = K.steps, plot.it = F, scale.x = F )
  
  # iterate for bagging
  if (NITER > 1){
    for (iter in seq(NITER)){
      
      print(paste("Permutation", perm, "Iteration", iter))
      
      #bootstrap
      mysample.b = sample(seq(length(y.train)), replace = TRUE)
      X.train.sample = X.train[mysample.b, nonzero]
      y.train.sample = y.train[mysample][mysample.b]
      
      # feature selection
      filtered.p = apply(X.train.sample, 2, 
                         function(x) cor.test(y.train.sample, 
                                              x)$p.value)
      filtered = filtered.p < filterthr
      if(sum(filtered) < maxcomp) {
        filtered = nonzero[order(filtered.p)[1:maxcomp]]
      }
      
      K.opt = min(mycv$K.opt, sum(filtered))
      mypls <- spls( X.train.sample[, filtered], 
                     y.train.sample, 
                     eta = mycv$eta.opt, 
                     K = K.opt, scale.x = F  )
      
      y.pred.iter = cbind(y.pred.iter, 
                          predict(mypls, X.test[, nonzero][, filtered]))
      y.train.pred.iter = cbind(y.train.pred.iter, 
                                predict(mypls, X.train.sample[, filtered]))
      
      nfeatures.iter[iter] = sum(filtered)
      
      K.iter = c(K.iter, K.opt)
      eta.iter = c(eta.iter, mycv$eta.opt)
      coefs.full = 0*nonzero
      coefs.full[nonzero][filtered] = coef.spls(mypls)
      coefs.iter = cbind(coefs.iter, coefs.full)
      offset.iter = c(offset.iter, as.numeric(mypls$mu))
      
    } # iter loop
    
    #average results over iterations
    y.pred.perm = rowMeans(y.pred.iter)
    y.train.pred.perm = rowMeans(y.train.pred.iter)
    nfeatures.perm = mean(nfeatures.iter)
    
  } else { # only 1 iteration
    
    print(paste("Permutation", perm))
    
    #bootstrap
    X.train.sample = X.train[, nonzero]
    y.train.sample = y.train[mysample]
    
    # feature selection
    filtered.p = apply(X.train.sample, 2, 
                       function(x) cor.test(y.train.sample, 
                                            x)$p.value)
    filtered = filtered.p < filterthr
    if(sum(filtered) < maxcomp) {
      filtered = nonzero[order(filtered.p)[1:maxcomp]]
    }
    
    K.opt = min(mycv$K.opt, sum(filtered))
    mypls <- spls( X.train.sample[, filtered], 
                   y.train.sample, 
                   eta = mycv$eta.opt, 
                   K = K.opt, scale.x = F  )
    
    K.iter = matrix(K.opt, 1, 1)
    eta.iter = matrix(mycv$eta.opt, 1, 1)
    coefs.full = 0*nonzero
    coefs.full[nonzero][filtered] = coef.spls(mypls)
    coefs.iter = matrix(coefs.full, ncol = 1)
    offset.iter = matrix(as.numeric(mypls$mu), 1, 1)
    
    #average results over iterations
    y.pred.perm = matrix(predict(mypls, 
                                 X.test[, nonzero][, filtered]), ncol = 1)
    y.train.pred.perm = matrix(predict(mypls, 
                                       X.train.sample[, filtered]), ncol = 1)
    nfeatures.perm = sum(filtered)
  } # iter if
  
  return(list(
    y.pred.perm = y.pred.perm, 
    y.train.pred.perm = y.train.pred.perm,
    nfeatures.perm = nfeatures.perm, 
    mysample = mysample, 
    coefs.iter = coefs.iter,
    K.iter = K.iter,
    eta.iter = eta.iter,
    offset.iter = offset.iter)
  )
}

###############################################################################
# cross-validated model fitting and testing with permutation tests
###############################################################################
do_crossvalidate_spls_covars_perm_par = function(fold_index, input, 
                                                 maxcomp, 
                                                 cluster, 
                                                 subject = NULL, 
                                                 group = NULL, 
                                                 NITER = 20, 
                                                 NPERM = 100, 
                                                 filterthr = 0.05,
                                                 savecoefs = ''){
  
  # selecting all features by default  
  X = input$X
  y = input$y
  
  covars = input$covars
  fold = input$folds[[fold_index]]
  
  # to generate permutations
  subject = input$subject[-fold]
  if (!is.null(group)) group = input$group[-fold]
  
  nonzero = apply(X, 2, sd) > 0
  
  X.train = as.matrix(X[-fold, , drop=FALSE])
  X.test = as.matrix(X[fold, , drop=FALSE])
  
  y.train = y.train.orig = y[-fold]  
  y.test = y.test.orig = y[fold]
  
  # deconfound
  myform = as.formula(paste("y.train ~ 1 + ", 
                            paste(colnames(covars), collapse=" + ")))  
  
  if(!is.null(covars)){
    covars.train = as.data.frame(covars[-fold, , drop=FALSE])
    covars.test = as.data.frame(covars[fold, , drop=FALSE])
    mycovarmod = lm(myform, data = covars.train)
    y.train = y.train.orig - predict(mycovarmod, covars.train)
    y.test = y.test.orig - predict(mycovarmod, covars.test)
  }
  
  # demean and scale
  mu = colMeans(X.train) 
  sigma = apply(X.train, 2, sd)
  
  for (j in seq(ncol(X.train))){
    X.train[, j] = (X.train[, j] - mu[j])/sigma[j]     
    X.test[, j] = (X.test[, j] - mu[j])/sigma[j]    
  }
  
  # run permutations
  # learn to predict the unconfounded data (y.train)
  results_list = foreach(perm = seq(NPERM), .packages=c('spls'), 
                         .export=c('doperm', 'generate_permutation')) %dopar% 
    doperm(perm, y.train, X.train, X.test, maxcomp, subject, group, 
           filterthr, nonzero, NITER)
  
  # collect results
  y.pred.perm = sapply(results_list, function(x) x$y.pred.perm)
  y.train.pred.perm = sapply(results_list, function(x) x$y.train.pred.perm)
  nfeatures.perm = sapply(results_list, function(x) x$nfeatures.perm)
  
  mysample.perm = sapply(results_list, function(x) x$mysample)
  K.iter = results_list[[1]]$K.iter
  eta.iter = results_list[[1]]$eta.iter
  y.pred = y.pred.perm[, 1]
  y.train.pred = y.train.pred.perm[, 1] 
  coefs.perm = sapply(results_list, function(x) rowMeans(x$coefs.iter))
  coefs = coefs.perm[, 1]
  offset = sapply(results_list, function(x) mean(x$offset.iter))
  
  preprocessing = list(mycovarmod = mycovarmod, 
                       mu = mu,
                       sigma = sigma)
  
  # save these to disk to spare memory
  if ( savecoefs != '') {
    coefs.perm.file = tempfile(pattern = paste0("coefs-", fold_index, "-"), 
                               tmpdir = savecoefs, fileext = ".rda") 
    save(coefs.perm, file = coefs.perm.file) 
    
    preprocessing.file = tempfile(pattern = paste0("preproc-", fold_index, "-"), 
                                  tmpdir = savecoefs, fileext = ".rda")
    save(preprocessing, file = preprocessing.file) 
  } else {
    preprocessing.file = ''
    coefs.perm.file = ''
  }
  
  return(list(fold = fold,
              y.pred = y.pred,
              mysample.perm = mysample.perm,
              y.pred.perm = y.pred.perm,
              y.train = y.train, 
              y.train.pred = y.train.pred, 
              y.test = y.test, 
              y.test.orig = y.test.orig,
              nfeatures.perm = nfeatures.perm,
              coefs = coefs,
              eta.iter = eta.iter, 
              K.iter = K.iter,
              preprocessing.file = preprocessing.file,
              coefs.perm.file = coefs.perm.file,
              coefs = coefs,
              offset = offset
  )
  )
}


