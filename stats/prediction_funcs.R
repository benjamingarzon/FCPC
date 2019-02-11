

library(spls)
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
MYPALETTE=brewer.pal(11, "Spectral")[seq(11, 1, -1)]# colorspace::diverge_hsv(15, power=2)

FINELINE = 0.5
GGTEXTSIZE2 = 10
GGTEXTSIZE3 = 10
NDIGITS = 3



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
  
plot_adj = function(adj, partition, lim=NULL, pal=MYPALETTE, nolegend=FALSE,main = NULL){
  
  ord = order(partition)
  partition.ord = partition[ord, ]
  
  adj = adj[ord, ord]
  
  rownames(adj) = colnames(adj) = seq(nrow(adj))
  
  #module_names = sort(unique(names(partition.ord)))
  breaks = which(diff(partition.ord)==1) + .5
  maxval = length(partition.ord) + 0.5
  midpoints = (c(0.5, breaks) + c(breaks, maxval )) / 2
  
  molten_adj = melt(adj)
  plot_data = rename(molten_adj, from = Var1, to = Var2) 
  
  if (nolegend) {
    
    theme_new = theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),   
      legend.text = element_text(size=GGTEXTSIZE2),
      legend.title = element_blank(), 
      aspect.ratio = 1,
      legend.position = "none"
    )
    
  } else {
    
    theme_new = theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),   
      legend.text = element_text(size=GGTEXTSIZE2 - 5),
      legend.title = element_blank(), 
      aspect.ratio = 1)
  }
  
  
  if (is.null(lim)) {
    #lim = max(abs(adj), na.rm=TRUE)
    lim = 1
    lim1 = quantile(abs(adj[adj!=0]), .95, na.rm=F)
    lim0 = -lim1 
  }
  
  if (lim == -1) {
    myplot <- ggplot(plot_data, aes(x = from, y = to, fill = factor(value))) +
      geom_raster() +
      theme_bw() +
      theme_new + 
      scale_fill_manual(values=pal) #+
    #scale_y_reverse()
  }
  
  if (lim == -2) {
    lim = 1
    lim0 = quantile(abs(adj), .05, na.rm=TRUE)
    lim1 = quantile(abs(adj), .95, na.rm=TRUE)
  }  
  
  if (lim > 0) {
    myplot <- ggplot(plot_data, aes(x = from, y = to, fill = value)) +
      geom_raster() +
      theme_bw() +
      theme_new + 
      scale_fill_gradientn(colours = pal, 
                           limits=c(lim0, lim1), 
                           oob = squish) #+, na.value = "black"
    #  scale_y_reverse()
  }
  
  
  # separate in modules
  
  for (i in seq(length(breaks))){
    myplot <- myplot + 
      geom_segment(x = 0.5, xend = maxval, y = breaks[i], yend = breaks[i],  size = FINELINE) + 
      geom_segment(y = 0.5, yend = maxval, x = breaks[i], xend = breaks[i],  size = FINELINE)
    
  }
  
  myplot <- myplot + 
    geom_segment(x = 0.5, xend = maxval, y = 0.5, yend = 0.5,  size = FINELINE) + 
    geom_segment(x = 0.5, xend = maxval, y = maxval, yend = maxval,  size = FINELINE) + 
    geom_segment(y = 0.5, yend = maxval, x = 0.5, xend = 0.5,  size = FINELINE) + 
    geom_segment(y = 0.5, yend = maxval, x = maxval, xend = maxval,  size = FINELINE) 
  
  if (nolegend) MYTEXTSIZE = GGTEXTSIZE3
  else MYTEXTSIZE = GGTEXTSIZE3 - 3 
  
  
  for (i in seq(length(midpoints))){
    myplot <- myplot + 
      annotate("text", x = -.5, y = midpoints[i], label = module_names[i], size = MYTEXTSIZE, hjust=1) +
      annotate("text", y = -.5, x = midpoints[i], label = module_names[i], size = MYTEXTSIZE, hjust=1, angle=90) 
  }
  
  myplot = myplot + ylim(-10, maxval) + xlim(-10, maxval) 
  if (!is.null(main)) myplot = myplot + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5, size = GGTEXTSIZE2 + 5))
  
  print(myplot)
}


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
    # adj[modules == to_select[1], modules == to_select[2]] = T
    # adj[modules == to_select[2], modules == to_select[1]] = T
    # adj[modules == to_select[1], modules == to_select[1]] = T
    # adj[modules == to_select[2], modules == to_select[2]] = T
    
    adj[modules == to_select[1], ] = T
    adj[, modules == to_select[1]] = T
    adj[modules == to_select[2], ] = T
    adj[, modules == to_select[2]] = T
    
  }  
  
  valid = which(squareform2(adj))
  print(paste("Connections between module(s) ", to_select[1], ", " , to_select[2], ": " , length(valid), "Inverted: ", invert))
  return(valid)
}


# run prediction for different partitions

do_predictions_loop = function(DIR_ICA_ME, # working dir
                               demo.pars.NSPN, # predictors
                               var, # independent variable
                               mycovars, # covariates
                               modules_file,
                               to_select_list,
                               maxcomp = 10, # max components for spls
                               N_FOLDS = 10, # number of xvalidation folds
                               NPROCS = 10, # processors to use
                               NITER = 50, # spls bagging iteration
                               invert = F)
{
  
  results = list()
  for (i in seq(length(to_select_list)))
  {
    to_select = to_select_list[[i]]
    print(i/length(to_select_list)*100)
    print(to_select)
    results[[i]] = do_prediction(DIR_ICA_ME, 
                                 demo.pars.NSPN, 
                                 var,
                                 mycovars,
                                 maxcomp,
                                 N_FOLDS,
                                 NPROCS,
                                 NITER,
                                 modules_file,
                                 to_select, 
                                 invert)
  }
  
  #  print(mean(results[[i]]$rho))
  return(results)
}



output_results = function(rDIR_ICA_ME, 
                          results,                       
                          RESULTS_DIR = file.path(DIR_ICA_ME, 'mediation/'), modules){
  var = results$var
  mycoefs.mat = squareform(as.numeric(results$coefs.correct.soft))
  
  plot_adj(mycoefs.mat, modules)
    
  
  write.table(mycoefs.mat, file = file.path(RESULTS_DIR, 
                                            paste0(var, "_coefs_mat.txt")
  ), 
  col.names = F, 
  row.names = F)
  
  #  coefs.degree = colSums(mycoefs.mat)
  coefs.degree.abs = colSums(abs(mycoefs.mat))
  coefs.degree.pos = colSums(mycoefs.mat * (mycoefs.mat > 0))  
  coefs.degree.neg = -colSums(mycoefs.mat * (mycoefs.mat < 0))
  coefs.degree.count = colSums(mycoefs.mat != 0)
  
  # plot degrees as maps
  melodic_IC = readNIfTI(file.path(DIR_ICA_ME, 'groupmelodic.ica', 'melodic_IC.nii.gz'))
  bad_nets_file = Sys.glob(file.path(DIR_ICA_ME, paste0("bad_nets", "*.txt")))
  bad_nets = unlist(read.table(bad_nets_file)) + 1
  
  netnames = seq(dim(melodic_IC)[4])[-bad_nets]
  colnames(mycoefs.mat) = rownames(mycoefs.mat) = netnames
  
  # apply threshold
  THR = 12
  melodic_IC = 1 * (melodic_IC[, , , -bad_nets] > THR)
  IC_counts = apply(melodic_IC, 4, sum)
  print(min(IC_counts))
  
  IC.pos = 0*melodic_IC[, , , 1]
  IC.neg = 0*melodic_IC[, , , 1]
  IC.abs = 0*melodic_IC[, , , 1]
  
  n_nets = dim(melodic_IC)[4]
  
  for (i in seq(n_nets)){
    IC.pos = IC.pos + melodic_IC[, , , i] * coefs.degree.pos[i]/n_nets 
    IC.neg = IC.neg + melodic_IC[, , , i] * coefs.degree.neg[i]/n_nets
    IC.abs = IC.abs + melodic_IC[, , , i] * coefs.degree.abs[i]/n_nets
  }
  
  file.pos = file.path(RESULTS_DIR, paste0(var,"_regions.pos"))
  file.neg = file.path(RESULTS_DIR, paste0(var,"_regions.neg"))
  file.abs = file.path(RESULTS_DIR, paste0(var,"_regions.abs"))
  
  writeNIfTI(IC.pos, filename = file.pos)
  writeNIfTI(IC.neg, filename = file.neg)
  writeNIfTI(IC.abs, filename = file.abs)
  
  IC.pos.h = quantile(IC.pos[IC.pos!=0], .95)
  IC.pos.l = quantile(IC.pos[IC.pos!=0], .05)
  IC.neg.h = quantile(IC.neg[IC.neg!=0], .95)
  IC.neg.l = quantile(IC.neg[IC.neg!=0], .05)
  IC.abs.h = quantile(IC.abs[IC.abs!=0], .95)
  IC.abs.l = quantile(IC.abs[IC.abs!=0], .05)
  
  setwd(RESULTS_DIR)  
  system(paste('fslcpgeom', file.path(DIR_ICA_ME, 'bg_image.nii.gz'), file.pos))
  file.copy(paste0(file.pos, '.nii.gz'), 'aux.nii.gz', overwrite = T )
  system('flirt -in aux.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out aux.nii.gz -applyisoxfm 1')
  system('fslmaths aux.nii.gz -s 5 aux.nii.gz')
  #  system(paste('overlay 1 1 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -a aux.nii.gz',
  #               IC.pos.l, IC.pos.h, 'aux.nii.gz'))
  #  system('slicesdir -S aux.nii.gz')
  #  file.copy('slicesdir/aux.png', paste0(file.neg, '.png') )
  
  system(paste('~/Software/FCPC/stats/slices_summary_2 aux.nii.gz', IC.pos.h/10, '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz slicesdir 1 -all'))
  file.copy('slicesdir/0000.png', paste0(file.pos, '.png'), overwrite = T )
  system('rm -r aux.nii.gz slicesdir')  
  
  system(paste('fslcpgeom', file.path(DIR_ICA_ME, 'bg_image.nii.gz'), file.neg))
  file.copy(paste0(file.neg, '.nii.gz'), 'aux.nii.gz', overwrite = T )
  system('flirt -in aux.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out aux.nii.gz -applyisoxfm 1')
  system('fslmaths aux.nii.gz -s 5 aux.nii.gz')
  #  system(paste('overlay 1 1 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -a aux.nii.gz',
  #               IC.neg.l, IC.neg.h , 'aux.nii.gz'))
  #    system('slicesdir aux.nii.gz')
  system(paste('~/Software/FCPC/stats/slices_summary_2 aux.nii.gz', IC.neg.h/10, '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz slicesdir 1 -all'))
  file.copy('slicesdir/0000.png', paste0(file.neg, '.png') , overwrite = T )
  system('rm -r aux.nii.gz slicesdir')  
  
  #  system(paste('fslcpgeom', file.path(DIR_ICA_ME, 'bg_image.nii.gz'), file.abs))
  # file.copy(paste0(file.abs, '.nii.gz'), 'aux.nii.gz')
  # system('flirt -in aux.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out aux.nii.gz -applyisoxfm 1')
  # system(paste('overlay 1 1 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -a aux.nii.gz',
  #              IC.abs.l, IC.abs.h , 'aux.nii.gz'))
  # system('fslmaths aux.nii.gz -s 3 aux.nii.gz')
  # system('slicesdir aux.nii.gz')
  # file.copy('slicesdir/aux.png', paste0(file.abs, '.png') )
  # system('rm -r aux.nii.gz slicesdir')  
}




do_prediction = function(DIR_ICA_ME, # working dir
                         demo.pars.NSPN, # predictors
                         var, # independent variable
                         mycovars, # covariates
                         maxcomp = 10, # max components for spls
                         N_FOLDS = 10, # number of xvalidation folds
                         NPROCS = 10, # processors to use
                         NITER = 50, # spls bagging iteration
                         modules_file = NULL,
                         to_select = NULL, 
                         invert =  F
)
{ 
  print(var)
  results = list()
  results$var = var
  results$empty = F
  # needed files
  FILE_ICA_ME = file.path(DIR_ICA_ME, 'netmats_ridge.csv') #  FILE_ICA_ME = file.path(DIR_ICA_ME, 'metrics_ridge.csv')
  FILE_BRAIN_VOL = file.path(DIR_ICA_ME, 'brain_volumes.txt')
  FILE_DOF = file.path(DIR_ICA_ME, 'meica_DOF_nom.txt')
  FILE_FD = file.path(DIR_ICA_ME, 'fd.txt')
  FILE_DVARS = file.path(DIR_ICA_ME, 'dvars.txt')
  FILE_SUBJECTS = file.path(DIR_ICA_ME, 'subjects.txt')
  
  # load and organize covariates
  brain_vol = read.table(FILE_BRAIN_VOL)
  meica_dof = read.table(FILE_DOF)
  fd = read.table(FILE_FD)
  fd = cbind(fd[1], rowMeans(fd[-c(1,2)]))
  dvars = read.table(FILE_DVARS)
  dvars = cbind(fd[1], rowMeans(dvars[-c(1,2)]))
  subjects = as.numeric(read.table(FILE_SUBJECTS))
  ica_data = read.table(FILE_ICA_ME, sep = ',', head=F)
  colnames(ica_data) = paste0("feature.", seq(ncol(ica_data)))
  ica_data$Subject = subjects
  
  # sync data and covariates
  mydata = merge(ica_data, demo.pars.NSPN, by="Subject")
  features = grep('feature', colnames(mydata))
  
  covars = merge(dvars, fd, by = "V1")
  covars = merge(covars, brain_vol, by = "V1")
  covars = merge(covars, meica_dof, by = "V1")
  colnames(covars) = c("Subject","dvars", "fd","brain_vol","meica_dof")
  mydata = merge(mydata, covars, by ="Subject")
  # print(file.path(DIR_ICA_ME, 'full_data.csv'))
  # write.table(mydata, file =  file.path(DIR_ICA_ME, 'full_data.csv'))  
  
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
  
  mydata.exp = cbind(
    mydata$Subject,
    mydata[features],
    mydata[c(var, unique(c(mycovars, "age", "sex.num")))]   
  ) 
  
  # add higher order terms and interactions
  mydata.exp = mydata.exp[complete.cases(mydata.exp), ]
  mydata.exp$age.1 = as.numeric(scale(mydata.exp$age))
  mydata.exp$age.2 = mydata.exp$age.1^2
  mydata.exp$age.1.sex = mydata.exp$age.1*mydata.exp$sex.num
  mydata.exp$age.2.sex = mydata.exp$age.2*mydata.exp$sex.num
  
  y = unlist(mydata.exp[var])
  
  #remove outliers if any
  valid.rows = (y - mean(y)) < 3*sd(y)
  
  mycovars = c(mycovars, "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex")
  covars.exp = mydata.exp[valid.rows, mycovars]
  
  features = grep('feature', colnames(mydata.exp))
  ica_data = as.matrix(mydata.exp[valid.rows, features])
  y = unlist(mydata.exp[var])[valid.rows]
  
  print(dim(ica_data))
  
  if(!is.null(N_FOLDS)){
    cl <- makeCluster(NPROCS)
    
    registerDoParallel(cl)    
    folds = createFolds(y, k = N_FOLDS)
    cv = foreach(fold = folds, .packages=c('spls'), 
                 .export=c('do_crossvalidate_spls')) %dopar% 
      do_crossvalidate_spls(fold, list(X = ica_data, y = y), 
                            maxcomp = maxcomp, NITER = NITER)
    stopCluster(cl)
    
    results$rho = unlist(lapply(cv, function(x){x$rho}))
    results$RMSE = unlist(lapply(cv, function(x){x$RMSE}))
    
    results$y.pred = unlist(lapply(cv, function(x){x$y.pred}))
    results$y.test = unlist(lapply(cv, function(x){x$y.test}))
    results$fold = unlist(lapply(cv, function(x){x$fold}))      
    results$eta = unlist(lapply(cv, function(x){x$eta}))      
    results$K = unlist(lapply(cv, function(x){x$K}))  
    results$coefs.cv = do.call("cbind", lapply(cv, function(x){x$coefs.iter}) )
    
    results$coefs.correct.soft = apply(results$coefs.cv, 1, 
                                       function(x) {
                                         q = quantile(x, c(.005, .995))
                                         ifelse( prod(q) < 0, 0, mean(x))
                                       })
    
    results$coefs.correct = apply(results$coefs.cv, 1,
                                  function(x) ifelse( prod(range(x)) < 0, 0, mean(x))
    )
    
    # average y.pred across folds
    results$y.pred.mean = unlist(Hmisc::summarize(results$y.pred, by = as.factor(results$fold), mean)[2])
    results$y.test.mean = unlist(Hmisc::summarize(results$y.test, by = as.factor(results$fold), mean)[2])
  } else {
    cv <- cv.spls( ica_data, y, eta = seq(0.1, 0.9, 0.1), K = c(1:maxcomp), plot.it = F  )
    mypls <- spls( ica_data, y, eta = cv$eta.opt, K = cv$K.opt)
    #   mypls <- spls( ica_data, y, eta = 0.3, K = 3  )
    
    coefs = coef.spls(mypls)
    ci = ci.spls(mypls)
    results$coefs.correct = correct.spls(ci, plot.it=TRUE)
    #y.pred = as.matrix(ica_data) %*% coefs.correct
    y.pred = predict(mypls, ica_data)
    results$y.pred.mean = results$y.pred = y.pred
    results$y.test.mean = results$y.test = y
    results$fold = seq(length(y))  
    
      }
  
  
  results$covars.exp = covars.exp
  results$nfeatures = ncol(ica_data)
  results$cor.test = cor.test(results$y.pred.mean, results$y.test.mean, alternative = "greater")
  results$cor.test.twotailed = cor.test(results$y.pred.mean, results$y.test.mean)
  results$Subject = mydata.exp$`mydata$Subject`[valid.rows]
  
  print("Coefs: ")
  print(results$cor.test)
  print(dim(results$coefs.cv))
  print(sum(results$coefs.correct != 0))
  print(sum(results$coefs.soft != 0))
  
  return(results)
}



do_mediation = function(DIR_ICA_ME, 
                        results, 
                        NSIMS = 2000,
                        write_coefs = F,
                        RESULTS_DIR = file.path(DIR_ICA_ME, 'mediation/')
)
{ 
  FIGS_DIR <<- RESULTS_DIR
  
  print(RESULTS_DIR)
  
  #  dir.create(RESULTS_DIR)
  
  data.M = results$covars.exp
  data.M$y = results$y.test.mean
  data.M$M = results$y.pred.mean
  data.M$sex = as.factor(data.M$sex.num)
  
  print(paste("Analyzing", nrow(results$covars.exp), "observations."))
  
  # 1) X -> Y
  model.0 = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU, data = data.M)  
  print(summary(model.0))
  
  # 2) X -> M
  model.M = lm(M ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU, data = data.M)
  print(summary(model.M))
  
  # 3) X + M -> Y
  model.Y = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU + M, data = data.M)  
  print(summary(model.Y))
  
  # mediation ?
  mediation.age = mediate(model.M, model.Y, treat = "age.1", mediator = "M", boot = T, sims = NSIMS)
  mediation.sex = mediate(model.M, model.Y, treat = "sex", mediator = "M", boot = T, sims = NSIMS)
  
  print(summary(mediation.age))
  print(summary(mediation.sex))
  
  results$mediation.age.p = mediation.age$d0.p
  results$mediation.sex.p = mediation.sex$d0.p
  
  #output the model coefficients as well
  if (write_coefs){
    mycoefs.mat = squareform(as.numeric(results$coefs.correct))
    write.table(mycoefs.mat, file = file.path(RESULTS_DIR, 
                                              paste0(var, "_coefs_mat.txt")
    ), 
    col.names = F, 
    row.names = F)
    
  }
  
  results$FIGS_DIR = FIGS_DIR
  results$RESULTS_DIR = RESULTS_DIR
  results$DIR_ICA_ME = DIR_ICA_ME
  results$model.0.results = round(summary(model.0)$coefficients[c('age.1', 'sex1'), c(1, 3, 4)], digits = NDIGITS)
  results$model.M.results = round(summary(model.M)$coefficients[c('age.1', 'sex1'), c(1, 3, 4)], digits = NDIGITS)
  results$model.Y.results = round(summary(model.Y)$coefficients[c('age.1', 'sex1', 'M'), c(1, 3, 4)], digits = NDIGITS)
  colnames(results$model.0.results) = colnames(results$model.M.results) = colnames(results$model.Y.results) = c("Estimate", "t", "p")
  
  mediation.age.sum = summary(mediation.age) # ACME d.avg, d0.ci, d0.p /ADE z... / TOTAL tau /prop n.avg
  mediation.sex.sum = summary(mediation.sex) # ACME d.avg, d0.ci, d0.p /ADE z... / TOTAL tau /prop n.avg
  
  results$mediation.age = rbind(unlist(mediation.age.sum[ c('d.avg', 'd.avg.ci', 'd.avg.p')]),
                                unlist(mediation.age.sum[ c('z.avg', 'z.avg.ci', 'z.avg.p')]),
                                unlist(mediation.age.sum[ c('tau.coef', 'tau.ci', 'tau.p')]),
                                unlist(mediation.age.sum[ c('n.avg', 'n.avg.ci', 'n.avg.p')]))
  
  results$mediation.sex = rbind(unlist(mediation.sex.sum[ c('d.avg', 'd.avg.ci', 'd.avg.p')]),
                                unlist(mediation.sex.sum[ c('z.avg', 'z.avg.ci', 'z.avg.p')]),
                                unlist(mediation.sex.sum[ c('tau.coef', 'tau.ci', 'tau.p')]),
                                unlist(mediation.sex.sum[ c('n.avg', 'n.avg.ci', 'n.avg.p')]))
  
  colnames(results$mediation.age) = colnames(results$mediation.sex) = c("Average", "Cil", "Cih", "p")
  rownames(results$mediation.age) = rownames(results$mediation.sex) = c("ACME", "ADE", "Total", "Proportion")
  
  return(results)
}


################################################################################################################## 

# old functions... 

do_mediation_old = function(DIR_ICA_ME, # working dir
                            demo.pars.NSPN, # predictors
                            var, # independent variable
                            mycovars, # covariates
                            NSIMS = 2000, # for mediation analysis
                            maxcomp = 10, # max components for spls
                            N_FOLDS = 10, # number of xvalidation folds
                            NPROCS = 10, # processors to use
                            NITER = 50 # spls bagging iteration
)
{ 
  print(var)
  results = list()
  results$var = var
  
  # create figs and results dir
  FIGS_DIR <<- file.path(DIR_ICA_ME, 'figs/')
  RESULTS_DIR = file.path(DIR_ICA_ME, 'mediation/')
  
  dir.create(RESULTS_DIR)
  
  results = do_prediction(DIR_ICA_ME, # working dir
                          demo.pars.NSPN, # predictors
                          var, # independent variable
                          mycovars, # covariates
                          maxcomp, # max components for spls
                          N_FOLDS, # number of xvalidation folds
                          NPROCS, # processors to use
                          NITER # spls bagging iteration
  )
  
  data.M = results$covars.exp
  data.M$y = results$y.test.mean
  data.M$M = results$y.pred.mean
  data.M$sex = as.factor(data.M$sex.num)
  
  print(paste("Analyzing", nrow(results$covars.exp), "observations."))
  
  # 1) X -> Y
  model.0 = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU, data = data.M)  
  print(summary(model.0))
  
  # 2) X -> M
  model.M = lm(M ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU, data = data.M)
  print(summary(model.M))
  
  # 3) X + M -> Y
  model.Y = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU + M, data = data.M)  
  print(summary(model.Y))
  
  mediation.age = mediate(model.M, model.Y, treat = "age.1", mediator = "M", boot = T, sims = NSIMS)
  mediation.sex = mediate(model.M, model.Y, treat = "sex", mediator = "M", boot = T, sims = NSIMS)
  
  print(summary(mediation.age))
  print(summary(mediation.sex))
  
  
  #mycoefs.mat = squareform(as.numeric(results$coefs.correct))
  results$mediation.age.p = mediation.age$d0.p
  results$mediation.sex.p = mediation.sex$d0.p
  
  results$FIGS_DIR = FIGS_DIR
  results$RESULTS_DIR = RESULTS_DIR
  results$DIR_ICA_ME = DIR_ICA_ME
  
  return(results)
}



do_spls_hierarchical = function(DIR_ICA_ME, # working dir
                                demo.pars.NSPN, # predictors
                                var, # independent variable
                                mycovars, # covariates
                                maxcomp = 10, # maximum number of components in SPLS 
                                N_FOLDS = 10, # number of folds
                                NPROCS = 10, # processors
                                TEST_PROP = 0.1, # 
                                NITER = 50 # bagging iterations
)
{ 
  print(var)
  results = list()
  results$var = var
  
  # create figs and results dir
  FIGS_DIR <<- file.path(DIR_ICA_ME, 'figs/')
  RESULTS_DIR = file.path(DIR_ICA_ME, 'results/')
  
  # load partition of ROIs
  dir.create(RESULTS_DIR)
  modules_file = file.path(DIR_ICA_ME, 'BrainMapPartition.txt')
  modules = read.table(modules_file)
  
  # needed files
  FILE_ICA_ME = file.path(DIR_ICA_ME, 'netmats_ridge.csv')
  #  FILE_ICA_ME = file.path(DIR_ICA_ME, 'netmats_full.csv')
  #  FILE_ICA_ME = file.path(DIR_ICA_ME, 'metrics_ridge.csv')
  FILE_BRAIN_VOL = file.path(DIR_ICA_ME, 'brain_volumes.txt')
  FILE_DOF = file.path(DIR_ICA_ME, 'meica_DOF_nom.txt')
  FILE_FD = file.path(DIR_ICA_ME, 'fd.txt')
  FILE_DVARS = file.path(DIR_ICA_ME, 'dvars.txt')
  FILE_SUBJECTS = file.path(DIR_ICA_ME, 'subjects.txt')
  
  # load and organize covariates
  brain_vol = read.table(FILE_BRAIN_VOL)
  meica_dof = read.table(FILE_DOF)
  fd = read.table(FILE_FD)
  fd = cbind(fd[1], rowMeans(fd[-c(1,2)]))
  dvars = read.table(FILE_DVARS)
  dvars = cbind(fd[1], rowMeans(dvars[-c(1,2)]))
  subjects = as.numeric(read.table(FILE_SUBJECTS))
  ica_data = read.table(FILE_ICA_ME, sep = ',', head=F)
  colnames(ica_data) = paste0("feature.", seq(ncol(ica_data)))
  ica_data$Subject = subjects
  
  # sync data and covariates
  mydata = merge(ica_data, demo.pars.NSPN, by="Subject")
  features = grep('feature', colnames(mydata))
  
  covars = merge(dvars, fd, by = "V1")
  covars = merge(covars, brain_vol, by = "V1")
  covars = merge(covars, meica_dof, by = "V1")
  colnames(covars) = c("Subject","dvars", "fd","brain_vol","meica_dof")
  mydata = merge(mydata, covars, by ="Subject")
  
  #  write.table(mydata, file =  file.path(DIR_ICA_ME, 'full_data.csv'))  
  features = grep('feature', colnames(mydata))
  
  mydata.exp = cbind(
    mydata$Subject,
    mydata[features],
    mydata[c(var, unique(c(mycovars, "age", "sex.num")))]   
  ) 
  
  # add higher order terms and interactions
  mydata.exp = mydata.exp[complete.cases(mydata.exp), ]
  mydata.exp$age.1 = as.numeric(scale(mydata.exp$age))
  mydata.exp$age.2 = mydata.exp$age.1^2
  mydata.exp$age.1.sex = mydata.exp$age.1*mydata.exp$sex.num
  mydata.exp$age.2.sex = mydata.exp$age.2*mydata.exp$sex.num
  
  y = unlist(mydata.exp[var])
  #outliers ?
  valid = (y - mean(y)) < 3*sd(y)
  
  mycovars = c(mycovars, "age.1", "age.2", "sex.num", "age.1.sex", "age.2.sex")
  covars.exp = mydata.exp[valid, mycovars]
  
  features = grep('feature', colnames(mydata.exp))
  ica_data = as.matrix(mydata.exp[valid, features])
  y = unlist(mydata.exp[var])[valid]
  
  cl <- makeCluster(NPROCS)
  
  registerDoParallel(cl)    
  folds = createFolds(y, k = N_FOLDS)
  
  cv = foreach(fold = folds, .packages=c('spls', 'caret', 'glmnet'), 
               .export=c('do_crossvalidate_spls_covars', 'do_crossvalidate_spls')) %dopar% 
    do_crossvalidate_spls_covars(fold, list(X = ica_data, y = y), 
                                 maxcomp = maxcomp, NITER = NITER, covars.exp)
  
  
  stopCluster(cl)
  results$rho.full = unlist(lapply(cv, function(x){x$full$rho}))
  results$rho.nets = unlist(lapply(cv, function(x){x$nets$rho}))
  results$rho.covars = unlist(lapply(cv, function(x){x$covars$rho}))
  results$RMSE.full = unlist(lapply(cv, function(x){x$full$RMSE}))
  results$RMSE.nets = unlist(lapply(cv, function(x){x$nets$RMSE}))
  results$RMSE.covars = unlist(lapply(cv, function(x){x$covars$RMSE}))
  
  y.pred.full = unlist(lapply(cv, function(x){x$full$y.pred}))
  y.pred.nets = unlist(lapply(cv, function(x){x$nets$y.pred}))
  y.pred.covars = unlist(lapply(cv, function(x){x$covars$y.pred}))
  
  y.test = unlist(lapply(cv, function(x){x$full$y.test}))
  
  results$full.cor = cor.test(y.pred.full, y.test)
  results$nets.cor = cor.test(y.pred.nets, y.test)
  results$covars.cor = cor.test(y.pred.covars, y.test)
  
  results$FIGS_DIR = FIGS_DIR
  results$RESULTS_DIR = RESULTS_DIR
  results$DIR_ICA_ME = DIR_ICA_ME
  
  results$RMSE.diff = t.test(results$RMSE.full - results$RMSE.covars)
  results$cv = cv
  return(results)
}


plot_results_hierarchical = function(results){ 
  var = results$var
  
  cv = results$cv
  
  eta = unlist(lapply(cv, function(x){x$nets$eta}))      
  K = unlist(lapply(cv, function(x){x$nets$K}))  
  coefs.cv = do.call("cbind", lapply(cv, function(x){x$nets$coefs}) )
  coefs.iter = do.call("cbind", lapply(cv, function(x){x$nets$coefs.iter}) )
  
  # t test
  alpha = 0.1
  coefs.correct = apply(coefs.iter, 1, 
                        function(x) {
                          q = quantile(x, c(alpha, 1-alpha))
                          ifelse( (q[1] > 0)|(q[2] < 0), mean(x), 0)
                        }
  )
  
  coefs.t = apply(coefs.iter, 1, function(x) t.test(x)$estimate)
  coefs.mean = rowMeans(coefs.iter)
  coefs.sd = apply(coefs.iter, 1, function(x) sd(x))
  #  plot(coefs.t, coefs.correct)
  #  plot(abs(coefs.mean), coefs.sd)
  #  points(abs(coefs.correct), coefs.sd, col = 'red', pch = 20)
  #  matplot(coefs.iter[coefs.correct!=0, ], pch = 20, col = 1)
  #  hist(coefs.mean, 100)
  mycoefs.mat = squareform(as.numeric(coefs.correct))
  
  #modules = read.table(file.path(results$DIR_ICA_ME, 'modules.txt'))
  modules = read.table(file.path(results$DIR_ICA_ME, 'BrainMapPartition.txt'))
  
  FIGS_DIR <<- results$FIGS_DIR
  dir.create(FIGS_DIR)
  setwd(results$RESULTS_DIR)
  bad_nets_file = Sys.glob(file.path(results$DIR_ICA_ME, paste0("bad_nets", "*.txt")))
  bad_nets = unlist(read.table(bad_nets_file)) + 1
  
  #--------------------------------
  # represent by modules
  #--------------------------------
  n_modules = max(modules)
  print(n_modules)
  coefs.modules.pos = coefs.modules.neg = coefs.modules = coefs.modules.mean = size = matrix(0, n_modules, n_modules)
  
  for (i in seq(n_modules)){
    for (j in seq(n_modules)){
      coefs.mat = mycoefs.mat[modules == i, modules == j]
      
      if (i == j){
        
        coefs.modules[i, j] = sum(coefs.mat != 0)/2
        coefs.modules.mean[i, j] = mean(abs(coefs.mat[coefs.mat != 0]))        
        coefs.modules.pos[i, j] = sum(coefs.mat > 0)/2
        coefs.modules.neg[i, j] = sum(coefs.mat < 0)/2
        size[i, j] = (length(coefs.mat) - sum(modules == i))/2
        
      } else {
        
        coefs.modules[i, j] = sum(coefs.mat != 0)
        coefs.modules.mean[i, j] = mean(abs(coefs.mat[coefs.mat != 0]))        
        coefs.modules.pos[i, j] = sum(coefs.mat > 0)
        coefs.modules.neg[i, j] = sum(coefs.mat < 0)
        size[i, j] = length(coefs.mat)
        
      }      
    }
  }     
  
  # fraction of features correcting for size
  frac = (coefs.modules/sum(mycoefs.mat!=0))/(size/length(squareform(mycoefs.mat)))
  coefs.modules.frac = melt(frac)
  coefs.modules.frac$type = 'fraction'
  coefs.modules.count = melt(coefs.modules) #/max(coefs.modules)
  coefs.modules.count$type = 'count'
  coefs.modules.mean = melt(coefs.modules.mean)
  coefs.modules.mean$type = 'mean'
  coefs.modules.size = melt(size/length(squareform(mycoefs.mat))*10)
  coefs.modules.size$type = 'size'
  coefs.modules.pos = melt(coefs.modules.pos)
  coefs.modules.pos$type = 'pos'
  coefs.modules.neg = melt(coefs.modules.neg)
  coefs.modules.neg$type = 'neg'
  
  save_fig(paste0(var, '_count_sum_modules')) 
  plot(colSums(coefs.modules, na.rm=T), type = "b", pch = 20)
  
  save_fig(paste0(var, '_size_sum_modules')) 
  plot(colSums(size, na.rm=T), type = "b", pch = 20)
  
  coefs.modules.count[coefs.modules.count == 0] = NA
  coefs.modules.frac[coefs.modules.frac == 0] = NA
  #coefs.modules = rbind(coefs.modules.count, coefs.modules.size)
  #, coefs.modules.count, coefs.modules.frac, coefs.modules.pos, coefs.modules.neg
  myplot = ggplot(data = coefs.modules.count, aes(x = sprintf("%02d", X1), y = sprintf("%02d", X2), fill=value)) + 
    geom_tile() + xlab("Module") + ylab("Module") + theme(text = element_text(size=20), axis.text.x = element_text(size = 20),
                                                          axis.text.y = element_text(size = 20))#  + facet_grid(.~ type)
  
  save_fig(paste0(var, '_coefsmat_modules')) 
  print(myplot)
  
  
  dev.off()
  #--------------------------------
  # represent by region
  #--------------------------------
  # need some cleaning, perhaps get highest 10% of connections for each net?
  
  coefs.degree = colSums(mycoefs.mat)
  coefs.degree.abs = colSums(abs(mycoefs.mat))
  coefs.degree.pos = colSums(mycoefs.mat * (mycoefs.mat > 0))  
  coefs.degree.neg = colSums(mycoefs.mat * (mycoefs.mat < 0))
  coefs.degree.count = colSums(mycoefs.mat != 0)
  coefs.degree.count.pos = colSums(mycoefs.mat > 0)  
  coefs.degree.count.neg = colSums(mycoefs.mat < 0)  
  coefs.degree.mean = apply(mycoefs.mat, 2, function(x) {mean(abs(x[x!=0])) })
  #coefs.degree.mean = colMeans(mycoefs.mat != 0)
  
  write.table(mycoefs.mat,
              file = paste0(var,"_coefs.mat.txt"), col.names = F, row.names = F)
  
  # plot degrees as maps
  melodic_IC = readNIfTI(file.path(results$DIR_ICA_ME, 'groupmelodic.ica', 'melodic_IC.nii.gz'))
  netnames = seq(dim(melodic_IC)[4])[-bad_nets]
  colnames(mycoefs.mat) = rownames(mycoefs.mat) = netnames
  
  # apply threshold
  THR = 12
  melodic_IC = 1 * (melodic_IC[, , , -bad_nets] > THR)
  IC_counts = apply(melodic_IC, 4, sum)
  print(min(IC_counts))
  
  IC.pos = 0*melodic_IC#[, , , 1]
  IC.neg = 0*melodic_IC#[, , , 1]
  IC.abs = 0*melodic_IC#[, , , 1]
  
  n_nets = dim(melodic_IC)[4]
  
  
  for (i in seq(n_nets)){
    #IC.pos = IC.pos + melodic_IC[, , , i] * coefs.degree.pos[i]#/n_nets 
    #IC.neg = IC.neg + melodic_IC[, , , i] * coefs.degree.neg[i]#/n_nets
    #IC.abs = IC.abs + melodic_IC[, , , i] * coefs.degree.count[i]#coefs.degree.abs[i]#/n_nets
    IC.pos[, , , i] = melodic_IC[, , , i] * coefs.degree.count.pos[i] 
    IC.neg[, , , i] = melodic_IC[, , , i] * coefs.degree.count.neg[i] #/n_nets
    IC.abs[, , , i] = melodic_IC[, , , i] * coefs.degree.count[i] #coefs.degree.abs[i]#/n_nets
    
  }
  IC.pos = apply(IC.pos, c(1,2,3), max)
  IC.neg = apply(IC.neg, c(1,2,3), max)
  IC.abs = apply(IC.abs, c(1,2,3), max)
  
  file.pos = paste0(var,"_regions.pos")
  file.neg = paste0(var,"_regions.neg")
  file.abs = paste0(var,"_regions.abs")
  writeNIfTI(IC.pos, filename = file.pos)
  writeNIfTI(IC.neg, filename = file.neg)
  writeNIfTI(IC.abs, filename = file.abs)
  
  IC.pos.h = quantile(IC.pos[IC.pos!=0], .95)
  IC.pos.l = quantile(IC.pos[IC.pos!=0], .05)
  IC.neg.h = quantile(IC.neg[IC.neg!=0], .95)
  IC.neg.l = quantile(IC.neg[IC.neg!=0], .05)
  IC.abs.h = quantile(IC.abs[IC.abs!=0], .95)
  IC.abs.l = quantile(IC.abs[IC.abs!=0], .05)
  
  
  # need to change the threshold in overlay!!
  # plot as overlay 
  system(paste('fslcpgeom ../bg_image', file.pos))
  file.copy(paste0(file.pos, '.nii.gz'), 'aux.nii.gz')
  system('flirt -in aux.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out aux.nii.gz -applyisoxfm 1')
  system(paste('overlay 1 1 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -a aux.nii.gz',
               IC.pos.l, IC.pos.h, 'aux.nii.gz'))
  system('slicesdir -S aux.nii.gz')
  file.copy('slicesdir/aux.png', file.path(FIGS_DIR, paste0(file.pos, '.png') ))
  system('rm -r aux.nii.gz slicesdir')  
  
  system(paste('fslcpgeom ../bg_image', file.neg))
  file.copy(paste0(file.neg, '.nii.gz'), 'aux.nii.gz')
  system('flirt -in aux.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out aux.nii.gz -applyisoxfm 1')
  system(paste('overlay 1 1 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -a aux.nii.gz',
               IC.neg.l, IC.neg.h , 'aux.nii.gz'))
  system('slicesdir -S aux.nii.gz')
  file.copy('slicesdir/aux.png', file.path(FIGS_DIR, paste0(file.neg, '.png') ))
  system('rm -r aux.nii.gz slicesdir')  
  
  system(paste('fslcpgeom ../bg_image', file.abs))
  file.copy(paste0(file.abs, '.nii.gz'), 'aux.nii.gz')
  system('flirt -in aux.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out aux.nii.gz -applyisoxfm 1')
  system(paste('overlay 1 1 /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -a aux.nii.gz',
               IC.abs.l, IC.abs.h , 'aux.nii.gz'))
  system('slicesdir -S aux.nii.gz')
  file.copy('slicesdir/aux.png', file.path(FIGS_DIR, paste0(file.abs, '.png') ))
  system('rm -r aux.nii.gz slicesdir')  
  
  
  # plot the prediction
  y.pred = unlist(lapply(cv, function(x){x$nets$y.pred}))
  y.test = unlist(lapply(cv, function(x){x$nets$y.test}))
  fold = unlist(lapply(cv, function(x){x$nets$fold}))
  
  model = lm(y.pred ~ 1 + y.test)
  print(summary(model))
  #print(cor.test(y.test, y.pred)) 
  plot(y.test, y.pred)
  abline(model, col = 'red')
  abline(0, 1)
  
  save_fig(paste0(var, '_prediction')) 
  print(myplot)
  dev.off()
  
  
  
}

addstar = function(x, fdr_correct) {
  
  fdrp = p.adjust(x[fdr_correct], method = "fdr")
  sig = ifelse(x < alpha, "* ", "  ")
  sig[fdr_correct][fdrp < alpha] = "**"
  return(paste0(x, sig))
}

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
  dt.Y.M = cbind(module_labels, model.Y.results.M) [myorder, ]
  dt.age = cbind(module_labels, model.M.results.age, model.Y.results.age, mediation.age.ACME)[myorder, ]
  dt.sex = cbind(module_labels, model.M.results.sex, model.Y.results.sex, mediation.sex.ACME)[myorder, ]
  
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


plot_combinations = function(myresults, level, module_names, measure){
  
  print("--------------------------")
  alpha = 0.05
  
  rho = data.frame(
    mean = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test$estimate, NA)), 
    low = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.twotailed$conf.int[1], NA)),
    high = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.twotailed$conf.int[2], NA)),
    p = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test$p.value, NA)),
    nfeatures = sapply(myresults, function(x) ifelse(!x$empty, x$nfeatures, NA)),
    from = as.factor(sapply(myresults, function(x) x$to_select[1])),
    to = as.factor(sapply(myresults, function(x) x$to_select[2])),
    mediation.age.p = sapply(myresults, function(x) x$mediation.age.p),
    mediation.sex.p = sapply(myresults, function(x) x$mediation.sex.p)
  )        
  
  rho$mediation.sig = 'none'
  rho$mediation.sig[rho$mediation.age.p < alpha] = 'age'
  rho$mediation.sig[rho$mediation.sex.p < alpha] = 'sex'
  rho$mediation.sig[rho$mediation.age.p < alpha & rho$mediation.sex.p < alpha] = 'age & sex'
  
  
  # check the prediction is not driven by the number of features
  print("Correlation between number of features and mean rho (abs and original)")
  print(cor.test(rho$nfeatures, abs(rho$mean)))  
  print(cor.test(rho$nfeatures, rho$mean, alternative = "greater"))  
  print(
    ggplot(rho, aes( x = nfeatures, y = mean, ymin = low, ymax = high )) + 
      geom_pointrange() + 
      xlab("Number of features") +
      ylab("Correlation between observed and predicted") + 
      theme_classic() 
  )
  
  #rho$mean[rho$p > alpha] = NA
  #rho$mediation.sig[rho$p > alpha] = NA
  #rho = subset(rho, p < alpha)
  rho$width = ifelse(rho$p < alpha, 1, 0.5)
  rho$size = ifelse(rho$p < alpha, 5, 4)
  
  print(
    ggplot(rho, aes(x = from, y = to, fill = mean)) +
      #  ggtitle('Correlation between true and predicted') +
      geom_tile() +
      labs(title = 'Correlation between true and predicted',
           x = 'Network',
           y = 'Network',
           fill = 'Correlation'
      ) + 
      scale_x_discrete(position = "bottom") +
      scale_y_discrete(position = "right") +
      scale_fill_gradient2(midpoint = .1, mid ="orange", high = "yellow", low = "red", limits = c(0, 0.3)) +
      geom_text(aes(x = from, y = to, label = round(mean, 2), size = size), color = 'black', fontface = "bold") +
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
  )
  
  # plot mediations
  print(
    ggplot(rho, aes(x = from, y = to, fill = mediation.sig)) +
      #  ggtitle('Correlation between true and predicted') +
      geom_tile() +
      labs(title = 'Correlation between observed and predicted',
           x = 'Network',
           y = 'Network',
           fill = 'Mediator'
      ) + 
      scale_x_discrete(position = "bottom") +
      scale_y_discrete(position = "right") +
      geom_text(aes(x = from, y = to, label = round(mean, 2), size = size), color = 'black', fontface = "bold") +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.3, 0.7),
        legend.direction = "vertical")) #+
  #      guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5) )
  #  )
  
  save_fig(figname = 'Accuracy', res = BWRES)
  rho.diag = subset(rho, from == to)
  rho.diag$network = as.factor(module_names[rho.diag$from]) 
  
  rho.diag$fdrp  = p.adjust(rho.diag$p, method = "fdr") 
  rho.diag$sig = ifelse(rho.diag$p < alpha, "*",  "")
  rho.diag$sig[rho.diag$fdrp < alpha] = "**"
  
  print(rho.diag)
  print(
    ggplot(rho.diag, aes(network, mean, label = rho.diag$sig)) +
      ggtitle(paste('Correlation between observed and predicted', measure, "\n")) + xlab('Network') + 
      geom_text(size = CEX_TEXT) +
      geom_hline(yintercept = level, lty = 2, size = CEX_LINE) + 
      geom_col(col = 'black', fill = 'gray', alpha = 0.5) + ylim(0, 0.3) + 
      theme_bw() + 
      xlab('\nModule') + 
      ylab('r') +
      
      theme(
        axis.text = element_text(size = CEX_AXIS, margin = margin(t = 100, r = 0, b = 100, l = 0)),
        axis.title = element_text(size = CEX_AXIS, margin = margin(t = 100, r = 0, b = 100, l = 0)),
        plot.title = element_text(hjust = 0.5, size = CEX_TITLE),
        plot.margin = margin(2, 0, 2, 0, "cm"),
        axis.line = element_line(size = 2)
      )
  ) 
  dev.off()
  
  return(rho)
}



do_predictions_loop_deprec = function(DIR_ICA_ME, # working dir
                                      demo.pars.NSPN, # predictors
                                      var, # independent variable
                                      mycovars, # covariates
                                      modules_file,
                                      to_select_list,
                                      maxcomp = 10, # max components for spls
                                      N_FOLDS = 10, # number of xvalidation folds
                                      NPROCS = 10, # processors to use
                                      NITER = 50 # spls bagging iteration
)
{
  
  results = list()
  for (i in seq(length(to_select_list)))
  {
    to_select = to_select_list[i]
    print(paste(i, to_select))
    results[[i]] = do_prediction(DIR_ICA_ME, 
                                 demo.pars.NSPN, 
                                 var,
                                 mycovars,
                                 maxcomp,
                                 N_FOLDS,
                                 NPROCS,
                                 NITER,
                                 modules_file,
                                 to_select)
  }
  
  results[[i]]$to_select = to_select
  print(mean(results[[i]]$rho))
  return(results)
}


select_modules_deprec = function(modules_file, to_select, internal_only = F){
  
  modules = unlist(read.table(modules_file))
  n.rois = length(modules)
  adj = matrix(F, n.rois, n.rois)
  for (i in to_select){
    if (internal_only){
      adj[modules == i, ] = T
      adj[modules == i, ] = T
    } else {
      adj[modules == i, modules == i] = T
    }
  }
  valid = which(squareform2(adj))
  print(paste("Connections in module(s) ", to_select, ": " , length(valid)))
  return(valid)
}


