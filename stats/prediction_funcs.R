######## Prediction functions ##########################

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
library(RVAideMemoire)
#library(ppcor)


MYPALETTE=brewer.pal(11, "Spectral")[seq(11, 1, -1)]# colorspace::diverge_hsv(15, power=2)

FINELINE = 0.5
GGTEXTSIZE2 = 10
GGTEXTSIZE3 = 10
NDIGITS = 3

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
    
  }  
  
  valid = which(squareform2(adj))
  print(paste("Connections between module(s) ", to_select[1], ", " , to_select[2], ": " , length(valid), "Inverted: ", invert))
  return(valid)
}

###############################################
# run prediction for different partitions
###############################################
do_predictions_loop = function(DIR_ICA_ME, # working dir
                               demo.pars, # predictors
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
                                 demo.pars, 
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



###############################################
# fit and test predictive model
###############################################
do_prediction = function(DIR_ICA_ME, # working dir
                         demo.pars, # predictors
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
  mydata = merge(ica_data, demo.pars, by="Subject")
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


###############################################
# mediation analysis
###############################################
do_mediation = function(DIR_ICA_ME, 
                        results, 
                        NSIMS = 2000,
                        write_coefs = F,
                        RESULTS_DIR = file.path(DIR_ICA_ME, 'mediation/'), 
                        controlIQ = F
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
  
  if (controlIQ){
    # 1) X -> Y
    model.0 = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU + IQcomp, data = data.M)  
    
    # 2) X -> M
    model.M = lm(M ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU + IQcomp, data = data.M)
    
    # 3) X + M -> Y
    model.Y = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU + M + IQcomp, data = data.M)  
  } else {
    # 1) X -> Y
    model.0 = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU, data = data.M)  
    
    # 2) X -> M
    model.M = lm(M ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU, data = data.M)
    
    # 3) X + M -> Y
    model.Y = lm(y ~ age.1*sex + brain_vol + fd + meica_dof + UCL + CBU + M, data = data.M)  
  }
  
  print(summary(model.0))
  print(summary(model.M))
  print(summary(model.Y))
  
  # mediation 
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

###############################################
# compute results but corrected for covariates
###############################################
do_correction = function(DIR_ICA_ME, 
                         results, 
                         demo.pars,
                         controlfor = NULL)
{ 
  data.M = results$covars.exp
  data.M$y.obs = results$y.test.mean
  data.M$y.pred = results$y.pred.mean
  data.M$Subject = results$Subject
  data.all = merge(data.M, demo.pars, by = "Subject", suffixes = c("", ".y"), no.dups = T)
  data.all$age.1 = scale(data.all$age)
  data.all$age.1.sex = data.all$age.1*data.all$sex.num
  
  results$cor.test.twotailed = pcor.test(data.M$y.obs, data.M$y.pred, data.all[c("age.1", "sex.num", "age.1.sex", controlfor)]) 
  # results$cor.test.twotailed$conf.int = c(0, 0)
  
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
# plot results for a module combination 
###############################################
plot_combinations = function(myresults, level, module_names, measure, myfigname = 'Accuracy'){
  
  print("--------------------------")
  alpha = 0.05
  
  rho = data.frame(
    mean = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.twotailed$estimate, NA)), 
    t = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.twotailed$statistic, NA)), 
    low = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.twotailed$conf.int[1], NA)),
    high = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.twotailed$conf.int[2], NA)),
    p = sapply(myresults, function(x) ifelse(!x$empty, x$cor.test.twotailed$p.value, NA)),
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
  print(cor.test(rho$nfeatures, rho$mean))#, alternative = "greater"))  
  
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
  
  save_fig(figname = myfigname, res = BWRES)
  rho.diag = subset(rho, from == to)
  rho.diag$network = as.factor(module_names[rho.diag$from]) 
  ps = c(level$p.value, rho.diag$p)
  ps.adjusted = p.adjust(ps, method = "fdr")
  print("Adjusted correlation all features")
  print(paste(sprintf("%0.3f", level$estimate), ps.adjusted[1], level$p.value))
  rho.diag$fdrp = ps.adjusted[-1] 
  rho.diag$sig = ifelse(rho.diag$p < alpha, "*",  "")
  rho.diag$sig[rho.diag$fdrp < alpha] = "**"
  
  print(rho.diag %>% arrange(network))
  # print(
  #   ggplot(rho.diag, aes(network, mean, label = rho.diag$sig)) +
  #     ggtitle(paste('Correlation between observed and predicted', measure, "\n")) + xlab('Network') + 
  #     geom_text(size = CEX_TEXT) +
  #     geom_hline(yintercept = level$estimate, lty = 2, size = CEX_LINE) + 
  #     geom_col(col = 'black', fill = 'gray', alpha = 0.5) + ylim(0, 0.3) + 
  #     theme_bw() + 
  #     xlab('\nModule') + 
  #     ylab('r') +
  #     
  #     theme(
  #       axis.text = element_text(size = CEX_AXIS, margin = margin(t = 100, r = 0, b = 100, l = 0)),
  #       axis.title = element_text(size = CEX_AXIS, margin = margin(t = 100, r = 0, b = 100, l = 0)),
  #       plot.title = element_text(hjust = 0.5, size = CEX_TITLE),
  #       plot.margin = margin(2, 0, 2, 0, "cm"),
  #       axis.line = element_line(size = 2)
  #     )
  
  print(
    ggplot(rho.diag, aes( x = network, y = mean, ymin = low, ymax = high, label = rho.diag$sig )) + 
      geom_errorbar(size = 2, linetype = 1, width = 0.5) + 
      ggtitle(paste('Correlation between observed and predicted', measure, "\n")) + xlab('Network') + 
      geom_text(size = CEX_TEXT, nudge_y = .12) +
      geom_point(size = 15, shape = 19) + 
      geom_hline(yintercept = level$estimate, lty = 2, size = CEX_LINE) + 
      ylim(-.3, .4) + 
      theme_bw() + 
      xlab('\nModule') + 
      ylab('r') +
      theme(
        axis.text = element_text(size = CEX_AXIS - 10),
        axis.text.x = element_text(vjust = -1.5),
        axis.text.y = element_text(hjust = -.2),
        axis.title = element_text(size = CEX_AXIS, vjust = 2),
        plot.title = element_text(hjust = 0.5, vjust = -5, size = CEX_TITLE - 10),
#        plot.margin = margin(1, 0, 2, 0, "cm"),
        axis.line = element_line(size = 2),
        axis.ticks.length = unit(20, "pt"),
        axis.ticks = element_line(size = 2)
      )
  ) 
  dev.off()
  return(rho)
}


