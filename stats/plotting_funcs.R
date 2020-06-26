
###############################################################################
# Functions for plotting results
# benjamin.garzon@gmail.com
###############################################################################

mypackages = c("ggplot2", "dplyr", "reshape2", "psych", "pracma", "mclust")
#library(knitr)
#library(xtable)
#library(dplyr)

lapply(mypackages, function (x) if (!require(x, character.only = T)) install.packages(x) else print(paste(x, "already installed!")))


###############################################################################
# Create data frame with results for all networks
###############################################################################

statstodf = function(myresults, group_name, check_pars = T, alpha = 0.05){
  print("--------------------------")
  if (check_pars)  check_parameters(myresults)  

  module_names = sapply(myresults, function(x) x$label)

  rho = data.frame(
      mean = sapply(myresults, 
                    function(x) ifelse(!x$empty, x$cor.test$estimate, NA)), 
      # low = sapply(myresults, function(x) ifelse(!x$empty, 
      #quantile(x$cor.test$perm, 0.025), NA)),
      # high = sapply(myresults, function(x) ifelse(!x$empty, 
      #quantile(x$cor.test$perm, 0.975), NA)),
      low = sapply(myresults, 
                    function(x) ifelse(!x$empty, min(x$cor.test$perm), NA)),
      high = sapply(myresults, 
                    function(x) ifelse(!x$empty, 
                                       quantile(x$cor.test$perm, 0.95), NA)),
      p = sapply(myresults, 
                 function(x) ifelse(!x$empty, x$cor.test$p.value, NA)),
      p.approx = sapply(myresults, 
                        function(x) ifelse(!x$empty, 
                                           x$cor.test$p.value.approx, NA)),
      nfeatures = sapply(myresults, 
                         function(x) ifelse(!x$empty, x$nfeatures, NA)),
      from = as.factor(sapply(myresults, 
                              function(x) ifelse(!x$empty, x$label, NA))),
      time = sapply(myresults, function(x) ifelse(!x$empty, x$time, NA))/3600
    ) 
  rho$group_name = group_name
  
  cor.perm = sapply(myresults, function(x) x$cor.test$perm)
  R2.perm = sapply(myresults, function(x) x$R2.test$perm)
  
  # get max statistic and corresponding p-values
  nullmax = apply(cor.perm, 1, max)
  rho$pvalues.FWE = sapply(cor.perm[1, ], 
                           function(x) sum(x <= nullmax)/length(nullmax))
  
  # use the approximate p-values
  z.perm = fisherz(cor.perm)
  nullmax = apply(z.perm, 1, max)
  
  mysigma = apply(z.perm[-1, ], 2, function(x) sqrt(sum(x^2)/length(x)))
  rho$low = -1 #fisherz2r(2*mysigma) 
  rho$high = fisherz2r(qnorm(0.95, sd = mysigma)) 
  rho$p.approx = 1 - pnorm(z.perm[1, ]/mysigma)
  rho$sigma = mysigma
  
  #browser()
  # fit a Gaussian mixture to the distribution of the max
  fit = Mclust(nullmax[-1], G=2, model="V")
  mymu = fit$parameters$mean
  mysigma = sqrt(fit$parameters$variance$sigmasq)
  mypro = fit$parameters$pro
  
  plot(fit, what="density", main="", xlab="r")
  hist(nullmax[-1], 10, add = T, col = rgb(1, 0, 0, 0.5), probability = T)
  rho$pvalues.FWE.approx = 
  0.5*(mypro[1]*(1 - pnorm((z.perm[1, ] - mymu[1])/mysigma[1])) + 
  mypro[2]*(1 - pnorm((z.perm[1, ] - mymu[2])/mysigma[2]))
  )
  
  qqnorm((nullmax[-1] - mean(nullmax[-1]))/sd(nullmax[-1]))
  abline(0, 1)
  
  # check the prediction is not driven by the number of features
  print("Correlation between number of features and mean rho (abs and original)")
  print(cor.test(rho$nfeatures, abs(rho$mean)))  
  print(cor.test(rho$nfeatures, rho$mean), method = "spearman")
  par(mfrow=c(1,1))
  plot(rho$nfeatures, rho$mean)
  
  rho$width = ifelse(rho$p < alpha, 1, 0.5)
  rho$size = ifelse(rho$p < alpha, 5, 4)
  
  rho$network = as.factor(rho$from) 
  ps = rho$p.approx
  ps.adjusted = p.adjust(ps, method = "fdr")
  print("Adjusted correlation all features")
  rho$fdrp = ps.adjusted 
  rho$sig = ifelse(rho$p.approx < alpha, "*",  "")
  rho$sig[rho$fdrp < alpha] = "**"
  #rho$sig.FWE = ifelse(rho$p.approx < alpha, "*",  "")
  #rho$sig.FWE[rho$pvalues.FWE.approx < alpha] = "**"
  
    #rho$sig[rho$pvalues.FWE < alpha] = "**"

  rho = rho %>% arrange(network)
  print(rho)

  # Check how long it took  
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

###############################################################################
# Plot correlation coefficients with intervals for null distribution
###############################################################################

plot_intervals = function(RESULTS_FILE, measure, myfigname, selection = NULL){
  
  load(RESULTS_FILE)
  var.network.results = results
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
  
  myplot = 
    ggplot(rho, aes( x = network, y = mean, ymin = high, ymax = high, 
                     label = rho$sig)) + 
    geom_errorbar(position = position_dodge(0.4), size = 2, linetype = 1, 
                  width = 0.4) + 
    geom_linerange(aes(ymin = -0.25, ymax = high), width = 0.4) + 
    ggtitle(paste('Correlation between observed and predicted', 
                  measure, "\n")) + 
    xlab('Network') + 
    geom_point(position = position_dodge(0.4), size = 15, shape = 19) + 
    geom_text(aes( x = network, y = .25, label = rho$sig), 
              position = position_dodge(0.4), size = CEX_TEXT) +
    ylim(-.25, .25) + 
    theme_bw() + 
    xlab('\nModule') + 
    ylab('r') +
    theme(
      legend.position = "none",
      axis.text = element_text(size = CEX_AXIS - 10),
      axis.text.x = element_text(angle = 45, hjust = 1),# vjust = -1.5, 
      axis.text.y = element_text(hjust = -.2),
      axis.title = element_text(size = CEX_AXIS, vjust = 2),
      plot.title = element_text(hjust = 0.5, vjust = -5, 
                                size = CEX_TITLE - 10),
      axis.line = element_line(size = 2),
      axis.ticks.length = unit(20, "pt"),
      axis.ticks = element_line(size = 2)
    )
  
  save_fig(figname = myfigname, res = BWRES)
  print(myplot)
  dev.off()
  return(rho)
}


###############################################################################
# Some parameter checks
###############################################################################

check_parameters = function(myresults){
  
  labels = as.factor(sapply(myresults, 
                            function(x) ifelse(!x$empty, x$label, NA)))
  nfeatures.perm = sapply(myresults, function(x) colMeans(x$nfeatures.perm))
  K.iter = sapply(myresults, function(x) x$K.iter)
  eta.iter = sapply(myresults, function(x) x$eta.iter)
  colnames(nfeatures.perm) = colnames(K.iter) = colnames(eta.iter) = labels
  
  # motion parameters
  print("Model parameters")
  
  print(ggplot(melt(K.iter), aes(x = value)) + geom_density() +
          facet_grid(. ~ Var2))
  print(ggplot(melt(eta.iter), aes(x = value)) + geom_density() + 
          facet_grid(. ~ Var2))
  print(ggplot(melt(nfeatures.perm), aes(x = Var1, y = value)) + geom_line() + 
          geom_point() + facet_grid(. ~ Var2))
  
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

###############################################################################
# Save figure
###############################################################################
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
    jpeg(figname, width = floor(res/2.54*width), 
         height = floor(res/2.54*height), 
         pointsize=POINTSIZE)
  } else {
    figname = file.path(FIGS_DIR, paste0( figname, '.png'))
    png(figname, width = floor(res/2.54*width), 
        height = floor(res/2.54*height), 
        pointsize=POINTSIZE)
  }
  
}

###############################################################################
# Older plotting functions
###############################################################################
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
  #  for (i in seq(length(var.network.results))) var.network.results[[i]]$label = paste0("net", sprintf("%02d", i))
  
  rho = statstodf(var.network.results, group_name  = "") %>% arrange(network) 
  
  cor.perm = sapply(var.network.results, function(x) x$cor.test$perm)[-1, ] 
  colnames(cor.perm) = sapply(var.network.results, function(x) x$label)
  
  cor.perm.melt = melt(cor.perm, value.name = "cor", varnames = c("permutation", "network")) %>% 
    mutate(network = factor(network, levels = rho$network)) %>% arrange(network)
  myplot = 
    ggplot(data = cor.perm.melt, aes(x = network, y = cor)) +
    geom_violin(adjust = 0.5, size = 2) +
    geom_point(data = rho, aes( x = network, y = mean), position = position_dodge(0.4), size = 15, shape = 19) +
    geom_text(data = rho, aes( x = network, y = .25, label = rho$sig), position = position_dodge(0.4), size = CEX_TEXT) +
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

