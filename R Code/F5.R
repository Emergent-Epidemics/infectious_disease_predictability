#SV Scarpino
#F5 for https://arxiv.org/abs/1703.07317

#set working dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to source file location

#libraries (not included in limits_acc_functions.R)
library(entropy)

###########
#acc funcs#
###########
source("limits_acc_functions.R")

######
#data#
######
static <- read.csv("~/Dropbox (EmergentEpidemicsLab)/Active Projects/Limits to prediction/exp_static_simulation_data_bulk_100.csv", header = TRUE)
data <- list()
for(i in 2:ncol(static)){
  data[[i]] <- static[,i]
}

t <- 0.01

##########
#Analysis#
##########

how_many_same <- rep(NA, length(data))
ents <- how_many_same
pb <- txtProgressBar(1, length(data), style=3)
for(i in 2:length(data)){
  x.i <- as.numeric(data[[i]])
  
  first.numb <- min(which(is.na(x.i) == FALSE))
  
  if(is.finite(first.numb) == FALSE){
    next
  }
  
  x.i <- x.i[first.numb:length(x.i)]
  last.numb <- max(which(is.na(x.i) == FALSE))
  x.i <- x.i[1:last.numb]
  
  if(length(x.i) < 10){
    next
  }
  NAs <- c(1, which(is.na(x.i) == TRUE), length(x.i))
  if(length(NAs) == 1 | length(NAs) == length(x.i)){
    next
  }
  
  if(length(x.i) < 10){
    next
  }
  
  perm_ent <- rel_ent(x = as.numeric(x.i), m.min = 4, m.max = 4, t.min = 1, t.max = 1)
  ents[i] <- perm_ent$ent
  
  do_test <- FALSE
  if(do_test == TRUE){
    test2 <- fit$results
    kls <- test2$kl
    kls[which(test2$kl.p > 0.05)] <- NA
    
    length(which(test2$kl.p > 0.05))/nrow(test2)
    
    #quartz()
    #plot(test2$raw.perm.entropy, type = "l", cex = 2)
    #par(new = TRUE)
    #plot(x.i, type = "l", lty = 3, col = "red", yaxt = "n", ylab = "")
    #axis(4)
    
    window <- 75
    fit <- full_lenth_pred_window(data = as.numeric(x.i), start = 1, end = length(x.i), window = window, d_1 = 4, d_2 = 4)
    
    if(length(fit$dists) == 0){
      next
      setTxtProgressBar(pb, i)
    }
    
    full_dist <- codebook(as.numeric(x.i), m = 4)
    
    ps <- rep(NA, length(fit$dists))
    for(dist_d in 1:length(fit$dists)){
      p.i <- comp_dist(dist_1 = full_dist, dist_2 = fit$dists[[dist_d]], size = window)
      ps[dist_d] <- p.i
    }
    how_many_same[i] <- length(which(ps > 0.05))/dist_d
  }
  setTxtProgressBar(pb, i)
}
stats <- read.csv("~/Dropbox (EmergentEpidemicsLab)/Active Projects/Limits to prediction/exp_static_simulation_stats_bulk_100.csv", header = TRUE)
#tc <- as.numeric(as.character(stats[1,-1]))/(as.numeric(as.character(stats[2,-1]))-as.numeric(as.character(stats[1,-1])))
tc <- as.numeric(as.character(stats[3,-1]))
#mu_ent <- by(ents[-1], tc, median, na.rm = TRUE)
#tc_plot <- as.numeric(names(mu_ent))
#plot(tc_plot, mu_ent)
#boxplot(ents[-1] ~ tc)
dat.plot <- data.frame(tc, ents[-1])
colnames(dat.plot) <- c("tc", "PE")
p <- ggplot(data = dat.plot, aes(x = tc - t, y = PE, colour = tc))

quartz(width = 7, height = 5)
p + geom_point() + geom_smooth(method = "loess")+ xlab("Distance from an epidemic") + ylab("Permutation entropy") + theme(legend.position = "none", panel.background = element_rect(fill = "#f0f0f0", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.005,0.005), limits = c(0.5, 0.85)) + scale_x_continuous(expand = c(0.001,0.001))
