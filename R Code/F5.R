#SV Scarpino
#F5 for https://arxiv.org/abs/1703.07317

#set working dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to source file location

#libraries (not included in limits_acc_functions.R)
library(EpiModel)
library(statcomp)
library(wesanderson)
library(gridExtra)
library(cowplot)

###########
#acc funcs#
###########
source("limits_acc_functions.R")

###############
#Global Params#
###############
n_iter <- 1000
run_new <- FALSE #set to TRUE to re-run the iterations
do_coarse_grain <- FALSE #set TRUE to coarse grain
grain_by <- 10 #time steps to coarse grain over
cols <- wes_palette(name = "Darjeeling2", n = 2)

######
#Data#
######
irr_files <- list.files("../Data/irr-sims2/", full.names = TRUE)
reg_files <- list.files("../Data/reg-sims2/", full.names = TRUE)

data <- list()

for(i in irr_files){
  irr.i <- read.csv(i, header = FALSE)
  trans.i <- strsplit(i, "clique-size")[[1]][2]
  trans.i <- strsplit(trans.i, ".txt")[[1]][1]
  name.i <- paste0("irr-", trans.i)
  
  if(do_coarse_grain == TRUE){
    irr.i <- course_grain(x = irr.i[,1], grain.by = grain_by)
  }else{
    irr.i <- irr.i[,1]
  }
  
  data[[name.i]] <- irr.i
}

for(i in reg_files){
  reg.i <- read.csv(i, header = FALSE)
  trans.i <- strsplit(i, "clique-size")[[1]][2]
  trans.i <- strsplit(trans.i, ".txt")[[1]][1]
  name.i <- paste0("reg-", trans.i)
  
  if(do_coarse_grain == TRUE){
    reg.i <- course_grain(x = reg.i[,1], grain.by = grain_by)
  }else{
    reg.i <- reg.i[,1]
  }
  
  
  data[[name.i]] <- course_grain(x = reg.i, grain.by = 10)
}


##########
#Analysis#
##########
if(run_new == TRUE){
  names <- unlist(strsplit(names(data), "-"))
  disease <- names[seq(1, length(names), by = 2)]
  transmissibility <- names[seq(2, length(names), by = 2)]
  pdc.perm.entropy <- rep(NA, length(data))
  n <- rep(NA, length(data))
  d <- rep(NA, length(data))
  tau <- rep(NA, length(data))
  raw.perm.entropy <- rep(NA, length(data))
  boot.d <- rep(NA, length(data))
  boot.tau <- rep(NA, length(data))
  filt.perm.entropy <- rep(NA, length(data))
  n.filt <- rep(NA, length(data))
  boot.n <- rep(NA, length(data))
  weighted.perm.entropy <- rep(NA, length(data))
  total.cases <- rep(NA, length(data))
  
  results <- data.frame(disease, transmissibility, pdc.perm.entropy, n, d, tau, boot.d, boot.tau, filt.perm.entropy, n.filt, boot.n, weighted.perm.entropy, total.cases)
  
  #looping
  pb <- txtProgressBar(1, length(data), style=3)
  for(iter in 1:n_iter){
    for(i in 1:length(data)){
      x.i <- as.numeric(data[[i]])
      #write.csv(x.i, file = paste0("data/", names(data)[[i]], ".csv"), quote = FALSE)
      
      # if(course.grain == TRUE){
      #   x.i <- course_grain(x = x.i, grain.by = grain.by)
      # }
      
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
      
      #if(iter == 1){
      # pdf(paste0(names(data)[i], ".pdf"))
      # plot(x.i, type = "l", main = names(data)[i], xlab = "Time", ylab = "Cases", bty = "n")
      # dev.off()
      #}
      
      if(length(x.i) < 10){
        next
      }
      
      start.i <- sample(1:(length(x.i)-10), 1)
      length.i <- sample(11:100, 1)
      if(length(x.i[-c(1:start.i)]) < length.i){
        x.i <- x.i[start.i:length(x.i)]
      }else{
        x.i <- x.i[start.i:(start.i+length.i)]
      }
      
      ent.i <- rel_ent(x = x.i, m.min = 4, m.max = 4, t.min = 1, t.max = 1, do_mc_ent = FALSE)
      if(is.na(ent.i$d) == TRUE) next
      ent_w.i <- entropy(weighted_ordinal_pattern_distribution(x.i, ent.i$d))/log(factorial(ent.i$d))
      ent_raw.i <- entropy(codebook(x.i, m = ent.i$d, t = 1))/log(factorial(ent.i$d))
      results[i,"pdc.perm.entropy"] <- ent.i$ent
      results[i,"raw.perm.entropy"] <- ent_raw.i
      results[i,"weighted.perm.entropy"] <- ent_w.i
      results[i,"n"] <- ent.i$n
      results[i,"d"] <- ent.i$d
      results[i, "tau"] <- ent.i$tau
      results[i, "boot.d"] <- ent.i$rel.ent.d
      results[i, "boot.tau"] <- ent.i$rel.ent.tau
      results[i, "total.cases"] <- sum(x.i, na.rm = TRUE)
    }
    results$disease <- as.character(results$disease)
    results$transmissibility <- as.character(results$transmissibility)
    if(iter == 1){
      RESULTS <- results
    }else{
      RESULTS <- rbind(RESULTS, results)
    }
    setTxtProgressBar(pb, iter)
  }
  
  time_stamp <- as.numeric(Sys.time())
  save(RESULTS, file = paste0(time_stamp, "_big_results_PE_SIM.RData"))
}else{
  load("../Results/1533154427.72181_big_results_PE_SIM.RData")
}

RESULTS$transmissibility <- as.numeric(as.character(RESULTS$transmissibility))
RESULTS$trans_factor <- round(RESULTS$transmissibility, 1)
RESULTS$trans_factor[which(RESULTS$trans_factor == 0.11)] <- 0.12
lvls <- unique(RESULTS$trans_factor)
lvls <- lvls[order(lvls, decreasing = FALSE)]
RESULTS$trans_factor <- factor(RESULTS$trans_factor, levels = lvls)

#######
#Model#
#######
mod <- lm((1-raw.perm.entropy) ~ n*disease*transmissibility, data = RESULTS)
summary(mod)

#########
#Figures#
#########
p0 <- ggplot(data = RESULTS, aes(x = n, y = 1-raw.perm.entropy, colour = disease))
p0 <- p0 + geom_point() + scale_colour_manual(values = paste0(cols, 75),  guide_legend(title = "")) + xlab("Number of weeks") + ylab("Predictability (1 - H)") + theme(legend.position = c(0.0, 0.85), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01)) + scale_x_continuous(expand = c(0.01,0.01), limits = c(0, 100)) + geom_smooth(inherit.aes = FALSE, aes(x = n, y = 1-raw.perm.entropy, group = disease, linetype = disease), color = "black", se = FALSE) +scale_linetype(guide_legend(title = " ")) 

#p0 + geom_smooth(method = "loess") + scale_colour_manual(values = paste0(cols, 75),  guide_legend(title = "Model")) + xlab("Number of weeks") + ylab("Predictability (1 - H)") + theme(legend.position = c(0.05, 0.9), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01)) + scale_x_continuous(expand = c(0.01,0.01), limits = c(0, 100)) 

use.plot <- which(RESULTS$transmissibility <= 1.6)
plot.df <- RESULTS[use.plot, ]
p1 <- ggplot(data = plot.df, aes(x = trans_factor, y = 1-raw.perm.entropy, fill = disease))
p1 <- p1 + geom_boxplot(coef = 100)+ scale_fill_manual(values = cols,  guide_legend(title = ""))  + xlab("Transmissibility") + ylab("Predictability (1 - H)") + theme(legend.position = c(0.0, 0.90), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(limits = c(0,1), expand = c(0.001,0.001)) 

quartz(width = 15, height = 6)
plot_grid(p0, p1, nrow = 1, labels = "AUTO")
