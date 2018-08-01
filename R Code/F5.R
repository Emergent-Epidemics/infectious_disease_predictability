#SV Scarpino
#F5 for https://arxiv.org/abs/1703.07317

#set working dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to source file location

#libraries (not included in limits_acc_functions.R)
library(EpiModel)
library(statcomp)

###########
#acc funcs#
###########
source("limits_acc_functions.R")

###############
#Global Params#
###############
n_iter <- 1000
run_new <- TRUE #set to TRUE to re-run the iterations
do_coarse_grain <- FALSE #set TRUE to coarse grain
grain_by <- 10 #time steps to coarse grain over

######
#Data#
######
irr_files <- list.files("../Data/irr-sims2/", full.names = TRUE)
reg_files <- list.files("../Data/reg-sims2/", full.names = TRUE)

data <- list()

for(i in irr_files){
  irr.i <- read.csv(i, header = FALSE)
  name.i <- paste0("irr-0.", strsplit(i, "[.]")[[1]][4])
  
  if(do_coarse_grain == TRUE){
    irr.i <- course_grain(x = irr.i[,1], grain.by = grain_by)
  }else{
    irr.i <- irr.i[,1]
  }
  
  data[[name.i]] <- irr.i
}

for(i in reg_files){
  reg.i <- read.csv(i, header = FALSE)
  name.i <- paste0("reg-0.", strsplit(i, "[.]")[[1]][4])
  
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
      
      ent.i <- rel_ent(x = x.i, m.min = 2, m.max = 7, t.min = 1, t.max = 1, do_mc_ent = FALSE)
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
  load("../Results/1533140782.44114_big_results_PE_SIM.RData")
}

RESULTS$transmissibility <- as.numeric(as.character(RESULTS$transmissibility))
RESULTS$trans_factor <- round(RESULTS$transmissibility, 2)
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
p0 <- ggplot(data = RESULTS, aes(x = n, y = 1-raw.perm.entropy, colour = transmissibility, group = disease))

quartz(width = 10, height = 6)

p0 + geom_point() + xlab("Number of weeks") + ylab("Predictability (1 - H)") + theme(legend.position = c(0.86, 0.73), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01)) + scale_x_continuous(expand = c(0.01,0.01), limits = c(0, 100)) + geom_smooth()

p1 <- ggplot(data = RESULTS, aes(x = trans_factor, y = 1-raw.perm.entropy, fill = disease))
p1+geom_boxplot() + xlab("Transmissibility") + ylab("Predictability (1 - H)") + theme(legend.position = c(0.86, 0.73), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3))
