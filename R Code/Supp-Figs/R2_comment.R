#SV Scarpino
#R2 comment for https://www.nature.com/articles/s41467-019-08616-0

#libraries (not included in limits_acc_functions.R)
library(EpiModel)
library(wesanderson)

###########
#acc funcs#
###########
source("../limits_acc_functions.R")

###############
#Global Params#
###############
n_iter <- 1000
run_new <- FALSE #set to TRUE to re-run the iterations
cols <- c("#2166ac", "#e0e0e0", "#4d4d4d","#b2182b")

######
#Data#
######
data <- list()

noise <- rnorm(n = 5000, mean = 0, sd = 1)
noise_10_0 <- noise
noise_10_100 <- noise
noise_10_m100 <- noise
noise_10_0[seq(1, length(noise_10_0), by = 10)] <- 0
noise_10_100[seq(1, length(noise_10_0), by = 10)] <- 100
noise_10_m100[seq(1, length(noise_10_0), by = 10)] <- -100

data[["noise-noise"]] <- noise
data[["noise-every_ten_0"]] <- noise_10_0
data[["noise-every_ten_100"]] <- noise_10_100
data[["noise-every_ten_m100"]] <- noise_10_m100

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
      ent_raw.i <- entropy(codebook(x.i, m = ent.i$d, t = ent.i$tau))/log(factorial(ent.i$d))
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
  save(RESULTS, file = paste0(time_stamp, "_R2_noise_sim.RData"))
}else{
  load("1539035610.8112_R2_noise_sim.RData")
}

RESULTS$transmissibility <- as.character(RESULTS$transmissibility)
lvls <- c("every_ten_m100", "noise", "every_ten_0", "every_ten_100")
RESULTS$trans_factor <- factor(RESULTS$transmissibility, levels = lvls)

p1 <- ggplot(data = RESULTS, aes(x = trans_factor, y = 1-weighted.perm.entropy, fill = trans_factor))
quartz(width = 10, height = 6)
p1 + geom_boxplot(coef = 100) + scale_fill_manual(values = cols,  guide_legend(title = "Model"))  + xlab("Transmissibility") + ylab("Predictability (1 - H)") + theme(legend.position = c(0.1, 0.85), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(limits = c(0,1), expand = c(0.001,0.001))
