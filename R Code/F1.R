#SV Scarpino
#F1 for https://www.nature.com/articles/s41467-019-08616-0

#libraries (not included in limits_acc_functions.R)
library(wesanderson)

###########
#acc funcs#
###########
source("limits_acc_functions.R")

######
#Data#
######

#measles
meas <- read.csv("../Data/MEASLES_Cases_1909-2001_20150923120449.csv", na.strings = "NA", header = TRUE, stringsAsFactors = FALSE)

#adding comparison time series
data <- list()
meas_TX <- meas$TEXAS[which(meas$YEAR < 1965)] #pre-vaccine
meas_TX_filt <- filt_lead_trail_NA(meas_TX)
data[["Measles-Texas"]] <- meas_TX_filt$x
years <- meas$YEAR[which(meas$YEAR < 1965)][meas_TX_filt$first:(meas_TX_filt$first+meas_TX_filt$last)]
n <- length(data[["Measles-Texas"]])
data[["noise-full"]] <- rnorm(n)
data[["noiseMissing-missing"]] <- data[["noise-full"]]
data[["noiseMissing-missing"]][which(is.na(data[["Measles-Texas"]]))] <- NA #exactly the noise model from the data
data[["sine-noise"]] <- sin(1:n) + rnorm(n, 0, 0.001)

########
#Params#
########
window <- 52 #10 times 3!, which is greatest dimension allowed

##########
#Analysis#
##########
results <- list()
plot_data <- list()
pb <- txtProgressBar(1, length(data), style=3)
for(i in 1:length(data)){
  x.i <- data[[i]]
  
  fit.i <- full_lenth_pred_window(data = x.i, start = 1, end = length(x.i), window = window, d_1 = 2, d_2 = 5)
  
  results[[names(data)[i]]] <- fit.i
  plot_data[[names(data)[i]]] <- x.i
  
  setTxtProgressBar(pb, i)
}

##########
#Plotting#
##########
cols <- c("#7A8A83", "#ECBBAD", "#9D5241", "#334B58", "#7A8A83")
quartz(width = 15, height = 8)
par(mar = c(5,4,4,4))
for(i in 1:length(data)){
  plot(results[[i]]$results$raw.perm.entropy, type = "l", yaxt = "n", xaxt = "n", ylim = c(0.5,1), lwd = 3, col = cols[i], bty = "n", xlab = "", ylab = "")
  par(new = TRUE)
}
plot(data[["Measles-Texas"]], col = cols[i+1], lty = 3, lwd = 2, type = "l", yaxt = "n", xaxt = "n", bty = "n", xlab = "Year", ylab = "Permutation Entropy")

at.x <- seq(0, length(data[["Measles-Texas"]]), length.out = 5)
axis(1, at = at.x, labels = years[at.x+1])
at.y.z <- seq(min(data[["Measles-Texas"]], na.rm = TRUE), max(data[["Measles-Texas"]], na.rm = TRUE), length.out = 4)
axis(2, at = at.y.z, labels = round(seq(0.5, 1, length.out = 4), 2), las = 2)
axis(4, at = at.y.z, labels = floor(at.y.z))
mtext("Weekly cases", side=4, line=3)
legend(-50, 4000, legend = c(names(data), "Measles-Texas"), col = c(cols, cols[i+1]), lwd = 3, lty = c(1,1,1,1,3), bty = "n")


ps <- rep(NA, length(results[["Measles-Texas"]]$dists))
for(dist_d in 2:length(ps)){
  comp.i <- comp_dist(dist_1 = results[["Measles-Texas"]]$dists[[(dist_d-1)]], dist_2 = results[["Measles-Texas"]]$dists[[dist_d]], size = window)
  ps[dist_d] <- comp.i$p
}
length(which(ps < 0.05))/length(na.omit(ps))
