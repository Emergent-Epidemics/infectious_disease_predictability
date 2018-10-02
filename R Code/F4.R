#SV Scarpino
#F4 for https://arxiv.org/abs/1703.07317

#set working dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to source file location

#libraries (not included in limits_acc_functions.R)
library(entropy)

###########
#acc funcs#
###########
source("limits_acc_functions.R")

###############
#Global Params#
###############
run_new_thresh <- FALSE #this will take a long time and isn't needed unless the data set or distance metric(s) change.
vaccine_switch_year <- 1963

D1 <- 3
D2 <- 3
W <- 52

######
#data#
######

#measles
meas <- read.csv("../Data/MEASLES_Cases_1909-2001_20150923120449.csv", na.strings = "NA", header = TRUE, stringsAsFactors = FALSE)

data <- list()
meas_TX <- meas$TEXAS
meas_TX_filt <- filt_lead_trail_NA(meas_TX)
data[["Measles-Texas"]] <- meas_TX_filt$x
years <- meas$YEAR[meas_TX_filt$first:(meas_TX_filt$first+meas_TX_filt$last)]

##########
#Analysis#
##########
fit <- full_lenth_pred(data = data[["Measles-Texas"]], start = 1, end = length(data[["Measles-Texas"]]), d_1 = D1, d_2 = D2)

fit_window <- full_lenth_pred_window(data = data[["Measles-Texas"]], start = 1, end = length(data[["Measles-Texas"]]), d_1 = D1, d_2 = D2, window = W)

sine_wave <- as.numeric(sin(1:length(data[["Measles-Texas"]])))
white_noise <- rnorm(length(data[["Measles-Texas"]]))

fit_window_sine <- full_lenth_pred_window(data =  sine_wave, start = 1, end = length(sine_wave), d_1 = D1, d_2 = D2, window = W)
fit_window_noise <- full_lenth_pred_window(data =  white_noise, start = 1, end = length(white_noise), d_1 = D1, d_2 = D2, window = W)

jls_ts_meas <- ts_JS(fit_window = fit_window)
jls_ts_sine <- ts_JS(fit_window = fit_window_sine)
jls_ts_noise <- ts_JS(fit_window = fit_window_noise)

dat.plot <- data.frame(c(jls_ts_meas$JSs, jls_ts_sine$JSs, jls_ts_noise$JSs), c(jls_ts_meas$time_dists, jls_ts_sine$time_dists, jls_ts_noise$time_dists), c(rep("measles-Texas", length(jls_ts_meas$JSs)), rep("sine-noNoise", length(jls_ts_sine$JSs)), rep("whiteNoise", length(jls_ts_sine$JSs))))
colnames(dat.plot) <- c("JS", "time", "group")

#run vaccination analysis
data_full <- list()
year_range <- c(1900,2017)
for(i in 3:ncol(meas)){
  filt.i <- which(meas$YEAR >= year_range[1] & meas$YEAR <= year_range[2])
  data_full[[paste0("Measles-",colnames(meas[i]))]] <- meas[filt.i,i]
}

fit_windows <- list()
years_full <- list()
pb <- txtProgressBar(1, length(data_full), style=3)
for(i in 1:length(data_full)){
  data.i <- filt_lead_trail_NA(data_full[[i]])
  
  if(is.na(data.i$x) == TRUE){
    next
  }
  
  years_full[[names(data_full)[i]]] <- meas$YEAR[data.i$first]:meas$YEAR[data.i$last]
  
  x.i <- diff(data.i$x)
  fit_window.i <- full_lenth_pred_window(data = x.i, start = 1, end = length(x.i), d_1 = D1, d_2 = D2, window = W)
  fit_windows[[names(data_full)[i]]] <- fit_window.i
  
  if(length(years_full[[names(data_full)[i]]]) < length(fit_window.i$dists)){
    stop()
  }
    
  setTxtProgressBar(pb, i)
}

for(i in 1:length(fit_windows)){
  JS.i <- c()
  comparison.i <- c()
  location.i <- c()
  
  for(j in 1:length(fit_windows[[i]]$dists)){
    if(j == length(fit_windows[[i]]$dists)){
      next()
    }
    loop.var <- (j+1):length(fit_windows[[i]]$dists)
    for(k in loop.var){
      fit_window.jk <- list()
      fit_window.jk$dists <- list(fit_windows[[i]]$dists[[j]], fit_windows[[i]]$dists[[k]])
      JS.jk <- ts_JS(fit_window = fit_window.jk)
      
      if(length(JS.jk$JSs) > 1){
        stop()
      }
      
      comp.i <- NA
      if(years_full[[i]][j] < vaccine_switch_year & years_full[[i]][k] < vaccine_switch_year){
        comp.i <- "pre-pre"
      }else{
        if(years_full[[i]][j] > vaccine_switch_year & years_full[[i]][k] > vaccine_switch_year){
          comp.i <- "post-post"
        }else{
          if(years_full[[i]][j] > vaccine_switch_year & years_full[[i]][k] < vaccine_switch_year | years_full[[i]][j] < vaccine_switch_year & years_full[[i]][k] > vaccine_switch_year){
          comp.i <- "pre-post"
        }
        }
      }

      JS.i <- c(JS.i, JS.jk$JSs)
      comparison.i <- c(comparison.i, comp.i)
      location.i <- c(location.i, names(fit_windows)[[i]])
    }
  }
  
  dat.i <- data.frame(JS.i, comparison.i, location.i)
  colnames(dat.i) <- c("JS", "Comparison", "Location")
  
  dat.i$JS <- as.numeric(as.character(dat.i$JS))
  dat.i$Comparison <- as.character(dat.i$Comparison)
  dat.i$Location <- as.character(dat.i$Location)
  
  if(i == 1){
    vaccine_switch <- dat.i
  }else{
    vaccine_switch <- rbind(vaccine_switch, dat.i)
  }
}

vaccine_switch <- na.omit(vaccine_switch)
vaccine_switch$Comparison <- factor(vaccine_switch$Comparison, levels = c("pre-pre", "pre-post", "post-post"))

#p_vacc <- ggplot(data = vaccine_switch, aes(x = Comparison, y = log(JS), fill = Location))
#p_vacc + geom_boxplot() + theme(legend.position = "none")

p_vacc <- ggplot(data = vaccine_switch, aes(x = Comparison, y = JS, fill = Comparison))
quartz(width = 7, height = 5)
p_vacc + geom_boxplot(notch = TRUE, coef = 100) + scale_fill_manual(values = c("#d4c5b2","#c8573b","#77674f"), guide_legend(title = "Comparison")) + xlab("Comparison") + ylab("Jensen-Shannon divergence (log2)") + theme(legend.position = c(0.7, 0.75), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "#f0f0f000", colour = "black"), axis.text.y = element_text(colour = "black", size = 12), axis.text.x = element_text(colour = "black", size = 12), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.015,0.015), limits = c(0, 0.65)) + scale_x_discrete(expand = c(0.01,0.01)) 

rm <- which(vaccine_switch$JS == 0)
test <- aov(log(JS[-rm])~Comparison[-rm]*Location[-rm], data = vaccine_switch)
Tukey_test <- TukeyHSD(test)
Tukey_test$`Comparison[-rm]`

test2 <- aov(JS~Comparison*Location, data = vaccine_switch)
Tukey_test2 <- TukeyHSD(test2)
Tukey_test2$`Comparison`

#compute distance threshold
if(run_new_thresh == TRUE){
  threshs <- c()
  samps <- 1:length(fit_window[[2]])
  for(n in 1:1000){
    samp.n.1 <- sample(samps, 1)
    samp.n.2 <- sample(samps[-samp.n.1], 1)
    dist.i <- try(comp_dist(dist_1 = fit_window[[2]][[samp.n.1]], dist_2 = fit_window[[2]][[samp.n.2]], size = 52), silent = TRUE)
    if(is(dist.i)[1] != "try-error"){
      threshs <- c(threshs, dist.i$cutoff)
    }
  }
  write.csv(threshs, file = paste0(as.numeric(Sys.time()), "_JS_cutoffs.csv"), row.names = FALSE, quote = FALSE)
}else{
  threshs <- read.csv("../Results/1532108814.58908_JS_cutoffs.csv")
}

#supp figure
dat.plot <- dat.plot[which(dat.plot$time<66), ] #truncating to 65 years, data are very sparse after
dat.plot$group <- factor(as.character(dat.plot$group), levels = c("measles-Texas", "whiteNoise", "sine-noNoise"))

p <- ggplot(dat.plot, aes(x = as.factor(time), y = JS, fill = group))

quartz(width = 15, height = 6)
p + geom_boxplot(coef = 100) + scale_fill_manual(values = c("#335b8d90", "#ddb868", "#948580"), guide_legend(title = "Time series")) + xlab("Years") + ylab("Jensen-Shannon divergence (log2)") + theme(legend.position = c(0.05, 0.88), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "#f0f0f000", colour = "black"), axis.text.y = element_text(colour = "black", size = 12), axis.text.x = element_text(colour = "black", size = 12), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.015,0.015), limits = c(0, 0.65)) + scale_x_discrete(expand = c(0.01,0.01)) + geom_hline(yintercept = max(threshs, na.rm = TRUE), col = "#9d302a") + geom_smooth(method = "loess", span = 0.5, aes(group = 1), color = "#335b8d")

#measles in texas figure
fit$results$raw.perm.entropy[1:50] <- NA

quartz(width = 15, height = 8)
par(mar = c(5,4,4,4))
plot(fit$results$raw.perm.entropy, type = "l", cex = 2, xlim = c(1,2740), xaxt = "n", yaxt = "n", bty = "n", lwd = 3, col = "#4d4d4d", xlab = "", ylab = "")
par(new = TRUE)
plot(data[["Measles-Texas"]], type = "l", lty = 3, xlim = c(1,2800), xaxt = "n", yaxt = "n", bty = "n", lwd = 3, col = "#b2182b", xlab = "Year", ylab = "Permutation Entropy")
at.x <- seq(0, length(data[["Measles-Texas"]]), length.out = 5)
axis(1, at = at.x, labels = years[at.x+1])
at.y.z <- seq(min(data[["Measles-Texas"]], na.rm = TRUE), max(data[["Measles-Texas"]], na.rm = TRUE), length.out = 4)
axis(2, at = at.y.z, labels = round(seq(0.5, 1, length.out = 4), 2), las = 2)
axis(4, at = at.y.z, labels = floor(at.y.z))
mtext("Weekly cases", side=4, line=3)
legend(2200, 3600, legend = c("Permutation Entropy", "Measles-Texas"), col = c("#4d4d4d", "#b2182b"), lwd = 3, lty = c(1,3))
abline(v = min(which(years == 1963)), lwd = 2, col = "#2b8cbe")
text(min(which(years == 1963))+220, 6795, "Measles vaccine licensed")
