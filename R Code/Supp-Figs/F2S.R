#SV Scarpino
#F2S for https://arxiv.org/abs/1703.07317

#set working dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to source file location

#libraries (not included in limits_acc_functions.R)
library(EpiModel)
library(statcomp)

###########
#acc funcs#
###########
source("../limits_acc_functions.R")

###############
#Global Params#
###############
n_iter <- 100
run_new <- TRUE #set to TRUE to re-run the iterations

######
#data#
######

#Zika Colombia
z_col_files <- list.files("../../Data/Colombia/Epidemiological_Bulletin/data", full.names = TRUE)
z_col_files <- z_col_files[-c(1:2)] #different locations

places <- c()
for(i in z_col_files){
  dat.i <- read.csv(i)
  places <- c(places, as.character(unique(dat.i$location)))
}

zika_col <- matrix(NA, nrow = length(z_col_files), ncol = length(unique(places))+1)
zika_col <- as.data.frame(zika_col)
names <- unique(places)
for(i in 1:length(z_col_files)){
  z.col.i <- read.csv(z_col_files[i])
  z.c.i <- by(z.col.i$value, z.col.i$location, sum, na.rm = TRUE)
  mt <- match(names, names(z.c.i))
  zika_col[i,-1] <- as.numeric(z.c.i)[mt]
  
  zika_col[i,1] <- as.character(unique(z.col.i$report_date))
}
dates_zc <- strptime(zika_col[,1], format = "%Y-%m-%d")
ord_zc <- order(dates_zc)
zika_col <- zika_col[ord_zc,]
use <- which(strptime(zika_col[,1], format = "%Y-%m-%d") < strptime("2017-01-01", format = "%Y-%m-%d"))
zika_col <- zika_col[use, -1]
colnames(zika_col) <- iconv(names, "latin1", "ASCII", sub="")

#influenza
flu_tycho <- read.csv("../../Data/INFLUENZA_Cases_1919-1951_20151014132002.csv", header = TRUE, skip = 2, na.strings = "-", stringsAsFactors = FALSE)

#dengue
den <- read.csv("../../Data/San_Juan_Training_Data.csv", header = TRUE, skip = 0, stringsAsFactors = FALSE)
den <- data.frame(den[,"week_start_date"], den[,"total_cases"])
colnames(den) <- c("date", "count")
den[,"date"] <- as.Date(den[,"date"], format = "%Y-%m-%d")

#polio
polio <- read.csv("../../Data/POLIOMYELITIS_Cases_1921-1971_20150923114821.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#measles
meas <- read.csv("../../Data/MEASLES_Cases_1909-2001_20150923120449.csv", na.strings = "NA", header = TRUE, stringsAsFactors = FALSE)

#whooping cough
whoop <- read.csv("../../Data/WHOOPING_COUGH_[PERTUSSIS]_Cases_1909-2014_20160607104756.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#hepatitis A
hepa <- read.csv("../../Data/HEPATITIS_A_Cases_1966-2014_20160707103116.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#mumps
mump <- read.csv("../../Data/MUMPS_Cases_1967-2014_20160707103045.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#Gonorrhea
gono <- read.csv("../../Data/GONORRHEA_Cases_1972-2014_20160707103202.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#Chlamydia
chlam <- read.csv("../../Data/CHLAMYDIA_Cases_2006-2014_20160707103149.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

static <- read.csv("../../Data/exp_static_simulation_data_bulk_100.csv", header = TRUE, stringsAsFactors = FALSE)

##########
#analysis#
##########

#building data set
data <- list("Dengue-Puerto Rico" = den$count)

for(i in 1:ncol(zika_col)){
  names_i <- gsub(pattern = "-", "", colnames(zika_col)[i])
  data[[paste0("Zika-", names_i)]] <- zika_col[,i]
}


data_sets <- list(flu_tycho, polio, meas, whoop, mump, gono, hepa, chlam)
names <- c("Influenza-", "Polio-", "Measles-", "Whooping Cough-", "Mumps-", "Gonorrhea-", "Hepatitis A-", "Chlamydia-")

for(d in 1:length(data_sets)){
  data <- load_data(data = data, dat = data_sets[[d]], name = names[d], col_start = 3, do_n_filt = FALSE)
}

for(i in 2:ncol(static)){
  data[[paste0("Sim-",colnames(static)[i])]] <- static[,i]
}

data[["noise-full"]] <- rnorm(1000)
data[["noiseMissing-missing"]] <- data[["noise-full"]]
data[["noiseMissing-missing"]][sample(1:1000, 100, replace = FALSE)] <- NA


#checking NAs
NAs <- rep(NA, length(data))
Ns <- rep(NA, length(data))
for(i in 1:length(data)){
  NAs[i] <- length(which(is.na(as.numeric(data[[i]])) == TRUE))
  Ns[i] <- length(as.numeric(data[[i]]))
}

if(run_new == TRUE){
  names <- unlist(strsplit(names(data), "-"))
  disease <- names[seq(1, length(names), by = 2)]
  location <- names[seq(2, length(names), by = 2)]
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
  
  results <- data.frame(disease, location, pdc.perm.entropy, n, d, tau, boot.d, boot.tau, filt.perm.entropy, n.filt, boot.n, weighted.perm.entropy, total.cases)
  
  #looping
  pb <- txtProgressBar(1, length(data), style=3)
  for(iter in 1:n_iter){
    for(i in 1:length(data)){
      # if(length(grep("influenza", names(data)[i], ignore.case = T))==0){
      #   next
      # }

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
      length.i <- sample(11:500, 1)
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
    results$location <- as.character(results$location)
    if(iter == 1){
      RESULTS <- results
    }else{
      RESULTS <- rbind(RESULTS, results)
    }
    setTxtProgressBar(pb, iter)
  }
  
  time_stamp <- as.numeric(Sys.time())
  save(RESULTS, file = paste0(time_stamp, "_big_results_PE.RData"))
}else{
  load("../../Results/1529387566.44103_big_results_PE.RData")
}

diseases <- unique(RESULTS$disease)
for(i in 1:length(diseases)){
  use.i <- which(RESULTS$disease == diseases[i])
  y.i <- RESULTS$weighted.perm.entropy[use.i]
  x.i <- RESULTS$n[use.i]
  loc.i <- RESULTS$location[use.i]
  
  rm.i <- which(is.finite(y.i) == FALSE)
  if(length(rm.i) > 0){
    y.i <- y.i[-rm.i]
    x.i <- x.i[-rm.i]
    loc.i <- loc.i[-rm.i]
  }
  
  if(length(unique(loc.i)) > 1){
    ent.i <- by(data = y.i, INDICES = x.i, FUN = mean, na.rm = TRUE)
    ent.low <- by(data = y.i, INDICES = x.i, FUN = function(x) return(summary(x)[2]))
    ent.high <- by(data = y.i, INDICES = x.i, FUN = function(x) return(summary(x)[5]))
    n.i <- as.numeric(names(ent.i))
  }else{
    ent.i <- by(data = y.i, INDICES = x.i, FUN = mean, na.rm = TRUE)
    ent.low <- by(data = y.i, INDICES = x.i, FUN = function(x) return(summary(x)[2]))
    ent.high <- by(data = y.i, INDICES = x.i, FUN = function(x) return(summary(x)[5]))
    
    n.i <- as.numeric(names(ent.i))
  }
  
  
  disease.i <- rep(diseases[i], length(ent.i))
  res.i <- data.frame(disease.i, 1-as.numeric(ent.i), n.i, 1-as.numeric(ent.low), 1-as.numeric(ent.high))
  colnames(res.i) <- c("disease", "PE", "n", "PEmin", "PEmax")
  
  if(i == 1){
    results <- res.i
  }else{
    results <- rbind(results, res.i)
  }
} 

#manuscript figures
R0s <- c("Dengue" = 4, "Influenza" = 1.5, "Measles" = 15, "Polio" = 6, "Whooping Cough" = 5.5, "Mumps" = 5.5, "Gonorrhea" = 1, "Hepatitis A" = 1.2, "Chlamydia" = 1)

####
#F1#
####
use.F1 <- which(results$disease %in% c("Chlamydia", "Gonorrhea", "Hepatitis A", "Influenza", "Measles", "Mumps", "Polio", "Whooping Cough", "Dengue", "Sim") & results$n > 4)

plot.dat <- results[use.F1, ]
plot.dat$disease <- as.character(plot.dat$disease)
plot.dat$disease  <- as.factor(plot.dat$disease)

pal <- c("#e7298a", "gray", "#7570b3","#666666", "#b2182b", "#2166ac" , "#1b9e77", "#E6AB02", "black","#ff7f00")

p0 <- ggplot(data = plot.dat, aes(x = n, y = PE, colour = disease))

quartz(width = 7, height = 6)

p0 + geom_line(size = 1) + scale_colour_manual(values = pal, guide_legend(title = "Disease")) + xlab("Number of weeks") + ylab("Predictability (1 - H)") + theme(legend.position = c(0.86, 0.73), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "#f0f0f0", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01), limits = c(0.0,1)) + scale_x_continuous(expand = c(0.01,0.01), limits = c(0, 505))+geom_ribbon(aes(x = n, ymin = PEmin, ymax = PEmax, fill = disease), alpha = 0.25, colour = "#00000000") + scale_fill_manual(values = pal, guide_legend(title = "Disease")) + geom_line(size = 1)
