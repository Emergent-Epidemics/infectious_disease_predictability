#SV Scarpino
#F3 for https://arxiv.org/abs/1703.07317

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

######
#data#
######

#Zika Colombia
z_col_files <- list.files("../Data/Colombia/Epidemiological_Bulletin/data", full.names = TRUE)
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
flu_tycho <- read.csv("../Data/INFLUENZA_Cases_1919-1951_20151014132002.csv", header = TRUE, skip = 2, na.strings = "-", stringsAsFactors = FALSE)

#dengue
den <- read.csv("../Data/San_Juan_Training_Data.csv", header = TRUE, skip = 0, stringsAsFactors = FALSE)
den <- data.frame(den[,"week_start_date"], den[,"total_cases"])
colnames(den) <- c("date", "count")
den[,"date"] <- as.Date(den[,"date"], format = "%Y-%m-%d")

#polio
polio <- read.csv("../Data/POLIOMYELITIS_Cases_1921-1971_20150923114821.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#measles
meas <- read.csv("../Data/MEASLES_Cases_1909-2001_20150923120449.csv", na.strings = "NA", header = TRUE, stringsAsFactors = FALSE)

#whooping cough
whoop <- read.csv("../Data/WHOOPING_COUGH_[PERTUSSIS]_Cases_1909-2014_20160607104756.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#hepatitis A
hepa <- read.csv("../Data/HEPATITIS_A_Cases_1966-2014_20160707103116.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#mumps
mump <- read.csv("../Data/MUMPS_Cases_1967-2014_20160707103045.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#Gonorrhea
gono <- read.csv("../Data/GONORRHEA_Cases_1972-2014_20160707103202.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

#Chlamydia
chlam <- read.csv("../Data/CHLAMYDIA_Cases_2006-2014_20160707103149.csv", skip = 2, na.strings = "-", header = TRUE, stringsAsFactors = FALSE)

static <- read.csv("../Data/exp_static_simulation_data_bulk_100.csv", header = TRUE, stringsAsFactors = FALSE)

#building data set
data <- list("Dengue-Puerto Rico" = den$count)

for(i in 1:ncol(zika_col)){
  names_i <- gsub(pattern = "-", "", colnames(zika_col)[i])
  data[[paste0("Zika-", names_i)]] <- zika_col[,i]
}

for(i in 2:ncol(static)){
  data[[paste0("Sim-",colnames(static)[i])]] <- static[,i]
}

data[["noise-full"]] <- rnorm(1000)
data[["noiseMissing-missing"]] <- data[["noise-full"]]
data[["noiseMissing-missing"]][sample(1:1000, 100, replace = FALSE)] <- NA

#uncomments each year_range makes things pre-vaccine
year_range <- c(1900,2017)
for(i in 3:ncol(flu_tycho)){
  filt.i <- which(flu_tycho$YEAR >= year_range[1] & flu_tycho$YEAR <= year_range[2])
  filt.sea.i <- which(flu_tycho$WEEK[filt.i]>=40|flu_tycho$WEEK[filt.i]<=20)
  data[[paste0("Influenza-",colnames(flu_tycho[i]))]] <- flu_tycho[filt.i[filt.sea.i],i]
}

#year_range <- c(1930,1950)
for(i in 3:ncol(polio)){
  filt.i <- which(polio$YEAR >= year_range[1] & polio$YEAR <= year_range[2])
  data[[paste0("Polio-",colnames(polio[i]))]] <- polio[filt.i,i]
}

#year_range <- c(1930,1965)
for(i in 3:ncol(meas)){
  filt.i <- which(meas$YEAR >= year_range[1] & meas$YEAR <= year_range[2])
  data[[paste0("Measles-",colnames(meas[i]))]] <- meas[filt.i,i]
}

#year_range <- c(1930,1955)
for(i in 3:ncol(whoop)){
  filt.i <- which(whoop$YEAR >= year_range[1] & whoop$YEAR <= year_range[2])
  data[[paste0("Whooping Cough-",colnames(whoop[i]))]] <- whoop[filt.i,i]
}

#year_range <- c(1930,1975)
for(i in 3:ncol(mump)){
  filt.i <- which(mump$YEAR >= year_range[1] & mump$YEAR <= year_range[2])
  data[[paste0("Mumps-",colnames(mump[i]))]] <- mump[filt.i,i]
}

#year_range <- c(2005,2017)
for(i in 3:ncol(gono)){
  filt.i <- which(gono$YEAR >= year_range[1] & gono$YEAR <= year_range[2])
  data[[paste0("Gonorrhea-",colnames(gono[i]))]] <- gono[filt.i,i]
}

#year_range <- c(1960,1990)
for(i in 3:ncol(hepa)){
  filt.i <- which(hepa$YEAR >= year_range[1] & hepa$YEAR <= year_range[2])
  data[[paste0("Hepatitis A-",colnames(hepa[i]))]] <- hepa[filt.i,i]
}

#year_range <- c(2006,2017)
for(i in 3:ncol(chlam)){
  filt.i <- which(chlam$YEAR >= year_range[1] & chlam$YEAR <= year_range[2])
  data[[paste0("Chlamydia-",colnames(hepa[i]))]] <- chlam[filt.i,i]
}

##########
#Analysis#
##########
NAs <- rep(NA, length(data))
Ns <- rep(NA, length(data))
for(i in 1:length(data)){
  NAs[i] <- length(which(is.na(as.numeric(data[[i]])) == TRUE))
  Ns[i] <- length(as.numeric(data[[i]]))
}

names <- unlist(strsplit(names(data), "-"))
disease <- names[seq(1, length(names), by = 2)]
location <- names[seq(2, length(names), by = 2)]
boot.perm.entropy <- rep(NA, length(data))
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

results <- data.frame(disease, location, boot.perm.entropy, n, d, tau, boot.d, boot.tau, filt.perm.entropy, n.filt, boot.n, weighted.perm.entropy)

#looping
pb <- txtProgressBar(1, length(data), style=3)
for(i in 1:length(data)){
  x.i <- filt_lead_trail_NA(as.numeric(data[[i]]))$x
  
  if(is.na(x.i) == TRUE){
    next
  }
  
  ent.i <- rel_ent(x = x.i, m.min = 2, m.max = 7, t.min = 1, t.max = 1, do_mc_ent = FALSE)
  
  if(is.na(ent.i$d) == TRUE) next
  ent_w.i <- entropy(weighted_ordinal_pattern_distribution(x.i, ent.i$d))/log(factorial(ent.i$d))
  
  results[i,"boot.perm.entropy"] <- ent.i$rel.ent
  results[i,"raw.perm.entropy"] <- ent.i$ent
  results[i,"weighted.perm.entropy"] <- ent_w.i
  results[i,"n"] <- ent.i$n
  results[i,"d"] <- ent.i$d
  results[i, "tau"] <- ent.i$tau
  results[i, "boot.d"] <- ent.i$rel.ent.d
  results[i, "boot.tau"] <- ent.i$rel.ent.tau
  
  setTxtProgressBar(pb, i)
}

results$disease <- as.character(results$disease)
results$location <- as.character(results$location)

######
#Plot#
######
R0s <- c("Influenza" = 1.47, "Measles" = 15.1, "Polio" = 5.36, "Whooping Cough" = 14.75, "Mumps" = 9.94, "Gonorrhea" = 1.34, "Hepatitis A" = 2.45, "Chlamydia" = 0.99, "Zika" = 2.7)

####
#F1#
####
use.F1 <- which(results$disease %in% c("Chlamydia", "Gonorrhea", "Hepatitis A", "Influenza", "Measles", "Mumps", "Polio", "Whooping Cough", "Dengue") & ! results$location %in% c("GUAM", "VIRGIN.ISLANDS", "PAC.TRUST.TERR", "NORTHERN.MARIANA.ISLANDS", "X", "AMERICAN.SAMOA"))

plot.dat <- results[use.F1, ]
rm.Inf <- which(plot.dat$raw.perm.entropy > 1)
if(length(rm.Inf) > 0){
  plot.dat <- plot.dat[-rm.Inf,]
}
plot.dat$d_trans <- plot.dat$n * 1/R0s[plot.dat$disease]
m.plot <- lmer(log(raw.perm.entropy)~log(d_trans)+(log(d_trans)|disease), data = plot.dat)

pal <- c("#e7298a", "gray", "#7570b3","#666666", "#b2182b", "#2166ac" , "#1b9e77", "#E6AB02", "#ff7f00")

p0 <- ggplot(data = plot.dat, aes(x = n, y = 1-raw.perm.entropy, colour = disease))
p1 <- ggplot(data = plot.dat, aes(x = log(n), y = log(raw.perm.entropy), colour = disease))
  
quartz(width = 7, height = 6)

p0 + geom_point() + scale_colour_manual(values = pal, guide_legend(title = "Disease")) + xlab("Time series length (weeks)") + ylab("Predictability") + theme(legend.position = c(0.84, 0.74), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01), limits = c(0,1)) + scale_x_continuous(expand = c(0.01,100))

quartz(width = 7, height = 6)
p1 + geom_point() + scale_colour_manual(values = pal, guide_legend(title = "Disease")) + xlab("Weeks/R0 (log)") + ylab("Permutation entropy (log)") + theme(legend.position = c(0.84, 0.28), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "#f0f0f0", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01)) + scale_x_continuous(expand = c(0.01,0.01))+geom_abline(aes(intercept = fixef(m.plot)[1], slope = fixef(m.plot)[2]))
