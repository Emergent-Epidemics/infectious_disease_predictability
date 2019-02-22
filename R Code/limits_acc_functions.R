#SV Scarpino
#Misc. functions for https://www.nature.com/articles/s41467-019-08616-0

#libraries
library(pdc)
library(ggplot2)
library(lme4)
library(xtable)
library(stats)
library(entropy)
library(philentropy)
library(statcomp)

###########
#acc funcs#
###########
course_grain <- function(x, grain.by){
  id <- 1:(length(x)/grain.by)
  id <- rep(id, each = grain.by)
  if(length(id) < length(x)){
    id <- c(id, rep(max(id)+1, length(x) - length(id)))
  }
  cg.x <- by(x, id, sum, na.rm = TRUE)
  cg.x <- as.numeric(cg.x)
  return(cg.x)
}

rel_ent <- function(x, m.min = 2, m.max = 6, t.min = 1, t.max = 3, do_mc_ent = FALSE){
  if(length(x) < 10){
    stop("x must be at least 10 time steps long")
  }
  test <- entropyHeuristic(X = x, m.min = m.min, m.max = m.max, t.min = t.min, t.max = t.max)
  ent <- min(test$entropy.values[,3], na.rm = TRUE)
  
  ent_check <- c()
  for(t in (length(x)-10):length(x)){
    test.t <- entropyHeuristic(X = x[1:t], m.min = m.min, m.max = m.max, t.min = t.min, t.max = t.max)
    ent.t <- min(test.t$entropy.values[,3], na.rm = TRUE)
    ent_check <- c(ent_check, ent.t)
  }
  
  if(is.na(mean(abs(diff(ent_check)), na.rm = TRUE)) == TRUE){
    ent <- NA
  }
  
  if(is.na(mean(abs(diff(ent_check)), na.rm = TRUE)) == FALSE){
    if(mean(abs(diff(ent_check)), na.rm = TRUE) > 1e-2){
      ent <- NA
    }
  }
  
  ents <- c()
  ms <- c()
  taus <- c()
  if(do_mc_ent == TRUE){
    for(i in 1:1000){
      test.i <- entropyHeuristic(X = sample(x, length(x), replace = TRUE), m.min = m.min, m.max = m.max, t.min = t.min, t.max = t.max)
      ent.i <- min(test.i$entropy.values[,3], na.rm = TRUE)
      m.i <- test.i$m
      tau.i <- test.i$t
      ents <- c(ents, ent.i)
      ms <- c(ms, m.i)
      taus <- c(taus, tau.i)
    }
    use.rel.ent <- which.max(ents)
    rel.ent <- ents[use.rel.ent]
    rel.ent.m <- ms[use.rel.ent]
    rel.ent.tau <- taus[use.rel.ent]
    if(is.finite(rel.ent) == FALSE){
      rel.ent <- NA
      rel.ent.m <- NA
      rel.ent.tau <- NA
    }
  }else{
    rel.ent <- NA
    rel.ent.m <- NA
    rel.ent.tau <- NA
  }
  m <- test$m
  if(length(m) == 0){
    m <- NA
  }
  
  tau <- test$t
  if(length(tau) == 0){
    tau <- NA
  }
  
  return(list("rel.ent" = rel.ent,"ent" = ent, "ent.wave" = NA, "n" = length(na.omit(x)), "d" = m, "tau" = tau, "rel.ent.d" = rel.ent.m, "rel.ent.tau" = rel.ent.tau))
}

filt_interp <- function(x){
  
  #strip leading and trailing NAs
  first.numb <- min(which(is.na(x) == FALSE))
  
  if(is.finite(first.numb) == FALSE){
    return(NA)
  }
  
  x <- x[first.numb:length(x)]
  last.numb <- max(which(is.na(x) == FALSE))
  x <- x[1:last.numb]
  
  #linear interp. 
  for(i in 1:length(x)){
    if(is.na(x[i]) == TRUE){
      #x[i] <- 0
      x[i] = mean(x = c(x[i-1], x[i+1]), na.rm = TRUE)
    }
  }
  
  #if there are still NAs, throw out the data
  if(length(which(is.na(x))) > 0){
    x <- NA
  }
  return(x)
}

load_data <- function(data, dat, name, col_start = 3, do_n_filt = FALSE, week_filt = numeric(0), year_filt = numeric(0)){
  for(i in col_start:ncol(dat)){
    if(do_n_filt == TRUE){
      filt.i <- which(dat[,i] > summary(dat[,i])[2])
    }else{
      if(length(week_filt) > 0){
        filt.i <- which(dat$WEEK >= week_filt[1] & dat$WEEK <= week_filt[2])
      }else{
        if(length(year_filt) > 0){
          filt.i <- which(dat$YEAR >= year_filt[1] & dat$YEAR <= year_filt[2])
        }else{
          filt.i <- 1:nrow(dat)
        }
      }
    }
    data[[paste0(name,colnames(dat[i]))]] <- dat[filt.i,i]
  }
  return(data)
}

#function for start to finish with window
full_lenth_pred_window <- function(data, start, end, window, d_1 = 2, d_2 = 5){
  begin <- start
  starts <- seq(begin, (end-window), by = window)
  ends <- starts + window
  ends[length(ends)] <- end
  
  pdc.perm.entropy <- rep(NA, length(starts))
  n <- pdc.perm.entropy
  d <- pdc.perm.entropy
  tau <- pdc.perm.entropy
  raw.perm.entropy <- pdc.perm.entropy
  weighted.perm.entropy <- pdc.perm.entropy
  boot.d <- pdc.perm.entropy
  boot.tau <- pdc.perm.entropy
  filt.perm.entropy <- pdc.perm.entropy
  n.filt <- pdc.perm.entropy
  boot.n <- pdc.perm.entropy
  dists <- list()
  
  results <- data.frame(raw.perm.entropy, pdc.perm.entropy, n, d, tau, weighted.perm.entropy, boot.d, boot.tau, filt.perm.entropy, n.filt, boot.n)
  
  for(i in 1:length(starts)){
    start.i <- starts[i]
    end.i <- ends[i]
    x.i <- data[start.i:end.i]
    
    if(length(x.i) < 10){
      next
    }
    
    ent.i <- rel_ent(x = x.i, m.min = 3, m.max = 3, t.min = 1, t.max = 1, do_mc_ent = FALSE)
    
    if(is.na(ent.i$d) == TRUE) next
    
    ent_w.i <- entropy(weighted_ordinal_pattern_distribution(x.i, ent.i$d))/log(factorial(ent.i$d))
    ent_raw.i <- entropy(codebook(x.i, m = ent.i$d, t = 1))/log(factorial(ent.i$d))
    
    perm_dist <- codebook(x.i, m = ent.i$d)
    dists[[i]] <- perm_dist
    
    results[i,"pdc.perm.entropy"] <- ent.i$ent
    results[i,"raw.perm.entropy"] <- ent_raw.i
    results[i,"weighted.perm.entropy"] <- ent_w.i
    results[i,"n"] <- ent.i$n
    results[i,"d"] <- ent.i$d
    results[i, "tau"] <- ent.i$tau
    results[i, "boot.d"] <- ent.i$rel.ent.d
    results[i, "boot.tau"] <- ent.i$rel.ent.tau
  }
  return(list("results" = results, "dists" = dists))
}

#start to finish, no window
full_lenth_pred <- function(data, start, end, d_1 = 2, d_2 = 5){
  starts <- rep(start, length(start:end))
  ends <- seq(start, end, by = 1)
  
  pdc.perm.entropy <- rep(NA, length(starts))
  n <- pdc.perm.entropy
  d <- pdc.perm.entropy
  tau <- pdc.perm.entropy
  raw.perm.entropy <- pdc.perm.entropy
  weighted.perm.entropy <- pdc.perm.entropy
  boot.d <- pdc.perm.entropy
  boot.tau <- pdc.perm.entropy
  filt.perm.entropy <- pdc.perm.entropy
  n.filt <- pdc.perm.entropy
  boot.n <- pdc.perm.entropy
  dists <- list()
  
  results <- data.frame(pdc.perm.entropy, n, d, tau, raw.perm.entropy, weighted.perm.entropy, boot.d, boot.tau, filt.perm.entropy, n.filt, boot.n)
  
  for(i in 1:length(starts)){
    start.i <- starts[i]
    end.i <- ends[i]
    x.i <- data[start.i:end.i]
    
    if(length(x.i) < 10){
      next
    }
    
    ent.i <- rel_ent(x = x.i, m.min = 3, m.max = 3, t.min = 1, t.max = 1, do_mc_ent = FALSE)
    
    if(is.na(ent.i$d) == TRUE) next
    
    ent_w.i <- entropy(weighted_ordinal_pattern_distribution(x.i, ent.i$d))/log(factorial(ent.i$d))
    ent_raw.i <- entropy(codebook(x.i, m = ent.i$d, t = 1))/log(factorial(ent.i$d))
    
    perm_dist <- codebook(x.i, m = ent.i$d)
    dists[[i]] <- perm_dist
    
    results[i,"pdc.perm.entropy"] <- ent.i$ent
    results[i,"raw.perm.entropy"] <- ent_raw.i
    results[i,"weighted.perm.entropy"] <- ent_w.i
    results[i,"n"] <- ent.i$n
    results[i,"d"] <- ent.i$d
    results[i, "tau"] <- ent.i$tau
    results[i, "boot.d"] <- ent.i$rel.ent.d
    results[i, "boot.tau"] <- ent.i$rel.ent.tau
  }
  return(list("results" = results, "dists" = dists))
}

#comparing distributions
comp_dist <- function(dist_1, dist_2, size){
  if(is(try(JSD(rbind(dist_1, dist_2))))[1] == "try-error"){
    return(list("dist" = NA, "p" = NA, "cutoff" = NA))
  }
  
  js.i <- JSD(rbind(dist_1, dist_2))
  draw.i <- rmultinom(n = 5000, size, dist_1)
  norm.i <- apply(draw.i, 2, function(x) return(x/sum(x)))
  diffs.i <- c()
  samps <- 1:ncol(norm.i)
  for(n in 1:1000){
    use.i1 <- sample(samps, 1)
    use.i2 <- sample(samps[-use.i1], 1)
    diff.i <- JSD(rbind(norm.i[,use.i1], norm.i[,use.i2]))
    diffs.i <- c(diffs.i, diff.i)
  }
  p.i <- length(which(diffs.i > js.i))/n
  cutoff <- diffs.i[order(diffs.i)][n*0.95]
  return(list("dist" = js.i, "p" = p.i, "cutoff" = cutoff))
}

filt_lead_trail_NA <- function(x){
  x.i <- as.numeric(x)
  
  #filtering leading/trailing NAs, which are common in the Tycho data
  first.numb <- min(which(is.na(x.i) == FALSE))
  last.numb <- max(which(is.na(x.i) == FALSE))
  if(is.finite(first.numb) == FALSE){
    return(list("x" = NA))
  }
  
  x.i <- x.i[first.numb:last.numb]
  
  if(length(x.i) < 10){
    return(list("x" = NA))
  }
  NAs <- c(1, which(is.na(x.i) == TRUE), length(x.i))
  if(length(NAs) == 1 | length(NAs) == length(x.i)){
    return(list("x" = NA))
  }
  
  if(length(x.i) < 10){
    return(list("x" = NA))
  }
  return(list("x" = x.i, "first" = first.numb, "last" = last.numb))
}


sub_sample_pe <- function(time_series, sub_sample = TRUE, N = 1000){
  ents <- c()
  ds <- c()
  taus <- c()
  rows <- c()
  starts <- c()
  Ns <- c()
  diff_ents <- c()
  diff_ents_ds <- c()
  weight_ents <- c()
  
  if(sub_sample == TRUE){
    N <- N
  }else{
    N <- 1
  }
  
  pb <- txtProgressBar(1, N, style=3)
  for(n in 1:N){
    x.i <- as.numeric(time_series)
    x.i <- x.i[-which(x.i < 1)]
    
    if(sub_sample == TRUE){
      start.i <- sample(1:(length(x.i)-11), 1)
      n.i <- sample(11:500, 1)
      if(start.i+n.i > length(x.i)){
        n.i <- length(x.i) - start.i
      }
      x.i <- x.i[start.i:(start.i+n.i)]
    }else{
      x.i <- x.i
      start.i <- 1
      n.i <- length(x.i)
    }
    
    if(length(x.i)<10){
      next
    }
    
    ent.i <- rel_ent(x = x.i, m.min = 2, m.max = 7, t.min = 1, t.max = 1)
    diff.ent.i <- rel_ent(x = diff(x.i), m.min = 2, m.max = 7, t.min = 1, t.max = 1)
    
    ent_w.i <- try(entropy(
      statcomp::weighted_ordinal_pattern_distribution(x.i, ent.i$d))/log(factorial(ent.i$d)), silent = TRUE)
    
    if(is(ent_w.i)[1] == "try-error"){
      ent_w.i <- NA
    } 
    
    d.i <- ent.i$d
    tau.i <- ent.i$tau
    pe.i <- ent.i$ent
    ents <- c(ents, pe.i)
    diff_ents <- c(diff_ents, diff.ent.i$ent)
    diff_ents_ds <- c(diff_ents_ds, diff.ent.i$d)
    ds <- c(ds, d.i)
    taus <- c(taus, tau.i)
    rows <- c(rows, 1)
    starts <- c(starts, start.i)
    Ns <- c(Ns, n.i)
    weight_ents <- c(weight_ents, ent_w.i)
    setTxtProgressBar(pb, n)
  }
  out <- data.frame(ents, diff_ents, diff_ents_ds, ds, taus, rows, starts, Ns, weight_ents)
  return(out)
}

start_end_ent <- function(time_series, start, end){
  ents_shift <- c()
  pb <- txtProgressBar(start, end, style=3)
  for(i in start:end){
    x.i <- as.numeric(time_series[c(start-10):i])
    rm0 <- which(x.i < 1)
    if(length(rm0) > 0){
      x.i <- x.i[-rm0]
    }
    
    if(length(x.i) < 10){
      ents_shift <- c(ents_shift, NA)
      next
    }
    
    ent.i <- rel_ent(x = x.i, m.min = 2, m.max = 7, t.min = 1, t.max = 1)
    
    ents_shift <- c(ents_shift, ent.i$ent)
    
    setTxtProgressBar(pb, i)
  }
  return(ents_shift)
}

ts_JS <- function(fit_window){
  JSs <- c()
  dists <- c()
  for(i in 1:(length(fit_window$dists)-1)){
    loop.var <- 1:length(fit_window$dists)
    if(length(fit_window$dists[[i]]) == 0 | sum(is.na(fit_window$dists[[i]])) != 0){
      JSs <- c(JSs, NA)
      dists <- c(dists, NA)
      next
    }
    for(j in loop.var[-i]){
      if(length(fit_window$dists[[j]]) == 0|sum(is.na(fit_window$dists[[j]])) != 0){
        JSs <- c(JSs, NA)
        dists <- c(dists, NA)
        next
      }
      
      JS.ij <- JSD(rbind(fit_window$dists[[j]], fit_window$dists[[i]]))
      dist.ij <- abs(i-j)
      JSs <- c(JSs, JS.ij)
      dists <- c(dists, dist.ij)
    }
  }
  return(list("JSs" = JSs, "time_dists" = dists))
}