# This script implements distance measures for Ordinal-pattern distributions:
# Sebastian Sippel
# 21.05.2015


### Hellinger distance
#% Calculates the Hellinger distance for two discrete distributions
#% given through the two vectors p and q of equal length.
#%p(i) >= 0, q(i) >= 0 for all i.
#%
#% Holger Lange 25.10.2013
# in the range 0...1

#' @title Distance measure between ordinal pattern distributions: Hellinger distance
#' @export
#' @description Compute the Hellinger Distance
#' @usage hellinger_distance(p, q)
#' @param p An ordinal pattern distribution
#' @param q A second ordinal pattern distribution to compare against p.
#' @details 
#' This function returns a distance measure.
#' @return A vector of length 1.
#' @references none
#' @author Sebastian Sippel
#' @examples
#' p = ordinal_pattern_distribution(rnorm(10000), ndemb = 5)
#' q = ordinal_pattern_distribution(arima.sim(model=list(ar=0.9), n= 10000), ndemb = 5)
#' hellinger_distance(p=p, q = q)
hellinger_distance = function(p,q) {
  h=0;
  if (length(which(p<0)) != 0 | length(which(q<0)) != 0) {
    print('positive values required'); return() }
  
  if (length(p) != length(q)) {
    print('Distributions must have the same length!') }
  
  #Normalization:
  ht=0;
  p = p/sum(p);
  q = q/sum(q);
  
  for (i in 1:length(p)) {
    ht=ht+(sqrt(p[i]) - sqrt(q[i]))^2
  }
  
  h=sqrt(ht/2);
  return(h)
}


#### Jensen-Shannon Distance:
### Jensen Shannon Divergence 
# normalized according to MPR (2006), q must be a uniform PDF!
# Sebastian Sippel
# Martin, M. T., A. Plastino, and O. A. Rosso. "Generalized statistical complexity measures: Geometrical and analytical properties." Physica A: Statistical Mechanics and its Applications 369.2 (2006): 439-462.

#' @title Generalized disequilibrium measure for ordinal pattern distributions based on the Jensen-Shannon Divergence
#' @export
#' @description Computes a normalized form of the Jensen-Shannon Divergence
#' @usage jensen_shannon_divergence(p, q="unif")
#' @param p An ordinal pattern distribution
#' @param q A second ordinal pattern distribution to compare against p, or a character vector q="unif" (comparison of p to uniform distribution)
#' @details 
#' This function returns a distance measure.
#' @return A vector of length 1.
#' @references Martin, M.T., Plastino, A. and Rosso, O.A., 2006. Generalized statistical complexity measures: Geometrical and analytical properties. Physica A: Statistical Mechanics and its Applications, 369(2), pp.439-462.
#' @author Sebastian Sippel
#' @examples
#' p = ordinal_pattern_distribution(rnorm(10000), ndemb = 5)
#' q = ordinal_pattern_distribution(arima.sim(model=list(ar=0.9), n= 10000), ndemb = 5)
#' jensen_shannon_divergence(p = p, q = q)
jensen_shannon_divergence = function(p,q="unif") {
  
  if ( any(is.na(p)) | any(is.na(q))) return(NA)
  
  j_s =0
  ## Assign uniform distribution to vector q if not specified otherwise
  if (q=="unif") {q <- rep(1, times=length(p))
                  print("Second vector not specified: Use Uniform Distribution in q!")}
  
  if (length(which(p<0)) != 0 | length(which(q<0)) != 0) {
    print('positive values required'); return() }
  
  if (length(p) != length(q)) {
    print('Distributions must have the same length!') }
  
  for (i in 2:length(q)) {
    if (q[i-1] != q[i])
      print('WARNING: q is not a uniform distribution! Generalized Disequilibrium might not be in the range [0,1]')
  } 
  
  alpha = 0.5
  pp0    = length(p)  #npdim
  pp1    = length(p) - 1  #npdim - 1 
  p00    = 1./pp0 
  aaa    = ( 1. - alpha ) / pp0
  bbb    = alpha + aaa
  aux1   = bbb * log( bbb )
  aux2   = pp1 * aaa * log( aaa )
  aux3   = ( 1. - alpha ) * log( pp0 )
  q03    =  -1. / ( aux1 + aux2 + aux3 )
  ssmax  = log( pp0 )
  
  #Normalization:
  p = p/sum(p);
  q = q/sum(q);
  
  for (i in 1:length(p)) {
    if(p[i] >= 1.e-30)
      S_p <- p[i] * log(p[i])
    else
      S_p = 0.  
    if(q[i] >= 1.e-30)
      S_q <- q[i] * log(q[i])
    else
      S_q = 0.  
    if(p[i] + q[i] >= 1.e-30)
      S_pq <- (p[i] + q[i])/2 * log((p[i] + q[i])/2)
    else
      S_pq = 0.  
    
    j_s <- j_s + S_p/2 + S_q/2 - S_pq
  }
  return(q03*j_s)
}



# -------------------------------------------------------------
# Distance measures within the Entropy-Complexity plane:
# -------------------------------------------------------------
# will follow



