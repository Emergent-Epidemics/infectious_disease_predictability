## Sebastian Sippel
# 05.01.2014


#' @title A function to compute ordinal pattern statistics
#' @export
#' @useDynLib statcomp
#' @description Computation of the ordinal patterns of a time series (see e.g. Bandt and Pompe 2002)
#' @usage ordinal_pattern_distribution(x, ndemb)
#' @param x A numeric vector (e.g. a time series), from which the ordinal pattern distribution is to be calculated  
#' @param ndemb Embedding dimension of the ordinal patterns (i.e. sliding window size). Should be chosen such as length(x) >> ndemb
#' @details 
#' This function returns the distribution of ordinal patterns using the Keller coding scheme, detailed in Physica A 356 (2005) 114-120. NA values are allowed, and any pattern that contains at least one NA value will be ignored.
#' (Fast) C routines are used for computing ordinal patterns.
#' @return A character vector of length factorial(ndemb) is returned.
#' @references Bandt, C. and Pompe, B., 2002. Permutation entropy: a natural complexity measure for time series. Physical review letters, 88(17), p.174102.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' ordinal_pattern_distribution(x = x, ndemb = 6)
ordinal_pattern_distribution = function(x, ndemb) {
  
  
  epsilon=1.e-10
  npdim=factorial(ndemb)
  
  #Berechnungs der Ordnungsstatistik nach der Kodierung von Karsten Keller:
  #Physica A 356 (2005) 114-120
  
  nfac=factorial(ndemb)
  
  ifrec=.C("ordinal_pattern_loop",
           as.double(x),
           as.integer(length(x)),
           as.integer(ndemb),
           integer(nfac),
           as.integer(nfac),as.integer(rep(0,length(x))),NAOK=T)[[4]]
  
  # ifrec is the ordinal pattern distribution in the Keller coding scheme!
  return(ifrec)
}


#' @title A function to compute time series of ordinal patterns
#' @export
#' @description Computation of the ordinal patterns of a time series (see e.g. Bandt and Pompe 2002)
#' @usage ordinal_pattern_time_series(x, ndemb)
#' @param x A numeric vector (e.g. a time series), from which the ordinal pattern time series is to be calculated  
#' @param ndemb Embedding dimension of the ordinal patterns (i.e. sliding window size). Should be chosen such as length(x) >> ndemb
#' @details 
#' This function returns the distribution of ordinal patterns using the Keller coding scheme, detailed in Physica A 356 (2005) 114-120. NA values are allowed, and any pattern that contains at least one NA value will be ignored.
#' (Fast) C routines are used for computing ordinal patterns.
#' @return A character vector of length(x) is returned.
#' @references Bandt, C. and Pompe, B., 2002. Permutation entropy: a natural complexity measure for time series. Physical review letters, 88(17), p.174102.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' ordinal_pattern_time_series(x = x, ndemb = 6)
ordinal_pattern_time_series = function(x, ndemb) {
  
  
  epsilon=1.e-10
  npdim=factorial(ndemb)
  
  #Berechnungs der Ordnungsstatistik nach der Kodierung von Karsten Keller:
  #Physica A 356 (2005) 114-120
  
  nfac=factorial(ndemb)
  
  ifrec_ts=.C("ordinal_pattern_loop",
           as.double(x),
           as.integer(length(x)),
           as.integer(ndemb),
           integer(nfac),
           as.integer(nfac),as.integer(rep(0,length(x))),NAOK=T)[[6]]
  ifrec_ts[(length(ifrec_ts) - (ndemb-2)):length(ifrec_ts)] <- NA
  
  # ifrec is the ordinal pattern distribution in the Keller coding scheme!
  return(ifrec_ts)
}


## Sebastian Sippel
# 05.01.2014
#' @keywords internal
ordinal_pattern_distribution_2 = function(x, ndemb) {
  
  ### Deal with gaps in the sliding window time series:
  # get indices to run through for calculation of complexity measures
  gapfree = stats::na.omit(sapply(1:(length(x)-ndemb + 1), FUN = function(y) if(!any(is.na(x[y:(y+ndemb-1)]))) return(y) else return(NA)))
  
  epsilon=1.e-10
  npdim=factorial(ndemb)
  
  #Berechnungs der Ordnungsstatistik nach der Kodierung von Karsten Keller:
  #Physica A 356 (2005) 114-120
  
  ifrec = numeric(length=npdim)  #ersetzt die for-schleife zum erstellen des Vektors, #for ip=1:npdim; ifrec( ip ) = 0; end;
  
  ## introduce for loop:
  for (nv in 1:(length(gapfree))) {
    
    xvec <- x[gapfree[nv]:(gapfree[nv] + ndemb - 1)]
    
    ## only if no gaps are in the "word" of the time series:
    ipa = matrix(data=0,nrow=ndemb, ncol=ndemb)  #Inversionsmatrix
    
    for (il in 2:ndemb) {
      for (it in il:ndemb) { 
        ipa[it, il] = ipa[it-1, il-1]
        if( (xvec[it] <= xvec[it - ( il - 1 ) ] ) || ( abs( xvec[it - ( il - 1)] - xvec[it]) < epsilon))
          ipa[ it, il ] = ipa[ it, il ] + 1;
      }
    }
    
    nd = ipa[ndemb,2] 
    for (il in 3:ndemb) {
      nd =il * nd + ipa[ndemb, il]
    }
    
    ifrec[nd + 1] = ifrec[nd+ 1] + 1;       
  }
  
  # ifrec is the ordinal pattern distribution in the Keller coding scheme!
  return(ifrec)
}