# require(zoo)

## Sebastian Sippel
# 04.08.2016
#' @title A function to compute bitflip statistics and time series
#' @export
#' @description Computation of bitflip statistics of a time series
#' @usage nbitflips(x, ndemb)
#' @param x A numeric vector (e.g. a time series), from which the ordinal pattern distribution is to be calculated  
#' @param ndemb Embedding dimension of the ordinal patterns (i.e. sliding window size) for which bitflips are to be calculated. Should be chosen such as length(x) >> ndemb
#' @details 
#' This function returns a histogram and time series of the number of bitflips occurring in the associated ordinal patterns. NA values are allowed, and any pattern that contains at least one NA value will be ignored. WARNING: Can be slow with very long time series (n > 10^7).
#' @return A list with two entries is returned.
#' @references Sippel, S., 2014. Evaluating the carbon dynamics of biogeochemical models using statistical complexity measures. Master Thesis, University of Bayreuth.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' nbitflips(x = x, ndemb = 6)
### Determine number of bitflips in a time series
nbitflips <- function(x, ndemb) {
  
  # get number of bitflips:
  bitflips = c(as.numeric(sign(diff(x)[1:(length(x)-2)]) != sign(diff(x)[2:(length(x)-1)])))
  
  # number of bitflips in ordinal patterns:
  bitflips_ts = c(zoo::rollapply(data=bitflips, width=ndemb-2, FUN=sum, na.rm=F), rep(NA, ndemb-1))
  # get histogram of bitflip time series:
  bins <- numeric(length=(ndemb-1))
  names(bins) <- as.character(0:(ndemb-2))
  
  for (i in 1:(length(bitflips)-ndemb+3)) {
    index <- sum(bitflips[i:(i+ndemb-3)])+1  # number of bitflips from position i to i+ndemb-3
    bins[index] <- bins[index] + 1
  }
  
  ret.list = list(bitflips_ts, bins)
  names(ret.list) = c("bitflips_ts", "bitflips_hist")
  
  return(ret.list)
}
