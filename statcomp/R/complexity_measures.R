## Sebastian Sippel
# 19.01.2015

# This script is to write functions that calculate complexity measures as a function of
# ANY ordinal pattern distribution!
# Any complexity measure can be computed also with a weighted OPD!

# Fisher Information: second argument should be a vector of the pattern coding



# ----------------------------------------------------
# Permutation Entropy:
# ----------------------------------------------------
#' @title A function to compute the permutation entropy
#' @export
#' @description Computation of the permutation entropy of a time series based on its ordinal pattern distribution (see Bandt and Pompe 2002).
#' Permutation entropy is a global information measure, hence insensitive to the permutation ordering scheme.
#' @usage permutation_entropy(opd)
#' @param opd A numeric vector that details an ordinal pattern distribution.
#' @details
#' This function calculates the permutation entropy as described in Bandt and Pompe 2002.
#' @return The normalized permutation entropy as a numeric value in the range [0,1].
#' @references Bandt, C. and Pompe, B., 2002. Permutation entropy: a natural complexity measure for time series. Physical review letters, 88(17), p.174102.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' permutation_entropy(opd)
permutation_entropy = function(opd) {
  # maximum Shannon Entropy:
  ssmax  = log( length(opd) )
  # compute Permutation entropy and return:
  PE = shannon_entropy(opd) / ssmax
  return(PE)
}

# ----------------------------------------------------
# MPR complexity using Jensen-Shannon Divergence:
# ----------------------------------------------------
#' @title A function to compute the MPR-complexity
#' @export
#' @description The function computes the MPR complexity, i.e. a generalized (global) complexity measure based on the Jenson-Shannon divergence.
#' @usage MPR_complexity(opd)
#' @param opd A numeric vector that details an ordinal pattern distribution.
#' @details
#' Generalized complexity measures combine an information measure (i.e. entropy) with the distance of the distribution from the uniform distribution ("disequilibrium").
#' As a global measure, MPR-complexity is insensitive to the permutation coding scheme.
#' @return The normalized MPR complexity measure in the range [0, 1].
#' @references Martin, M.T., Plastino, A. and Rosso, O.A., 2006. Generalized statistical complexity measures: Geometrical and analytical properties. Physica A: Statistical Mechanics and its Applications, 369(2), pp.439-462.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' MPR_complexity(opd)
MPR_complexity = function(opd) {
  
  alpha = 0.5
  p00  = 1./length(opd)
  pp0    = length(opd)
  pp1    = length(opd) - 1
  aaa    = ( 1. - alpha ) / pp0
  bbb    = alpha + aaa
  aux1   = bbb * log( bbb )
  aux2   = pp1 * aaa * log( aaa )
  aux3   = ( 1. - alpha ) * log( pp0 )
  q03    =  -1. / ( aux1 + aux2 + aux3 )
  
  # convert to probabilities:
  opd.prob = opd/sum(opd)
  
  # Jensen-Shannon Divergence (Disequilibrium measure, e.g. in Martin, Plastino and Rosso, 2006):
  opd.prob2 = alpha * (opd.prob)  + (1 - alpha) * p00 # Shannon Entropie von S[(P+P_e)/2] f?r ein bestimmtes pattern (nicht ausummiert), siehe auch: Rosso et al (2007), Eq. 3
  
  # Shannon Entropy of divergence measure:
  S_p_pe = shannon_entropy(opd.prob2)
  
  # Shannon Entropy of opd:
  S_p = shannon_entropy(opd=opd.prob)
  
  # Shannon Entropy of Equilibrium distribution:
  S_pe  = log( length(opd) )
  
  # Shannon entropy of OPD:
  H_s = S_p/S_pe
  
  # compute MPR-complexity:
  C_js = q03 * (S_p_pe - 0.5 * S_p - 0.5 * S_pe) * H_s
  return(C_js)
}





# ----------------------------------------------------
# Make high-level functions for complexity measures:
# ----------------------------------------------------
# complexity is a (fast) high-level function to calculate various complexity measures:
# x is a time series
# ndemb is the embedding dimension
# x = rnorm(10000)
# ndemb = 4


#' @title A function to compute global information and complexity measures for time series
#' @export
#' @description This is a high-level function that calculates global complexity measures directly from a given time series or ordinal pattern distribution.
#' @usage global_complexity(x = NA, opd = NA, ndemb)
#' @param opd A numeric vector that details an ordinal pattern distribution in a user-specified permutation coding scheme.
#' @param x (OPTIONAL) If opd is not specified, a time series vector x must be specified
#' @param ndemb (OPTIONAL) If x is given, the embedding dimension (ndemb) is required.
#' @details
#' This function calculates the following global measures of complexity and information:
#' \itemize{
#'  \item Permutation Entropy (PE, cf. Bandt and Pompe, 2002)
#'  \item Permutation Statistical complexity (MPR complexity, cf. Martin, Plastino and Rosso, 2006)
#'  \item Number of "forbiden patterns" (cf. Amigo 2010)
#' }
#' @return A named vector containing the three global complexity measures.
#' @references
#' Bandt, C. and Pompe, B., 2002. Permutation entropy: a natural complexity measure for time series. Physical review letters, 88(17), p.174102.
#' Martin, M.T., Plastino, A. and Rosso, O.A., 2006. Generalized statistical complexity measures: Geometrical and analytical properties. Physica A: Statistical Mechanics and its Applications, 369(2), pp.439-462.
#' Amigo, J., 2010. Permutation complexity in dynamical systems: ordinal patterns, permutation entropy and all that. Springer Science & Business Media.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' global_complexity(x = x, ndemb = 6)
#' # or:
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' global_complexity(opd = opd, ndemb = 6)
global_complexity = function(x = NA, opd = NA, ndemb) {
  # if no ordinal pattern distribution is specified:
  if (is.na(opd))  opd = ordinal_pattern_distribution(x=x, ndemb = ndemb)
  
  # calculate global complexity measures:
  H_s = permutation_entropy(opd)
  C_js = MPR_complexity(opd)
  nforbid = nforbiden(opd)
  
  # put global complexity measures into one vector:
  global_complexity = c(H_s, C_js, nforbid)
  names(global_complexity) = c('PE', 'MPR_Cjs', 'nforbiden')
  
  return(global_complexity)
}



# ----------------------------------------------------
# Number of forbiden patterns:
# ----------------------------------------------------
# returns the number of forbiden patterns (in %) from any ordinal pattern distribution:
nforbiden = function(opd) {
  return(length(which(opd == 0)) / length(opd))
}


# ----------------------------------
# AUXILIARY FUNCTIONS:
# ----------------------------------
# compute non-normalized Shannon Entropy:
#' @keywords internal
shannon_entropy = function(opd) {
  opd.prob = opd / sum(opd)
  H_s = (-1) * sum(sapply(opd.prob, FUN=function(prob) if (prob >= 1.e-30) return(prob * log(prob)) else return(0)))
  return(H_s)
}

# compute fisher information for ordinal pattern distribution:
#' @title A (low-level) function to compute the Fisher-information
#' @export
#' @description The function computes the Fisher information, i.e. a local information measure based on two different discretizations.
#' @usage fis(opd, discretization)
#' @param opd A numeric vector that details an ordinal pattern distribution in a user-specified permutation coding scheme.
#' @param discretization The discretization scheme to use, either 'Olivares.2012' or 'Ferri.2009'
#' @details
#' The Fisher information is a local information and complexity measure, computed based on the ordinal pattern distribution.
#' The Fisher information is based on local gradients, hence it is sensitive to the permutation coding scheme.
#' Options for discretization: 'Olivares.2012' or 'Ferri.2009', following Fisher Information discretization schemes in the respective publications.
#' @return The normalized Fisher information measure in the range [0, 1].
#' @references Olivares, F., Plastino, A. and Rosso, O.A., 2012. Ambiguities in Bandt-Pompe's methodology for local entropic quantifiers. Physica A: Statistical Mechanics and its Applications, 391(8), pp.2518-2526.
#' Ferri, G.L., Pennini, F. and Plastino, A., 2009. LMC-complexity and various chaotic regimes. Physics Letters A, 373(26), pp.2210-2214.
#' @author Sebastian Sippel
#' @examples
#' x = arima.sim(model=list(ar = 0.3), n = 10^4)
#' opd = ordinal_pattern_distribution(x = x, ndemb = 6)
#' fis(opd = opd)
fis = function(opd, discretization = 'Olivares.2012') {
  # convert to probabilities:
  opd.prob = opd / sum(opd)
  
  if (discretization == 'Olivares.2012') {
    ## calculate Fisher information according to Olivares et al 2012
    fis = sum(sapply(X=1:(length(opd.prob)-1), FUN=function(i) return((sqrt(opd.prob[i+1]) - sqrt(opd.prob[i]))^2))) / 2
  } else if (discretization == 'Ferri.2009') {
    fis = 0.5 * sum(sapply(X=1:(length(opd.prob)-1), FUN = function(i) (opd.prob[i+1] - opd.prob[i]) ^ 2 / (opd.prob[i+1] + opd.prob[i])), na.rm=T)
  }
  
  return(fis)
}

# test = rnorm(1000)
# fis(opd=ordinal_pattern_distribution(x=test, ndemb=4))
# fis(opd=ordinal_pattern_distribution(x=test, ndemb=4), discretization='Ferri.2009')
