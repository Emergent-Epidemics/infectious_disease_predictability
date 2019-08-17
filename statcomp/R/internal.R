# -------------------------------------------------------
# AUXILIARY FUNCTIONS
# -------------------------------------------------------
# This is a collection of package-internal, auxiliary functions:


# var.fun: calculates the variance, such as in Fadlallah et al (2013)
#' @keywords internal
var.fun = function(x) stats::var(x) * (length(x) - 1) / length(x)


# ------------------------------------------------------
# Internal functions to compute pattern coding schemes:
# all permutation codings are defined as ranks (i.e. in contrast to e.g. Olivares et al 2012)
# for comparability to Olivares et al, those matrices can still be generated


#' @keywords internal
Olivares_sorting <- function(x, ndemb) {
  
  # sorts based on the Olivares et al (2012) paper:
  pi = seq(0, ndemb-1, 1)
  epsilon = 10^-8
  
  for (l in 1:(ndemb-1)) {
    for (k in 1:(ndemb - l)) {
      if ( x[k] > x[k+1] | (abs(x[k] - x[k+1] ) <= epsilon)) {
        auxp = x[k]
        auxpi = pi[k]
        x[k] = x[k+1]
        pi[k] = pi[k+1]
        x[k+1] = auxp
        pi[k+1] = auxpi      
      }    
    }
  }
  return(pi)
}



#' @keywords internal
njumps <- function(permutations_matrix) {
  D = dim(permutations_matrix)[2]
  njumps <- numeric(length=factorial(D))
  
  for (i in 1:factorial(D)) {
    for (j in 1:(D-1)) {
      njumps[i] = njumps[i] + abs(permutations_matrix[i,j] - permutations_matrix[i,j+1]) -1
    }
  }
  return(njumps)
}

#' @keywords internal
findinversions <- function(patterns) {
  d <- dim(patterns)[2]
  inver <- NULL
  patdiff <- NULL
  for (i in 1:dim(patterns)[1]) {
    inver[i] = 0
    for (k in (1:(d-1))) {
      patdiff[k] = sign(patterns[i,k+1]-patterns[i,k])
    }
    inver[i] = 0
    for (l in 1:(length(patdiff)-1)) {
      inver[i] = inver[i] + abs(sign(patdiff[l+1]-patdiff[l]))
    }
  }
  return(inver)
}


## FACTORADIC function:
#' @keywords internal
factoradic <- function(M,N) {
  f = array(data = 0, dim=c(1,N))
  jj = 2
  
  while (M != 0) {
    f[,N-jj+1] <- M%%jj
    M <- floor(M/jj)
    jj <- jj+1
  }
  
  return(f)
}




# lehmerperm:
#' @keywords internal
lehmerperm = function(N,M) {
  # lehmerperm - obtain the M-th permutation of (1:N)
  #   perm = lehmerperm(N,M) returns the m-th permutation of the sorted list
  #   of all permutations from PERMS, where M=1 corresponds to identity 
  #   permutation. N, M are non-negative scalar, perm has size 1-by-N.
  #
  #   See also PERMS
  #        and NPERMUTEK, RECPERMS, NEXTPERM, PERMS1 on the File Exchange
  #
  # Algorithm: For given N and M, where 1 <= M <= N!, the M-th 
  #            permutation of N objects is closely related to the 
  #            factoradic of M; see factoradic on the File Exchange.
  #            To convert the factoradic into a permutation follow these
  #            steps
  #            
  #            For decreasing i
  #                If element(j)>=element(i) ; where j>i
  #                    element(j) increase by one.
  #
  #            The result will be a permutation of (0:N-1).
  #            Add 1 to yield the permutation of (1:N).
  
  # for Matlab (should work for most versions)
  # version 1.0 (Feb 2009)
  # (c) Darren Rowland
  # email: darrenjrowland@hotmail.com
  #
  # Keywords: single permutation
  
  # MatLab-code translated into R by Sebastian Sippel, 11.10.2013
  
  #error(nargchk(2,2,nargin));
  nargin <- length(as.list(match.call()))-1
  if(nargin != 2) stop("Wrong number of input elements")
  
  if(length(N) != 1 | N <= 0 | round(N, digits=0) != N) stop("The first input has to be a non-negative integer")
  
  if(length(M) != 1 | M <= 0 | round(M, digits=0) != M) stop("The second input has to be a non-negative integer")
  
  if(M > factorial(N)) stop("M should not exceed N!")
  
  # convert M to zero-based
  M = M-1;
  perm = factoradic(M,N)
  
  for (ii in (N-1):1) {
    
    ##Alternative Schleife:
    #     for jj = ii+1:N
    #         if(perm(jj)>=perm(ii))
    #             perm(jj) = perm(jj) + 1;
    #         end
    #     end
    
    #for (ii in (N-1):1) {
    #for (jj in (ii+1):N) {
    #  if(perm[,jj] >= perm [,ii]) perm[,jj] = perm[,jj] +1
    #}
    #}
    
    #perm(ii+1:N) = perm(ii+1:N) + (perm(ii+1:N)>=perm(ii)); -> original vector
    
    perm[,(ii+1):N] <- ifelse (perm[,(ii+1):N] >= perm[,ii], perm[,(ii+1):N]+1, perm[,(ii+1):N])
  }
  # convert permutation to one-based (from zero-based)
  perm = perm + 1
  
  return(perm)
}
