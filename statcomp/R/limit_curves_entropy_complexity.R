# This script implements the calculation of the limiting curves in the Entropy-Complexity plane:
# Sebastian Sippel
# 21.05.2015

#' @keywords internal
shannon <- function(P) {
  S_p <- 0
  for (i in 1:length(P)) {
    if (P[i] >= 10^(-30)) {
      S_p = S_p - P[i] * log(P[i])
    }
  }
  return(S_p)
}

# determine minimum complexity:
#' @keywords internal
minimum_limit_curve <- function(ndemb) {
  nk = 10000 ##number of k's to calculate (i.e. length of vector k_vary)
  N = factorial(ndemb)
  H_s <- numeric(length=nk)
  C_js <- numeric(length=nk)
  pk_vary <- seq(from=1/N, to=1, length.out=nk)
  Q0 <- -2*(((N+1)/N)*log(N+1)-2*log(2*N)+log(N))^(-1)
  Pe <- rep(1/N,N)
  
  for (k in 1:nk) {
    P <- array(data=NA,dim=c(N))
    P[1] <- pk_vary[k] #one state, varies between 1/N and 1
    P[2:length(P)] <- (1-pk_vary[k])/(N-1)
    H_s[k] <- shannon(P)/log(length(P))
    C_js[k] <- H_s[k]*Q0*(shannon((P+Pe)/2) - shannon(P)/2 - shannon(Pe)/2)
  }
  
  return(list(H_s,C_js))
}


### Determine maximum complexity:
#' @keywords internal
maximum_limit_curve <- function(ndemb) {
  nk=10000
  N <- factorial(ndemb)
  #N-1 Probability Distributions with Dimensions from 2 to N
  
  H_s <- array(data=NA, dim=c(N-1,nk))
  C_js <- array(data=NA, dim=c(N-1,nk))
  
  for (i in 1:(N-1)) {
    P <- array(data=0, dim=c(N))
    Q0 <- -2*(((N+1)/(N))*log(N+1)-2*log(2*(N))+log(N))^(-1)
    Pe <- rep(1/(N),N)
    
    pk_vary <- seq(from=0, to=1/(N), length.out=nk)
    
    for (k in 1:length(pk_vary)) {
      P[1] <- pk_vary[k]
      for (j in 1:(N-i)) {
        P[j+1] <- (1-pk_vary[k])/(N-i)
      }
      H_s[i,k] <- shannon(P)/log(N)
      C_js[i,k] <- H_s[i,k]*Q0*(shannon((P+Pe)/2) - shannon(P)/2 - shannon(Pe)/2)
    }
    print(i)
  }
  
  extract <- seq(0,1,by=0.0001)
  H_s_neu <- NULL
  C_js_neu <- NULL
  for (i in 1:(length(extract)-1)) {
    H_s_neu[i] <- H_s[which(H_s> extract[i] & H_s< extract[i+1])[1]]
    C_js_neu[i] <- C_js[which(H_s> extract[i] & H_s< extract[i+1])[1]]
  }
  
  #plot(stats::na.omit(H_s_neu),stats::na.omit(C_js_neu),type='l', xlim=c(0,1),ylim=c(0,0.5))
  return(list(stats::na.omit(H_s_neu),stats::na.omit(C_js_neu)))
}
  
  
  
#' @title Limit curves in the Entropy-Complexity plane
#' @export
#' @description Compute the limit curves in the Entropy Complexity plane
#' @usage limit_curves(ndemb, fun = "min")
#' @param ndemb Embedding dimension
#' @param fun Whether the upper (max) or lower (min) limit curve should be computed
#' @details 
#' This function returns the respective limit curve.
#' @return A list with two entries
#' @references none
#' @author Sebastian Sippel
limit_curves <- function(ndemb, fun = "min") {
  if (fun == "min") {
    xyz = (minimum_limit_curve(ndemb = ndemb))
  } else if (fun == "max") {
    xyz = (maximum_limit_curve(ndemb = ndemb))
  }
  return(xyz)
}

