## Create time series of chaotic maps to test Information Theory Quantifiers
# ---------------------------------------------

### Logistic Map
#' @title A function to generate a time series from the logistic map
#' @export
#' @description Generates a time series from the logistic map
#' @usage logistic_map(N, r, start="rand", disregard_N=0)
#' @param N length of the time series that is to be generated
#' @param r logistic map parameter, must be in the range [0,4]
#' @param start start value. Default is to random.
#' @param disregard_N Number of values at the beginning of the series to disregard  
#' @return A vector of length N
#' @references May, R.M., 1976. Simple mathematical models with very complicated dynamics. Nature, 261(5560), pp.459-467.
#' @author Sebastian Sippel
#' @examples
#' logistic_map(N = 10^4, r=4)
logistic_map <- function(N, r=4, start="rand", disregard_N=0) {
  
  # check startvalue start
  if (start == "rand")  {
    start = stats::runif(1, min=0, max=1) }
  if (start >= 1 | start <= 0) {
    print("Start value is not in the range [0,1]")
    return()
  }
  
  # check r
  if (r > 4.000000001 | r < 0) {
    print("r must be in the range [0,4]")
    return()
  }
  
  x <- numeric(length=N+disregard_N)
  x[1] <- start
  for (i in 1:(length(x)-1)) {
    x[i+1] <- r*x[i]*(1-x[i])  
  }
  x_ts <- x[(disregard_N+1):length(x)]
  return(x_ts)
}


### Tent Map (see e.g. Feldman et al 2008 or Wikipedia):
#' @title A function to generate a time series from the logistic map
#' @export
#' @description Generates a time series from the logistic map
#' @usage tent_map(N, mu, start="rand", disregard_N=0)
#' @param N length of the time series that is to be generated
#' @param mu Tent map parameter, must be in the range [0,2]
#' @param start start value. Default is to random.
#' @param disregard_N Number of values at the beginning of the series to disregard  
#' @return A vector of length N
#' @references Feldman, D.P., McTague, C.S. and Crutchfield, J.P., 2008. The organization of intrinsic computation: Complexity-entropy diagrams and the diversity of natural information processing. Chaos: An Interdisciplinary Journal of Nonlinear Science, 18(4), p.043106.
#' @author Sebastian Sippel
#' @examples
#' tent_map(N = 10^4, mu=1.8)
tent_map <- function(N, mu = 2, start="rand", disregard_N=0) {
  
  # check startvalue start
  if (start == "rand")  {
    start = stats::runif(1, min=0, max=1) }
  if (start >= 1 | start <= 0) {
    print("Start value is not in the range [0,1]")
    return()
  }
  
  # check mu
  if (mu > 2 | mu < 0) {
    print("mu must be in the range [0,2]")
    return()
  }
  
  x <- numeric(length=N+disregard_N)
  x[1] <- start
  for (i in 1:(length(x)-1)) {
    if (x[i] < 0.5) {
      x[i+1] <- x[i] * mu
    } else {
      x[i+1] <- mu * (1 - x[i])
    } 
  }
  x_ts <- x[(disregard_N+1):length(x)]
  return(x_ts)
}


### Skew Tent Map
#' @title A function to generate a time series from the logistic map
#' @export
#' @description Generates a time series from the Skew-Tent map
#' @usage skew_tent_map(N, a, start="rand", disregard_N=0)
#' @param N length of the time series that is to be generated
#' @param a Skew-Tent map parameter, must be in the range [0,1]
#' @param start start value. Default is to random.
#' @param disregard_N Number of values at the beginning of the series to disregard  
#' @return A vector of length N
#' @references Schuster, H.G., 1988. Deterministic chaos. An Introduction.
#' @author Sebastian Sippel
#' @examples
#' skew_tent_map(N = 10^4, a=0.1847)
skew_tent_map <- function(N, a=0.1847, start="rand", disregard_N=0) {
  
  # check startvalue start
  if (start == "rand")  {
    start = stats::runif(1, min=0, max=1) }
  if (start >= 1 | start <= 0) {
    print("Start value is not in the range [0,1]")
    return()
  }
  
  # check a
  if (a > 1 | a < 0) {
    print("a must be in the range [0,1]")
    return()
  }
  
  x <- numeric(length=N+disregard_N)
  x[1] <- start
  for (i in 1:(length(x)-1)) {
    if (x[i] < a) {
      x[i+1] <- x[i]/a
    } else {
      x[i+1] <- (1-x[i])/(1-a)
    } 
  }
  x_ts <- x[(disregard_N+1):length(x)]
  return(x_ts)
}

### Henon Map (see e.g. Olivares et al 2012)
# N length of the time series that should be generated
# a parameter a;
# b parameter b;
# start: startvalue
#disregard_N: number of values after thestart of the series that should be disregarded:

### Henon Map
#' @title A function to generate a time series from the Henon Map
#' @export
#' @description Generates a time series from the Henon map
#' @usage henon_map(N, a, b, startx="rand", starty="rand", disregard_N=0)
#' @param N length of the time series that is to be generated
#' @param a Henon map parameter a
#' @param b Henon map parameter b
#' @param startx start value in x direction. Default is to random.
#' @param starty start value in y direction. Default is to random.
#' @param disregard_N Number of values at the beginning of the series to disregard  
#' @return A vector of length N
#' @references Schuster, H.G., 1988. Deterministic chaos. An Introduction.
#' @author Sebastian Sippel
#' @examples
#' henon_map(N = 10^4, a=1.4, b=0.3)
henon_map <- function(N, a = 1.4, b = 0.3, startx="rand", starty="rand", 
                  disregard_N=0) {
  
  # check startvalue startx
  if (startx == "rand")  {
    startx = stats::runif(1, min=0, max=1) }
  if (startx >= 1 | startx <= 0) {
    print("Start value x is not in the range ]0,1[")
    return()
  }
  
  # check startvalue starty
  if (starty == "rand")  {
    starty = stats::runif(1, min=0, max=1) }
  if (starty >= 1 | starty <= 0) {
    print("Start value y is not in the range ]0,1[")
    return()
  }
  
  x <- numeric(length=(N+disregard_N))
  y <- numeric(length=(N+disregard_N))
  x[1] <- startx
  y[1] <- starty
  
  for (i in 1:(length(x)-1)) {
    x[i+1] <- 1 - a * (x[i])^2 + y[i]
    y[i+1] <- b * x[i]
  }
  x_ts <- x[(disregard_N+1):length(x)]
  y_ts <- y[(disregard_N+1):length(y)]
  return(data.frame(x_ts,y_ts))
}







### Schuster Map (see Rosso et al 2007)
# N length of the time series that should be generated
# z exponent parameter z; 
## Rosso et al (2007) use z=5/2, z=2, and z=3/2
## Olivares et al (2012) use 1.25 =< z =< 2 with dz = 0.25
# start: startvalue
#disregard_N: number of values after thestart of the series that should be disregarded:
# Olivares et al (2012) disregard 10^5 iterations before they generate the chaotic series

### Schuster Map
#' @title A function to generate a time series from the Schuster Map
#' @export
#' @description Generates a time series from the Schuster map
#' @usage schuster_map(N, z, start="rand", disregard_N=0)
#' @param N length of the time series that is to be generated
#' @param z Schuster map parameter
#' @param start start value. Default is to random.
#' @param disregard_N Number of values at the beginning of the series to disregard  
#' @return A vector of length N
#' @references Schuster, H.G., 1988. Deterministic chaos. An Introduction.
#' @author Sebastian Sippel
#' @examples
#' schuster_map(N = 10^4, z=2)
schuster_map <- function(N, z=2, start="rand", disregard_N=0) {
  
  # check startvalue start
  if (start == "rand")  {
    start = stats::runif(1, min=0, max=1) }
  if (start >= 1 | start <= 0) {
    print("Start value is not in the range ]0,1[")
    return()
  }
  
  # check z
  if (z > 5/2 | z < 1.25) {
    print("z should be in the range [1.25,2.5]")
    return()
  }
  
  
  x <- numeric(length=N+disregard_N)
  x[1] <- start
  for (i in 1:(length(x)-1)) {
    x[i+1] <- (x[i] + (x[i])^z)%%1
  }
  x_ts <- x[(disregard_N+1):length(x)]
  return(x_ts)
}





### Quadratic Map
#' @title A function to generate a time series from the Quadratic map
#' @export
#' @description Generates a time series from the Quadradtic map
#' @usage quadratic_map(N, k, start="rand", disregard_N=0)
#' @param N length of the time series that is to be generated
#' @param k Quadratic map parameter
#' @param start start value. Default is to random.
#' @param disregard_N Number of values at the beginning of the series to disregard  
#' @return A vector of length N
#' @references Grebogi, C., Ott, E. and Yorke, J.A., 1983. Crises, sudden changes in chaotic attractors, and transient chaos. Physica D: Nonlinear Phenomena, 7(1-3), pp.181-200.
#' @author Sebastian Sippel
#' @examples
#' quadratic_map(N = 10^4, k=1.4)
quadratic_map <- function(N, k, start="rand", disregard_N=0) {
  x = numeric(length=N+disregard_N)
  x[1] = stats::runif(1)
  for (i in 2:N) {
    x[i] = k - x[i-1]^2
  }
  return(x[(disregard_N+1):(length(x))])
}


