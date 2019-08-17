# -----------------------------------
## Generate k-noise by Holger Lange
# -----------------------------------

## This program is thoroughly checked with the knoise_fft Matlab function. 
#  [knoise_fft in Matlab is thoroughly checked with knoise_periodogram, which
#   is a slightly altered version of knoise_erzeugen (changed to produce 501
#   instead of 513 output frequencies.]
#  The first deviation (between the R and knoise_fft program) occurs in 
#  the eight decimal place. However, the first element of the output vector psp 
#  (spectral density estimates) remains to be checked with MatLAB!! (psp[1])
#  Sebastian Sippel, 09. November 2013

#' @title A function to generate k-noise
#' @export
#' @description Generates samples of power law noise. 
#' @usage powernoise(k, N)
#' @param k Power law scaling exponent
#' @param N number of samples to generate
#' @details 
#' Generates samples of power law noise. 
#' The power spectrum of the signal scales as f^(-k). The R function uses fft(),
#'  similarly to the knoise_fft Matlab function.
#' @return A named list with three entries is returned. 
#' x - N x 1 vector of power law samples
#  psp - periodogram values (squared amplitudes)
#  f - corresponding frequencies (units 1/time step)
#' @author Sebastian Sippel and Holger Lange
#' @examples
#' powernoise_series = powernoise(k=2, N=10000)
powernoise <- function(k, N) {
  #  Generate samples of power law noise. The power spectrum
  #  of the signal scales as f^(-k). The R function uses fft(),
  #  similarly to the knoise_fft Matlab function.
  # 
  #  Useage:
  #  x = powernoise(alpha, N)
  # 
  #  Inputs:
  #   k - power law scaling exponent
  #   N     - number of samples to generate
  # 
  #  Output:
  #  x     - N x 1 vector of power law samples
  #  psp - periodogram values (squared amplitudes)
  #  f - corresponding frequencies (units 1/time step)
  # 
  #  This version H. Lange Oct. 2010 
  #  Sebastian Sippel, November 2013
  
  N2 = floor(N/2)-1;
  f = (2:(N2+1))
  A2 = 1/(f^(k/2))
  
  p2 = complex(real=matrix(stats::rnorm(N2*1,0,1),N2,1), imaginary=matrix(stats::rnorm(N2*1,0,1),N2,1))  #p2 = randn(N2,1) + i * randn(N2,1);
  d2 = A2*p2;  # d2 = A2.*p2;
  
  d = c(1,d2,1/((N2+2)^k),rev(Conj(d2)))     #d = [1; d2; 1/((N2+2)^k); flipud(conj(d2))];
  
  x = Re(stats::fft(d, inverse=T))  #x = real(ifft(d));
  
  xdft = stats::fft(x)
  xdft = xdft[1:(N/2+1)]
  psd_fft = (1/(2*pi*N))*abs(xdft)^2
  psd_fft[2:(length(psd_fft)-1)] = 2*psd_fft[2:(length(psd_fft)-1)]
  f_fft = seq(from=0, to=pi, by=(2*pi)/N)
  
  
  # make series zero mean and var = 1, and rescale the power spectrum accordingly
  psp=psd_fft/stats::var(x);
  x_scaled = (x - mean(x))/stats::sd(x)
  #     x = ((x - min(x))/(max(x) - min(x)) - 0.5) * 2;
  
  ret.list = list(x_scaled,psp,f_fft)
  names(ret.list) = c("x", "psp", "f_fft")
  
  return(ret.list)
}
