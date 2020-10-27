#' Function for estimating pdfs (kernel smoothing).
#' 
#' @param  ret A list that contains data points.  For example, ret[[i]] is the data points sampled from the ith density.
#' @param  m The number of grid points over which densites are evaluated.
#' @param  band_choice Bandwidth selection for kernel density estimation, "Silverman" is the Rule-of-Thumb bandwidth (zero-stage plug-in), "DPI" is the Direct Plug-In bandwidth (two-stage plug-in).
#' @param  kernel "gaussian" or "epanechnikov".
#' 
#' @return A list that contains a matrix of estimated densities and the density support.
#' 
#' @export
#' 
#######################################################

density_estimation <- function(ret, m = 500,
                               band_choice = c("Silverman", "DPI"),
                               kernel = c("gaussian", "epanechnikov")) {
  ### sample size
  N <- length(ret)
  
  ### choice of kernel and bandwidths
  kernel <- match.arg(kernel)
  band_choice <- match.arg(band_choice)
  
  if(band_choice == "Silverman") {
    if(kernel == "gaussian") {
      h.hat <- sapply(1:N, function(t) 1.06*sd(ret[[t]])*(length(ret[[t]])^(-(1/5))))
    }
    else {
      h.hat <- sapply(1:N, function(t) 2.34*sd(ret[[t]])*(length(ret[[t]])^(-(1/5))))
    }
  }
  else {
    if(kernel == "gaussian") {
      h.hat <- sapply(1:N, function(t) dpik(ret[[t]], kernel = "normal"))
    }
    else {
      h.hat <- sapply(1:N, function(t) dpik(ret[[t]], kernel = "epanech"))
    }
  }
  
  ### density support
  a <- min(sapply(1:N, function(ik) min(ret[[ik]])))
  b <- max(sapply(1:N, function(ik) max(ret[[ik]])))
  u <- seq(from = a, to = b, length = m)
  du <- u[2] - u[1]
  
  ### creating an (m x n) matrix which represents the observed densities
  # Y[j,t] is the density at day t evaluated at u[j]
  if(kernel == "gaussian") {
    Y <- sapply(1:N, function(t) density(ret[[t]], bw = h.hat[t], kernel = 'gaussian', from = a, to = b, n = m)$y)
  }
  if(kernel == "epanechnikov") {
    Y <- sapply(1:N, function(t) density(ret[[t]], bw = h.hat[t], kernel = 'epanechnikov', from = a, to = b, n = m)$y)
  }
  
  for(t in 1:N) {
    Y[,t] = Y[,t]/(sum(Y[,t]) * du)
  }
  
  ### Dealing with zero values
  # inflate densities to single out extreme values
  Y.inflate <- Y*1e16
  epsilon = sapply(1:N, function(xx) max(Y.inflate[,xx] - round(Y.inflate[,xx], 2)))
  
  Y.inflate.matrix = matrix(NA, m, N)
  for(i in 1:N) {
    index = which(round(Y.inflate[,i], 2) == 0)
    if(epsilon[i] == 0) {
      bumper <- 1e-3
    }
    else {
      bumper <- epsilon[i]
    }
    Y.inflate.matrix[,i] <- replace(Y.inflate[,i], index, bumper)
    Y.inflate.matrix[-index,i] <- Y.inflate[-index, i] * (1 - (length(index) * bumper)/(10^16))
  }

  dens = sapply(1:N, function(xx) Y.inflate.matrix[,xx]/trapzRcpp(X = u, Y = Y.inflate.matrix[,xx]))
  
  return(list(dens, u))
}
#######################################################
#######################################################

