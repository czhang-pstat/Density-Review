#' Function for implementing the Exponential mapping under the Fisher-Rao geometry.
#' 
#' @param  dSup A numeric vector, the density support.
#' @param  psi_1 The Frechet mean under the Fisher-Rao geometry.
#' @param  v A numeric vector, an element in the tangent space of psi_1.
#' 
#' @return A numeric value, the squared-root of the density that corresponds to v.
#' 
#' @export
#' 
#######################################################

# exponential map of square-root density functions
# v is 
fr_Exp_Map = function(dSup, v, psi1) {
  norm = sqrt(trapzRcpp(dSup, v^2))
  
  # A threshold to prevent numerical instability:
  # if the norm of the tangent vector is small, then
  # the square root of the reference density is returned.
  if (norm < 1e-5) {
    result_scaled = psi1
  }
  
  else {
    result = cos(norm)*psi1 + sin(norm)*v/norm
    result_scaled = sqrt(result^2/trapzRcpp(dSup, result^2))
  }
  return(result_scaled)
}

#######################################################
#######################################################

#' Function for implementing the inverse Exponential map of square-root density functions.
#' 
#' @param  dSup A numeric vector, the density support.
#' @param  psi_1 The Frechet mean under the Fisher-Rao geometry.
#' @param  psi_2 A numeric vector, a square-root density function.
#' 
#' @return A numeric vector, the tangent element that corresponds to psi_2.
#' 
#' @export
#' 
#######################################################
# 
# this function lifts \psi_2 to the tangent space of \psi_1
fr_Log_Map = function(dSup, psi1, psi2) {
  u = psi2 - trapzRcpp(dSup, psi1*psi2)*psi1
  
  # A threshold to prevent numerical instability:
  # if the norm of u is small, then return the 0-vector
  # in the tangent space.
  if(trapzRcpp(dSup, u^2) > 1e-5) {
  result = u*acos(trapzRcpp(dSup, psi1*psi2))/sqrt(trapzRcpp(dSup, u^2))
  }
  
  else {result = rep(0, length(u))}
  return(result)
}

#######################################################
#######################################################

#' Function for implementing the Karcher mean.
#' 
#' @param  obs A matrix that contains the observed density functions.
#' @param  dSup A numeric vector, the density support.
#' @param  tolerance A numeric value, error threshold for algorithm convergence.
#' @param  step_size A numeric value, the descending step size.
#' 
#' @return A numeric vector, the Karcher mean.
#' 
#' @export
#' 
#######################################################

# function for finding Karcher mean
# obs: observations, rows are observations
Karcher_mean = function(obs, dSup, tolerance, step_size) {
  obs = sqrt(obs)
  start = obs[1,] # Random initialization.
  tngt_mean = apply(apply(obs, 1, function(xx) fr_Log_Map(dSup, start, xx)), 1, mean)
  error = sqrt(trapzRcpp(dSup, tngt_mean^2))
  counter = 1
  
  print(paste("Initial error:", error))
  
  while(error > tolerance & counter < 5001) {
    start = fr_Exp_Map(dSup, step_size*tngt_mean, start)
    tngt_mean = apply(apply(obs, 1, function(xx) fr_Log_Map(dSup, start, xx)), 1, mean)
    error = sqrt(trapzRcpp(dSup, tngt_mean^2))
    counter = counter + 1
  }
  print(paste("Terminating error:", error))
  return(start)
}

