### Fisher-Rao Geometry ###

# exponential map of square-root density functions
# v is an element in the tangent space of \psi_1
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

# inverse exponential map of square-root density functions
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

