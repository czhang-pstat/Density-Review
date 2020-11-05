#' Function for numerical differentiation for a function y = f(x).
#' 
#' @param  x A numeric vector, the grid over which f'(x) is evaluated.
#' @param  y A numeric vector, the values of f(x) corresponding to the x grid.
#' 
#' @return A numeric vector f'(x).
#' @export
#' 
#######################################################

Finite_Diff = function(x, y) {
  n = length(x)
  dx = vector(length = n)
  
  # finite difference
  for (i in 2:n) {
    dx[i-1] = (y[i] - y[i-1]) / (x[i] - x[i-1])
  }
  
  # set the last derivative equal to the adjacent value
  dx[n] = dx[n-1]
  
  return(dx)
}

#######################################################
#######################################################
