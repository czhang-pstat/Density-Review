#' Function for deregularizing density functions (GPCA analysis compatible).
#' 
#' This is a modified version of the function fdadensity::DeregulariseByAlpha().
#' The output density is normalized in a way that is compatible with the GPCA analysis
#' in this project.
#' 
#' @param  x A numeric vector, the support of a density function.
#' @param  y A numeric vector, the density values f(x) corresponding to the support x.
#' @param  alpha The minimum density value allowed, by default alpha = 0.
#' 
#' @return A numeric vector, deregularized density function.
#' @export
#' 
#######################################################
DeregulariseByAlpha_GPCA <- function(x,y,alpha=0){
  
  dx = diff(x)[1]
  l = length(y)
  
  if(min(y) < 0){
    stop('y must be a density!')
  }
  
  if(min(y) <= alpha || all(diff(y) == 0)){
    y
  }else {
    
    totalBumpArea = dx*l*min(y)
    if(totalBumpArea >= 1){
      stop('alpha deregularisation does not result in valid density.')
    }else {
      gam = (min(y) - alpha)/(1 - totalBumpArea)
      return((1 + gam*dx*l)*y - gam)
    }
  }
}

#######################################################
#######################################################

#' Function for regularizing density functions (GPCA analysis compatible).
#' 
#' This is a modified version of the function fdadensity::RegulariseByAlpha().
#' The output density is normalized in a way that is compatible with the GPCA analysis
#' in this project.
#' 
#' @param  x A numeric vector, the support of a density function.
#' @param  y A numeric vector, the density values f(x) corresponding to the support x.
#' @param  alpha The minimum density value allowed, by default alpha = 0.01.
#' 
#' @return A numeric vector, regularized density function.
#' @export
#'
#######################################################
RegulariseByAlpha_GPCA <- function(x,y,alpha=0.01){
  
  dx = diff(x)[1]
  l = length(y)
  
  if(min(y) < 0){
    stop('y must be a density!')
  }
  
  if(min(y) >= alpha || all(diff(y) == 0)){
    y
  }else {
    
    totalBumpArea = dx*l*alpha
    if(totalBumpArea >= 1){
      stop('alpha regularisation does not result in valid density.')
    }else {
      gam = (min(y) - alpha)/(totalBumpArea-1)
      return((y + gam)/(1 + gam*dx*l))
    }
  }
}
#######################################################
#######################################################
