#' Function for discrete clr transformation.
#' 
#' @param  X A numeric vector, the grid over which `dens` is evaluated.
#' @param  dens A numeric vector, the values of the density evaluated over `X`.
#' 
#' @return A numeric number, the clr-transformed density.
#' @export
#' 
#######################################################
dens2clr = function(X, dens) {
  return(log(dens)-trapzRcpp(X = X, Y = log(dens))/diff(range(X)))
}

#######################################################
#######################################################

#' Function for inverse clr transformation.
#' 
#' @param  X A numeric vector, the grid over which `clr` is evaluated.
#' @param  clr A numeric vector, the values of the clr-transformed density evaluated over `X`.
#' 
#' @return A numeric number, the original density.
#' @export
#' 
#######################################################
# back-transform a clr to a density
clr2dens = function(X, clr) { 
  if(is.fd(clr))
    return(exp(eval.fd(X,clr))/trapzRcpp(X,exp(eval.fd(X,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzRcpp(X,exp(clr)))
}

#######################################################
#######################################################
