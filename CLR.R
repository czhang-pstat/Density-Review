# transform a density to a clr
dens2clr = function(X, dens) {
  return(log(dens)-trapzRcpp(X = X, Y = log(dens))/diff(range(X)))
}

# back-transform a clr to a density
clr2dens = function(X, clr) { 
  if(is.fd(clr))
    return(exp(eval.fd(X,clr))/trapzRcpp(X,exp(eval.fd(X,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzRcpp(X,exp(clr)))
}
