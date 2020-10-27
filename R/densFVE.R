#' Function for calculating the Wasserstein distance between two distributions.
#' 
#' @param  dens A matrix of sample densities.
#' @param  fitted_dens A matrix of fitted densities.
#' @param  metric A string indicates the metric to evaluate FVE: "L2", "Wasserstein", "FisherRao", "Aitchison".
#' @param  dSup A numeric vector, the support of the density.
#' @param  mean_psi A numeric vector, the Karcher mean.  It is important to feed the Karcher mean into the function
#' instead of calculate it from within to improve efficiency (especially in the split sample FVE analysis).
#' 
#' @return A numeric number, the FVE.
#' @export
#' 
#######################################################
densFVE = function(dens, fitted_dens, metric, dSup, mean_psi) {
  if (metric == "L2") {
    mu = colMeans(dens)
    vtotal = mean(apply(dens, 1, function(xx) trapzRcpp(X = dSup, Y = (xx-mu)^2)))
    
    vK <- mean(apply((dens -  fitted_dens)^2, 1, function(xx) trapzRcpp(X = dSup, Y = xx)))
    FVE = (vtotal - vK)/vtotal
  }
  
  if (metric == "Wasserstein") {
    qSup = seq(0,1, length.out = length(dSup))
    quantile = t(apply(dens, 1, function(xx) dens2quantile(xx, dSup = dSup)))
    mu_quantile = colMeans(quantile)
    
    vtotal = mean(apply(quantile, 1, function(xx) trapzRcpp(X = qSup, Y = (xx - mu_quantile)^2)))
    
    fitted_quantile = t(apply(fitted_dens, 1, function(xx) dens2quantile(xx, dSup = dSup)))
    vK = mean(apply((quantile -  fitted_quantile)^2, 1, function(xx) trapzRcpp(X = qSup, Y = xx)))
    FVE = (vtotal - vK)/vtotal
  }
  
  if (metric == "FisherRao") {
    psi = sqrt(dens)
    fitted_psi = sqrt(fitted_dens)

    #mean_psi = Karcher_mean(dens, dSup, tolerance = 1e-6, step_size = 1e-2)
    vtotal = mean(apply(psi, 1, function(xx) (acos(trapzRcpp(X = dSup, Y = xx*mean_psi)))^2))
    
    inner_prod = matrix(NA, nrow = nrow(fitted_dens), ncol = ncol(fitted_dens))
    
    for(i in 1:nrow(fitted_dens)) {
      inner_prod[i,] = unlist(psi[i,]*fitted_psi[i,])
    }
    
    vK = mean(apply(inner_prod, 1, function(xx) (acos(trapzRcpp(X = dSup, Y = xx)))^2))
    FVE = (vtotal - vK)/vtotal
  }
  
  if (metric == "Aitchison") {
    clr = t(apply(dens, 1, function(xx) dens2clr(dSup, xx)))
    fitted_clr = t(apply(fitted_dens, 1, function(xx) dens2clr(dSup, xx)))
    
    mu = colMeans(clr)
    vtotal = mean(apply(clr, 1, function(xx) trapzRcpp(X = dSup, Y = (xx-mu)^2)))
    
    vK <- mean(apply((clr -  fitted_clr)^2, 1, function(xx) trapzRcpp(X = dSup, Y = xx)))
    FVE = (vtotal - vK)/vtotal
  }
  
  return(FVE)
}

