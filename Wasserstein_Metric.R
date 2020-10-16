#' Function for calculating the Wasserstein distance between two distributions.
#' 
#' @param  quantile_1 A numeric vector, the values of input quantile 1 evaluated at quantile_grid.
#' @param  quantile_2 A numeric vector, the values of input quantile 2 evaluated at quantile_grid.
#' @param  quantile_grid A numeric vector, the common grid of the input quantiles.
#' 
#' @return A numeric number, the Wasserstein distance between quantile 1 and quantile 2.
#' @export
#' 
#######################################################

Wasserstein_Metric <- function(quantile_1, quantile_2, quantile_grid){
  dist_pw <- (quantile_1 - quantile_2)^2
  dist_f <- splinefun(quantile_grid*1e6, y = dist_pw, method = "natural")
  dist <- sqrt((integrate(dist_f, 0, 1e6, stop.on.error = FALSE)$value))
  dist <- dist/(1e3)
  return(dist)
}

#######################################################
#######################################################














