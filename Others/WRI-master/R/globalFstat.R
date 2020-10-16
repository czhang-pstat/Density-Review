#' An internal function used in bootstrap global F test
#' @name globalFstat
#' @keywords internal
globalFstat <- function(xfit, Qobs, t) {

        Qmat_lm <- stats::lm(Qobs ~ xfit)
        Qfitted = fitted(Qmat_lm)
        n = nrow(Qobs)
        Qmean = colMeans(Qobs)
        ssr = sum(sapply(1:n, function(i) fdapace::trapzRcpp(t, (Qmean - Qfitted[i, ])^2)))
        return(ssr)
}
