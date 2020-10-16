#' convert density function to quantile and quantile density function
#' @export
#' @description
#' @param densityCurves n-by-m matrix of density curves
#' @param dSup length m vector contains the common support grid of the density curves
#' @param bandwidth smoothing parameter used in kernel smooth
#' @param t the desired grid points in [0, 1] for the quantile function.
den2Q_qd <- function(densityCurves, dSup, bandwidth = 0.05, t = seq(0, 1, by = 0.01)) {

        n = nrow(densityCurves)
        m = length(t)
        res = list()
        m_dSup = length(dSup)

        if(min(densityCurves) < 0) {
                stop("Please correct negative or zero probability density estimates.")
        }

        if (length(dSup) < 25) {
                stop("please give densely observed density curves, with length(dSup) >= 25 ")
        }

        if (bandwidth == 0)(
                # set the bandwidth small enough so that ksmooth function does nothing to densityCurves
                bandwidth = 1e-20
        )

        # test t is equally spaced or not
        min_max <- range(diff(t)) / mean(diff(t))
        equal_logi = isTRUE(all.equal(min_max[1], min_max[2]))

        res_list_dSup = sapply(1:n, function(i) {

                smooth_fit = ksmooth(dSup, densityCurves[i, ], kernel = "normal", bandwidth = diff(range(dSup))*bandwidth, x.points = dSup)
                dens = smooth_fit$y

                if (abs(fdapace::trapzRcpp(X = dSup, dens) - 1) > 1e-05) {
                        # warning("Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.")
                        dens = dens/fdapace::trapzRcpp(X = dSup, Y = dens)
                }

                cdf_raw = fdapace::cumtrapzRcpp(X = dSup, Y = dens)
                spline_fit = splinefun(x = dSup, y = cdf_raw, method = "hyman")
                fitted_cdf = spline_fit(dSup, deriv = 0)
                fitted_pdf = spline_fit(dSup, deriv = 1)
                fitted_pdf_prime = spline_fit(dSup, deriv = 2)
                ### avoid pdf is 0
                if(any(fitted_pdf==0)) {
                        index = which(fitted_pdf == 0)
                        #for (i in index) {
                        #        fitted_pdf[i] = min(fitted_pdf[i-1], fitted_pdf[i+1], na.rm = TRUE)
                        #}
                        fitted_pdf[index] = min(fitted_pdf[!index])
                        fitted_pdf = fitted_pdf/fdapace::trapzRcpp(X = dSup, Y = fitted_pdf)
                        fitted_cdf = fdapace::cumtrapzRcpp(X = dSup, Y = fitted_pdf)
                }

                return(cbind(cdf = fitted_cdf, pdf = fitted_pdf, pdf_prime = fitted_pdf_prime))
        }) # the first m_dSup rows are CDF, middle m_dSup rows are PDF, last m_dSup rows are PDF'

        Qobs_t = t(apply(res_list_dSup[1:m_dSup, ], MARGIN = 2, FUN = function(fitted_cdf) {
                approx(x = fitted_cdf, y = dSup, xout = t,
                       rule = c(2, 2))[[2]]
        }))

        qobs_t = t(sapply(1:n, FUN = function(i) {
                approx(x = res_list_dSup[1:m_dSup, i], y = 1/res_list_dSup[(m_dSup+1):(2*m_dSup), i], xout = t,
                       rule = c(2, 2))[[2]]
        }))

        qobs_prime_t = t(sapply(1:n, FUN = function(i) {
                approx(x = res_list_dSup[1:m_dSup, i], y = -1/res_list_dSup[(m_dSup+1):(2*m_dSup), i]^3 * res_list_dSup[(2*m_dSup+1):(3*m_dSup), i], xout = t,
                       rule = c(2, 2))[[2]]
        }))

        res$Qobs_t = Qobs_t
        res$qobs_t = qobs_t
        res$qobs_prime_t = qobs_prime_t
        res$fobs_t = 1/res$qobs_t

        #########  convert the nonequally spaced density support points into equally spaced    #########

        if (!equal_logi) {

        t_equal =  seq(min(t), max(t), length.out = m)
        Qobs = t(apply(res_list_dSup[1:m_dSup, ], MARGIN = 2, FUN = function(fitted_cdf) {
                approx(x = fitted_cdf, y = dSup, xout = t_equal,
                       rule = c(2, 2))[[2]]
        }))

        qobs = t(sapply(1:n, FUN = function(i) {
                approx(x = res_list_dSup[1:m_dSup, i], y = 1/res_list_dSup[(m_dSup+1):(2*m_dSup), i], xout = t_equal,
                       rule = c(2, 2))[[2]]
        }))

        res$Qobs = Qobs
        res$qobs = qobs
        res$fobs = 1/res$qobs

        } else {

                res$Qobs = res$Qobs_t
                res$qobs = res$qobs_t
                res$fobs = res$fobs_t

        }

        return(res)

}
