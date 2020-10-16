#' Perform Frechet Regression with the Wasserstein Distance
#' @export
#' @description
#' @param rightside_formula a right-side formula
#' @param Xfit_df n-by-p matrix (or dataframe) of predictor values for fitting (do not include a column for the intercept)
#' @param Ymat one of the following matrices:
#' \itemize{
#' \item{if Ytype = 'quantile'} an n-by-m matrix of the observed quantile functions. Ymat[i, :] is a 1-by-m vector of quantile function values on grid \code{Sup}.
#' \item{if Ytype = 'density'} an n-by-m matrix of the observed density functions. Ymat[i, :] is a 1-by-m vector of density function values on grid \code{Sup}.
#' }
#' @param bandwidth the smoothing parameter in kernel smoothing applied on the density function, a value in (0, 1). (default: 0.01)
#' @param Ytype 'quantile' or 'density'
#' @param Sup a length m vector - common grid for all density functions in Ymat. (default: seq(0, 1, length.out = ncol(Ymat)))
#' @param t a length m vector - common grid for quantile functions.
#' @return a list containing the following ojects:
#' \item{qfit}{n-by-m matrix of fitted quantile density functions.}
#' \item{ffit}{n-by-m matrix of fitted density functions.}
#' \item{Qfit}{n-by-m matrix of fitted quantile functions.}
#' \item{call}{function call.}
#' \item{predictor_names}{names of predictors as the colnames given in the xfit matrix or dataframe.}
#' \item{xfit}{design matrix in quantile fitting.}
#' \item{lm_res}{the fitting model}
#' \item{rformula}{\code{rightside_formula}}
#' \item{Ysmooth}{a list containing the following matrices:
#' \itemize{
#' \item{Qobs:} {n-by-m matrix of the observed quantile functions.}
#' \item{qobs:} {n-by-m matrix of the observed quantile density functions.}
#' \item{fobs:} {n-by-m matrix of the observed density functions.}
#' }}
#' \item{t_equal}{a length m vector - common grid for all quantile functions in Qobs.}
#' @examples
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#'
#' res1 = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, bandwidth = 0.01, Ytype = 'density', Sup = dSup)
#' res2 = wass_regress(rightside_formula = ~ log_b_vol * weight, Xfit_df = predictor, Ymat = densityCurves, bandwidth = 0.01, Ytype = 'density', Sup = dSup)
#' @references
#'   \cite{Wasserstein F-tests and confidence bands for the Frechet regression of density response curves, Alexander Petersen, Xi Liu and Afshin A. Divani, 2019}
wass_regress <- function(rightside_formula, Xfit_df, Ymat, bandwidth = 0.01, Ytype = 'density', Sup = NULL, t = NULL) {

        ### ======================  input check  ====================== ###

        if (nargs() < 3) {
                stop("Not enough arguments: rightside_formula, Xfit_df and Ymat are required")
        }
        if (is.null(Sup)) Sup = seq(0, 1, length.out = ncol(Ymat))

        if (ncol(Ymat) != length(Sup)) {
                stop("ncol(Ymat) and length(Sup) do not match")
        }
        if (Ytype != 'density' & (min(Sup) != 0 | max(Sup) != 1)) {
                stop("Input Sup should be an increasing grid beginning at 0 and ending at 1")
        }
        if (!Ytype %in% c('density', 'quantile')) {
                stop("Ytype must be \"density\" or \"quantile\" ")
        }

        if (Ytype != 'density' & !all(diff(t(Ytype)) >= 0)) {
                stop("Each row of Ymat should be nondecreasing")
        }

        # Added formula input check
        if (!inherits(rightside_formula, 'formula')) {
                stop("formula argument must be of class \"formula\" ")
        }

        ### ====================== begin function  ====================== ###
        cl <- match.call()
        if (is.null(t)) t = seq(0, 1, length.out = length(Sup))

        ### convert Ymat
        if (Ytype == 'density') {
                smoothY = den2Q_qd(densityCurves = Ymat, dSup = Sup, bandwidth = bandwidth, t = t)
        } else {
                stop('currently no support for quantile input')
                smoothY = quan2den_qd(QMat_raw = Ymat, t = Sup, bandwidth = 0.1)
                t = seq(0, 1, length.out = length(Sup))
        }

        Qmat = smoothY$Qobs

        ### OLS fit of Qmat
        twoside_formula <- as.formula(paste('Qmat', deparse(rightside_formula)))

        Qmat_lm <- stats::lm(twoside_formula, data = Xfit_df)
        xfit <- model.matrix(Qmat_lm)[ ,-1]
        Qmat_fitted = fitted(Qmat_lm)

        ### quadratic program to make qfitted > 0
        m = length(t)
        t_equal = seq(0, 1, length.out = m)
        Qmat_fitted = quadraticQ(Qmat_fitted, t_equal)

        ### ======================   return list   ====================== ###

        res = structure(list(
                   call = cl,
                   rformula = rightside_formula,
                   predictor_names = colnames(xfit),

                   Qfit = Qmat_fitted,
                   xfit = xfit,
                   Xfit_df = Xfit_df,
                   Ysmooth = smoothY,
                   t = t,
                   t_equal = t_equal),
                   class = 'WRI')
}
