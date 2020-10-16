#' Summary Function of Wasserstein Regression Model
#' @export
#' @description
#' @param wass_regress_res an object returned by the \code{wass_regress} function
#' @return a list containing the following fields:
#' \item{call}{function call of the Wasserstein regression}
#' \item{method}{methods used to compute p value, see details}
#' \item{statistic}{test statistics}
#' \item{critical_value}{critical value}
#' \item{p_value}{p value of global F test}
#' \item{wasserstein.F_stat}{Wasserstein F statistic value from the Satterthwaite method}
#' \item{chisq_df}{degrees of freedom of the null chi-square distribution}
#' @examples
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#'
#' res <- wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, Ytype = 'density', Sup = dSup)
#' summary(res)
summary.WRI <- function(WRIobject, digits = 3) {

        ans = list()
        ans$call = WRIobject$call
        ### ================== compute Wasserstain R^2 ===================== ###
        ans$r.square = round(wass_R2(WRIobject), digits = digits)

        ### ==================   compute global F test ===================== ###
        globalFres = globalFtest(WRIobject, alpha = 0.05, permutation = FALSE, numPermu = 200, bootstrap = FALSE, numBoot = 200)
        ans$global_F_pvalue = globalFres$summary_df$p_value[2]
        ans$global_wasserstein_F_stat = round(globalFres$wasserstein.F_stat, digits = digits)
        ans$global_wasserstein_F_df =  round(globalFres$chisq_df, digits = digits)

        ### ==================   compute partial F test ==================== ###
        p = ncol(WRIobject$xfit)
        coef_table = sapply(1:p, function(i) {
                reduced_model = short_wass_regress(Xfit_df = WRIobject$xfit[, -i], smoothY = WRIobject$Ysmooth)
                partialFres = partialFtest(WRIobject, reduced_res = reduced_model, alpha = 0.05)
                return(c(partialFres$statistic[1], partialFres$p_value))
        })
        coef_table = t(coef_table)
        rownames(coef_table) = WRIobject$predictor_names
        colnames(coef_table) = c('F-stat ', ' p-value(truncated) ', ' p-value(satterthwaite)')
        ans$partial_F_table = round(coef_table, digits = digits)

        class(ans) = 'summary.WRI'
        return(ans)
}
