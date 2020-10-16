#' print the summary of WRI object
#' @param summaryWRI a 'summary.WRI' object
#' @param digits integer indicating the number of decimal places in \code{round} function. (default: 3)
#' @export
print.summary.WRI <- function(summaryWRI){

        cat('Call:\n')
        print(summaryWRI$call)

        cat('\nPartial F test for individual effects:\n\n')
        printCoefmat(summaryWRI$partial_F_table)

        cat(paste0('\nWasserstein R-squared: ', summaryWRI$r.square, '\n'))

        cat('F-statistic (by Satterthwaite method): ')

        cat(paste0(summaryWRI$global_wasserstein_F_stat, ' on ', summaryWRI$global_wasserstein_F_df,
                   ' DF,', ' p-value: ', formatC(summaryWRI$global_F_pvalue, format = "e", digits = 3), '\n'))

}
