#' Returns the pvalue corresponding to the f-statistic of the complete versus the null model (intercept only)
#'
#' @param res result returned by getSignifActivities
#' @return p-vaue for overall significance test of the linear model
#' @export
pval_from_fstat <- function(res) {
    pval = pf(q = res$fstatistic['value'], df1 = res$fstatistic['numdf'], df2 = res$fstatistic['dendf'], lower.tail = F)
    return(pval)
}
