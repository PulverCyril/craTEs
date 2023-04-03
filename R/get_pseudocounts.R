#'  Adds a pseudocount to the count table
#'
#' Ccomputes a per-sample pseudocount based on a quantile of non-zero values (by default 5%) 
#' and adds them to count table (genes are rows, samples are columns).
#' Useful to replace zero values before logging
#'
#' @param count_table gene count table with rows as genes, columns as samples
#' @param q quantile of non-zero values used as a pseudocount in each sample (default: q = 0.05)
#' @return count table with added pseudocounts

get_pseudocounts <- function(count_table, q=0.05) {
        stopifnot(is.matrix(count_table))
    
        if(typeof(count_table)!='integer') {
            warning("The count table provided to get_pseudocounts() is not an integer matrix, use at your own risks")
        }

        pseudocounts = apply(X = count_table, MARGIN = 2, FUN = function(x) quantile(x[x!=0], q))
        
        # Broadcasting only works correctly by rows in R, so when broadcasting by columns
        # one has to transpose to broadcast by rows (instead of cols), do the operation
        # and then tranpose again to get back to the original conformation. 
        return(t(t(count_table) + pseudocounts))
}
