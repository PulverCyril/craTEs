#' split a cis-regulatory susceptibility matrix N into functional and non functional fractions of TE subfamilies
#' 
#' starting from a specific N matrix that was already filtered for specific genes (e.g defined through `preprocess_E_N_for_activities`), split subfamilies into functional vs non_functional fractions according to a source functional N matrix defined for all genes. Subfamilies are only split if both their functional and non functional sum of cis-regulatory weights are above the threshold.
#'
#' @param N pre-filtered input N matrix
#' @param N_functional unfiltered functional matrix corresponding to N (e.g same bandwidth parameter L)
#' @param threshold minimum sum of cis-regulatory weights over genes at both the non-funtional and funtional fractions for a subfamily to be split
#' @param subfams_protected will be split irrespective of `threshold`
#' @param suffix character string denoting functional subfamilies (default is "_functional")
#' @return functional equivalent of the input matrix `N`
#' @export
split_N_to_functional <- function(N,
                                  N_functional,
                                  threshold = 100,
                                  subfams_protected = c(),
				  suffix = "_functional") {
    # subsetting N_functional to match the genes selected in N due to preprocessing
    N_functional = N_functional[rownames(N), ]
    
    subfams_kept = colnames(N)
    subfams_kept_functional = paste0(subfams_kept, suffix)
    subfams_kept_functional = subfams_kept_functional[subfams_kept_functional %in% colnames(N_functional)]
    subfams_exact_split = c(subfams_kept, subfams_kept_functional)

    t = 100 # threshold for minimum number of counts
    cols_N_functional_kept = c()
    cols_N_standard_kept = c()
    for (s in subfams_kept) {
        # Are there actually functional integrants ?
        s_functional = paste0(s, suffix)
        if(s_functional %in% colnames(N_functional)) {
            if (s_functional %in% subfams_kept_functional) {
                # should we keep the functional and non functional fractions ?

                # we keep it if the subfam is "protected" or
                # if both fractions sum up to at least the threshold

                protected = (s %in% subfams_protected)
                functional_split_ok = (sum(N_functional[, s]) >= t) &
                    (sum(N_functional[, s_functional]) >= t)

                if (protected | functional_split_ok) {
                    cols_N_functional_kept = c(cols_N_functional_kept, s, s_functional)
                }
                else {
                    cols_N_standard_kept = c(cols_N_standard_kept, s)
                }
            }
        }
        else {
            cols_N_standard_kept = c(cols_N_standard_kept, s)
        }
    }
    return(cbind(N[, cols_N_standard_kept],
                              N_functional[, cols_N_functional_kept]))
}
