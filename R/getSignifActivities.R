#' Estimates and returns the activities of TE subfamilies by multiple linear regression
#'
#' Each activity coefficient is tested agains the null hypothesis that it is equal to zero
#' A Benjamini-Hochberg procedure to control the FDR is employed.
#'
#' @param E_centered matrix of gene expression values with genes as rows and samples as columns
#' @param N predictor matrix with genes as rows and TE subfamilies as columns
#' @param treatment_group vector with the name of treated samples (correspond to columns in E_centered)
#' @param control_group same for control samples
#' @param exclude_genes vector with the name of genes to be exluded from E and N, if their expression is assumed to be driven by TE-independent means (e.g overexpression / knockdown / knockout)
#' @return list object with estimated activities and useful diagnostics (e.g R2)
#' @export
getSignifActivities <- function(E_centered, N, treatment_group, control_group, exclude_genes = NULL) {
    print('Starting new sample')

    # defensive programming
    stopifnot(dim(E_centered)[1]==dim(N)[1])

    stopifnot(all(sapply(treatment_group, function(x) x %in% colnames(E_centered))))
    stopifnot(all(sapply(control_group, function(x) x %in% colnames(E_centered))))

    stopifnot(all(sapply(control_group, function(x) ! x %in% treatment_group)))

    # excluding artificially manipulated genes
    genes_excluded = c()
    if(! is.null(exclude_genes)) {
	    genes_excluded = exclude_genes[exclude_genes %in% rownames(E_centered)]
	    if (length(genes_excluded) > 0) {
		    to_delete = which(rownames(E_centered) %in% genes_excluded)
		    E_centered = E_centered[-to_delete, ]
		    N = N[-to_delete, ]
	    }
    }
    stopifnot(all(rownames(E_centered) == rownames(N)))

    # generating the difference in expression vector

    control = c(sapply(control_group, (function(x) E_centered[, x])))
    treatment = c(sapply(treatment_group, (function(x) E_centered[, x])))
    delta = treatment-control
               
    # generating the augmented N matrix to match the dimensions of delta
    n_replicates = length(control_group)
    N_augmented = N

    i = 1
    while (i < n_replicates) {
        N_augmented = rbind(N_augmented, N)  
        i = i+1
    }

                       
    # assembling the design matrix, regression delta as a function of all columns in N_augmented
    design_matrix = as.data.frame(cbind(delta, N_augmented))
    fit = lm(delta ~ ., data = design_matrix)
    coefs = as.data.frame(summary(fit)$coefficients)

    # FDR adjustment

    # ignoring the intercept (-1)
    p_adj = p.adjust(coefs[-1, 4], method='BH')

    # adding an NA to match the dimensions of `coefs` and appending
    p_adj = append(p_adj, NA, 0)
    coefs$p_adj = p_adj

    res = list()
    res$coefs = coefs
    res$control_group = control_group
    res$treatment_group = treatment_group
    res$r.squared = summary(fit)$r.squared
    res$ajd.r.squared = summary(fit)$adj.r.squared
    res$aliased = summary(fit)$aliased
    res$fstatistic = summary(fit)$fstatistic                                
    res$rss = deviance(fit)
    res$p = ncol(N)
    res$N = nrow(E_centered)
    res$genes_excluded = genes_excluded
    return(res)
}
