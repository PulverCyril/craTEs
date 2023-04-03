# +
#' Parallelized RMSE estimation using cross-validation
#'
#' Computes the RMSE in the context of a parallelizable cross-validation scheme, 
#' where the gene space is split according to the number of folds k.
#' k-1 folds are used for training, and the remaining fold is used as a validation fold to estimate the RMSE
#'
#' @param preprocessed data produced by preprocess_E_N_for_activities()
#' @param control_group vector containing the names of the control samples (column names in preprocessed$E_centered)
#' @param treatment_group vector containing the sames of the treatment samples (column names in preprocessed$E_centered)
#' @param N_list list of the different N matrices. Rows must already correspond to those of preprocessed$E_centered
#' @param N_metadata dataframe with a number of rows equal to the number of elements in N_list
#' @param folds vector of the same length as the number of observations in E_centered, with integers from 1:nFolds, each defining a fold
#' @param f fold index for which to estimate the RMSE
#' @return the dataframe `N_metadata`, modified with an additional column containing the estimated RMSEs
#' @export

estimate_rmse_crossval <- function(preprocessed, control_group, treatment_group,
                                  N_list, N_metadata,
                                  folds, f) {

    idx_training = folds==f
    idx_test = !idx_training
    # initializing for temporary results
    rmse = N_metadata
    rmse$fold = as.integer(0)
    rmse$rmse = 0       
    rmse['fold'] = f
    for(i in 1:nrow(N_metadata)) {
        rmse[[i, 'rmse']] = estimateRMSE(preprocessed$E_centered, N=N_list[[i]], treatment_group = treatment_group, control_group=control_group, idx_training = idx_training, idx_val = idx_test)
    }
    return(rmse)
}
