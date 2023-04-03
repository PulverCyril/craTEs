#' Estimates the RMSE for gene expression predicted from the TE subfamily activities
#'
#' Estimates the differential activities of TE subfamilies by multiple regression,
#' and uses the activities estimated on the training observations to estimate the
#' root mean squared error (RMSE) on the validation observations
#' 
#' @param E_centered matrix of normalized gene expression values with genes as rows and samples as columns
#' @param N TE cis-regulatory predictor matrix with genes as rows and TE subfamilies as columns
#' @param treatment_group vector with the name of treated samples (must correspond to columns in E_centered)
#' @param control_group vector with the name of control samples (must correspond to columns in E_centered)
#' @param idx_training indexes the training observations
#' @param idx_val indexes the validation observations
#' @return the estimated RMSE for the difference in gene expression predicted from the TE activities
#' @export
estimateRMSE <- function(E_centered, N, treatment_group, control_group, idx_training, idx_val) {
    print('Starting new sample')

    # defensive programming
    stopifnot(dim(E_centered)[1]==dim(N)[1])

    stopifnot(all(sapply(treatment_group, function(x) x %in% colnames(E_centered))))
    stopifnot(all(sapply(control_group, function(x) x %in% colnames(E_centered))))

    stopifnot(all(sapply(control_group, function(x) ! x %in% treatment_group)))
    stopifnot(all(idx_training != idx_val))

    # separating the training from the test set
    E_training = E_centered[idx_training, ]
    E_test = E_centered[idx_val, ]
    
    N_training = N[idx_training, ]
    N_test = N[idx_val, ]
    
                         
    # training set
    # generating the difference in expression vector
    control = c(sapply(control_group, (function(x) E_training[, x])))
    treatment = c(sapply(treatment_group, (function(x) E_training[, x])))
    delta = treatment-control    
                                           
    # generating the augmented N matrix to match the dimensions of delta
    n_replicates = length(control_group)
    N_augmented = N_training

    i = 1
    while (i < n_replicates) {
        N_augmented = rbind(N_augmented, N_training)  
        i = i+1
    }

                       
    # assembling the design matrix, regression delta as a function of all columns in N_augmented
    design_matrix = as.data.frame(cbind(delta, N_augmented))
    fit = lm(delta ~ ., data = design_matrix)
                                           
    
    # test set
    # generating the difference in expression vector
    control = c(sapply(control_group, (function(x) E_test[, x])))
    treatment = c(sapply(treatment_group, (function(x) E_test[, x])))
    delta = treatment-control    
                                           
    # generating the augmented N matrix to match the dimensions of delta
    n_replicates = length(control_group)
    N_augmented = N_test

    i = 1
    while (i < n_replicates) {
        N_augmented = rbind(N_augmented, N_test)  
        i = i+1
    }

                       
    # assembling the design matrix, regression delta as a function of all columns in N_augmented
    design_matrix = as.data.frame(N_augmented)
    pred = predict(fit, newdata=design_matrix)

    rmse = sqrt(mean((delta-pred)^2))
    return(rmse)
}
