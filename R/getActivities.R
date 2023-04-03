#' Compute the max glmnet lambda penalization across several samples
#'
#' Computes the max value of lambda such that all coefficients
#' remain 0 when running a penalized linear regression of Y onto X with the
#' glmnet package. Based on the coordinate descent udpate formula. 
#' 
#' @param Y observation/response vector
#' @param X predictor matrix
#' @param alpha mixing parameter for the penalty term (0 -> ridge, 1 -> lasso)
#' @return the smallest value of lambda for which all coefficients are set to 0 across all samples
get_lambda_max <- function(Y, X, alpha) {
    # compute the standard deviation with the 1/N formula (not the 1/(N-1)) formula:
    mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))

    stopifnot(length(Y)==dim(X)[1])

    return ((1/(length(Y)*alpha))*max(abs(apply(scale(X, center = T, scale = apply(X, 2, mysd))*(Y), MARGIN=2, sum))))
}

#' Estimate activities (transcription factors, TE subfamilies) using the elastic net regularization 
#' 
#' Estimates activities by regressing each
#' column of the expression matrix (1 col per sample) onto the per promoter
#' motif site counts. Multiple linear regression with the
#' elastic net regularization is performed by the glmnet package.
#' We aim to find the value of lambda which minimizes the error
#' over all samples. Therefore, we must provide the same grid of lambda values
#' to each regression problem (one per sample). The max value of lambda is 
#' defined as the maximum value of lambda_max across samples.
#' lambda_max for a single sample is given by the formula:
#' lambda_max = 1/(nobs*alpha)*max_s(abs(dot product(Z(x_s), y_s))) where Z
#' is a function that centers AND standardizes its input.
#'
#' @param expr row-centered TPM values in a promoter/gene (rows) x sample (columns) matrix.
#' @param X per promoter/gene (rows) motif (columns) counts. Column-centered. Also referred to as N.
#' @param k number of folds for cross-validation
#' @param alpha trades off between ridge (alpha = 0) and lasso (alpha = 1) penalties
#' @param ncores number of cores for multiprocessing (1 -> no multiprocessing)
#' @param nlambda number of lambda values in the grid from lambda max to lambda min
#' @return a list object with the matrix of estimated activities, and diagnostics (e.g r2 etc)
#' @export
getActivities <- function(expr, X, alpha = 0.1, k = 5, nlambda=100,
                           ncores = 1, log_id = '') {
    print("Going parallel !")
    cl <- makeCluster(ncores, type="FORK", outfile=paste(c("./getActivities_log", log_id, ".txt"), collapse=''))

    n = ncol(expr)


    cat("Computing lambda max over ", n, " samples\n")
    lambda_max_list = parApply(cl, expr, MARGIN = 2, FUN = get_lambda_max, X, alpha = 0.1)
    # stopping the cluster to communicate the grid of lambda to the workers after computing it
    stopCluster(cl)

    cat("Building the lambda grid")
    lambda_max = max(lambda_max_list)
    log_lambda_grid = seq(from=log(lambda_max), to = log(lambda_max*1e-4), length.out = nlambda)
    lambda_values = exp(log_lambda_grid)

    cat("Regressing activities of TFs for ", n, " samples\n")
    # regularized least squares on single expression vector
    f <- function(y) {
        print("started new sample")

        return(glmnet::cv.glmnet(x = X, y = y, alpha = alpha, 
                           nfolds = k, lambda = lambda_values, type.measure='mse'))
    }

    cl <- makeCluster(ncores, type="FORK", outfile=paste(c("./getActivities_log", log_id, ".txt"), collapse=''))

    res = list()

    res <- parApply(cl, expr, MARGIN = 2, FUN = f)

    # shutting down parallel clusters
    stopCluster(cl)

    cat("Finding the optimal lambda value across samples...\n")

    cat("\t... 1) Computing the per lambda MSE across samples\n")
    # mse estimates: rows as lambda values, and columns as samples
    errors = do.call('cbind', Biobase::subListExtract(res, 'cvm'))

    # average over mse estimates across samples, one value per lambda in the grid 
    lambda_mse = rowSums(errors)/ncol(errors)
    # keeping sd(mse) to plot lambda vs mse plot later
    lambda_mse_sd = apply(errors, MARGIN=1, sd)
    
    #retrieving the value of lambda that yields the minimum across sample mse
    lambda_min_idx = which.min(lambda_mse)
    lambda_min = lambda_values[lambda_min_idx]

    cat("\t... 2) retrieving the coefficients, crossval MSE, r_squared for the optimal lambda value\n")
    coefs = lapply(X = Biobase::subListExtract(res, 'glmnet.fit'), FUN = glmnet::predict.glmnet, type='coefficients', s=lambda_min)

    # r_squared at the overall best lambda
    temp = Biobase::subListExtract(res, 'glmnet.fit')
    r_squared_temp = Biobase::subListExtract(temp, 'dev.ratio')
    r_squared = lapply(r_squared_temp, FUN = function(x) x[lambda_min_idx])
    rm(list = 'temp')

    # r_squared max for each sample (to see how much we lose):
    r_squared_max = lapply(r_squared_temp, FUN = max)

    # best lambda for each sample (for the record)
    best_lambda = Biobase::subListExtract(res, 'lambda.min')

    # reshaping results into a five-list object: ans$bestlambda, ans$r2,
    # ans$test_err_est, ans$test_err_se and ans$activities

    ans <- list()
    ans$r2 <- unlist(r_squared)
    ans$r2_max <- unlist(r_squared_max)
    ans$lambdas <- lambda_values 
    ans$lambda_min_overall <- lambda_min
    ans$lambdas_best <- unlist(best_lambda)
    ans$lambda_mse <- unlist(lambda_mse)
    ans$lambda_mse_sd <- unlist(lambda_mse_sd)
    
    # Building the TF activities (rows) by samples (columns) table
    ans$activities <- do.call(cbind, coefs)
    colnames(ans$activities) <- names(coefs)
    return(ans)
}

#' Scatter plot of TF activities per sample 
#'
#' @param act TF x samples matrix, with row and col names set
#' @param TF string for TF of interest, e.g "MYC"
#' @param metadata samples x metadata dataframe. Metadata should be set as factors
#' @param sort_by string corresponding to a column of metadata
#' @param decreasing bool to sort samples by metadata[sort_by] in decr. order.
#' @param palette_type string, one of c("sequential", "diverging", "categorical")
#' @param palette_add_black bool to enable setting a color to black (useful for baselines)
#' @param palette_black_idx category to plot in black
#' @param title title if it is to be different than the TF name
#' @param act_lim limits of activities to plot, useful to get comparable y axes
#' @export
plotActivities <- function(act, TF, metadata, sort_by, decreasing=F, 
                           palette_type='sequential', 
                           palette_add_black = T, palette_black_idx=1,
                           title=NULL, act_lim=NULL) {
    # to reset after fun
    old_pal = palette()
    old_par = par(no.readonly=T)

    # matching TF name to idx
    tf_idx = which(rownames(act)==TF)

    # sorting RNA seq samples by selected metadata
    sample_order_idx = seq(1, dim(act)[2])
    if(!is.null(sort_by)){
        sample_order_idx = order(metadata[sort_by], decreasing=decreasing)
    }

    # generating a palette
    new_pal = palette("default")
    # sorting by categorical factor -> n colors in the palette = n factors
    if(is.factor(metadata[[sort_by]])) {
        if (palette_type=='sequential'){
                new_pal = RColorBrewer::brewer.pal(n = length(levels(metadata[[sort_by]])), name = "OrRd")
            } else if (palette_type=='diverging') {
            new_pal = RColorBrewer::brewer.pal(n = length(levels(metadata[[sort_by]])), name = "Spectral")
            } else if (palette_type=='qualitative') {
            new_pal = RColorBrewer::brewer.pal(n = length(levels(metadata[[sort_by]])), name = "Set1")
        }
    }

    # In case of ref level chosen as black
    if (palette_add_black){
        new_pal[palette_black_idx] = "black"
    }
    palette(new_pal)

    
    # enlarging margins for legend out of plot, forcing square plot area
    par(mai=c(0.2, 0.8, 0.4, 1.2), pty='s')
        

    # title
    main_title=title
    if (is.null(title)) {
        main_title=TF
    }

    plot(x = seq(from = 1, to = ncol(act)),
         y = act[tf_idx, sample_order_idx],
         col = metadata[[sort_by]][sample_order_idx],
         main = main_title,
         xlab = paste("Samples (ordered by ", sort_by, ")", sep=""),
         ylab = "Activity (deviation from average sample activity)",
         ylim=act_lim,
         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
    legend("right", 
           legend = levels(metadata[[sort_by]]), 
          col = new_pal, pch = 'o', inset=c(-0.26, 0), 
            box.lwd = 0, bty = "n", xpd=T, cex=1.2)
    abline(a = 0, b = 0, col = "blue", lty=2, lwd = 2)

    # resetting palette and par
    #palette(old_pal)
    #par(old_par)
}
