#' Volcano plot with TE subfamily activities on the x axis, activity coefficient significance on the y axis
#'
#' Black dots are significant activities, grey dots are non significant activities. 
#' The 95% confidence interval are plotted for all significant activities as grey bars.
#'
#' @param res result file returned by getSignifActivities.R
#' @param alpha BH-adjusted threshold for significance (black circles)
#' @param n_label number of subfamilies to label on the plot, starting from the most significant subfamily
#' @param main title of the plot
#' @param n_subplots when plotting several plots on the same figure (e.g par(mar=c(2, 2))), we only reset par() after the last subplot is drawn
#' @param idx belongs to the interval [1, n_subplots] and triggers the reset of par() when it equals n_subplots
#' @param plot_signif_thresh whether to explicitely plot the significance threshold as a dotted line
#' @export 
plot_signif_subfams <- function(res, alpha, n_label=NULL, main=NULL, n_subplots=1, idx=1, plot_signif_thresh=FALSE) { 
   
    stopifnot(alpha <= 1 & alpha >= 0)
    stopifnot(idx <= n_subplots)
    stopifnot(n_subplots >= 0)
    
    # saving current plot margins
    old_par = par(no.readonly = TRUE)
    
    if(idx == n_subplots) {
        on.exit(par(old_par))
    }
    # changing margin sizes
    par(mar=c(4.1, 4.8, 2.5, 2.5))
    
    # when doing subplots, no need for space for the x and y labels
    if(n_subplots > 1) {
        par(mar=c(3.5, 3.5, 2.5, 2.5))
    }
    
    # identifying significant coefficients (will be plotted in black, not grey)
    coefs_non_signif = res$coefs
    signif = which(res$coefs$p_adj <= alpha)
    coefs_signif = res$coefs[signif,]
    if (nrow(coefs_signif)>0) {
        coefs_non_signif = res$coefs[-signif,]
    }
    
    coef_labels = NULL
    largest_label_length = 0
    if (! is.null(n_label) & nrow(coefs_signif) > 0) {
        coef_labels = coefs_signif[head(order(coefs_signif$p_adj), n_label), ]
        largest_label_length = max(sapply(rownames(coef_labels), nchar))

    }
    
    # correcting the plot margins for the error bars and the labels to fit in
    std_offset = rep(0, nrow(coefs_signif))
    std_offset_margin = 0
    if (nrow(coefs_signif)> 0) {
        std_offset = 1.96*coefs_signif$`Std. Error`
        std_offset_margin = max(std_offset)
    }
    
    # empty plot with limits
    x_min = min(res$coefs$Estimate)
    x_max = max(res$coefs$Estimate)
    x_margin = 0
    if (! is.null(n_label)) {
        x_margin = max(std_offset_margin, largest_label_length/280)
    } else {x_margin = std_offset_margin}
    
    y_min = 0
    y_max = max(na.omit((-log10(res$coefs[, 'p_adj']))))
    if(plot_signif_thresh & y_max < -log10(alpha)) {
        y_max = -log10(alpha)
        
    }
    y_margin = 0.04*(y_max-y_min)
    
    xlabel = 'Activity'
    ylabel = '-log10(adj. p-value)'
    
    # when doing subplots, x and y labels should be added independently
    if(n_subplots > 1) {
        xlabel = ''
        ylabel = ''
    }
    plot(NULL, type='n', xlim = c(x_min - x_margin, x_max + 2*x_margin), ylim = c(y_min, y_max + y_margin),
        xlab = xlabel, ylab = ylabel, cex.lab = 2.5, cex.axis = 2.5)

    # moving main title up
    title(main, line=1, cex.main = 2)
    
    # plotting the significance threshold
    if(plot_signif_thresh) {
        abline(h = -log10(alpha), col = 'red', lty = 2, lwd = 2)
    }
    
    # plotting the 0.95 confidence interval for the estimated coefficient
    if(nrow(coefs_signif > 0)) {
        arrows(x0=coefs_signif$Estimate - std_offset, y0=-log10(coefs_signif$p_adj), x1=coefs_signif$Estimate + std_offset, y1=-log10(coefs_signif$p_adj), 
               code=3, angle=90, length=0.05,
               col="grey", lwd=1.5)
    }

    
    # non significant activities
    points(coefs_non_signif$Estimate, -log10(coefs_non_signif$p_adj), col = 'grey')
    
    # significant activities
    points(coefs_signif$Estimate, -log10(coefs_signif$p_adj), col = 'black', cex=1.3)
    

    if (! is.null(n_label) & nrow(coefs_signif) > 0) {
        # adding the name of the top n TE subfams to the plot
        text(x = coef_labels$Estimate, y=-log10(coef_labels$p_adj), labels = rownames(coef_labels), adj = c(-0.1, -0.4), cex = 1.8) 
    }
}
