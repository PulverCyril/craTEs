#' Volcano plot with TE subfamily activities on the x axis, activity coefficient significance on the y axis, with functional vs non functional activities for the selected subfamilies
#'
#' Black dot: activity of the entire subfamily. Red dot: activity of the functional fraction of the subfamily. Blue dot: activity of the non-functional fraction of the subfamily. 
#'
#' @param res result file returned by getSignifActivities.R
#' @param res_functional result file returned by `getSignifActivities()` using the functional version of N
#' @param subfams_to_plot subfamilies for which to plot the functional vs non functional activities
#' @param main title of the plot
#' @param n_subplots when plotting several plots on the same figure (e.g par(mar=c(2, 2))), we only reset par() after the last subplot is drawn
#' @param idx belongs to the interval [1, n_subplots] and triggers the reset of par() when it equals n_subplots
#' @param subfams_to_label subfamilies for which a label is written on the plot, next to the functional activity point (red). Defaults to `subfams_to_plot`
#' @param suffix character string denoting the functional fractions of subfamilies
#' @export 
plot_signif_subfams_functional <- function(res,
                                           res_functional,
                                           subfams_to_plot, 
                                           main=NULL,
                                           n_subplots=1,
                                           idx=1,
                                          subfams_to_label = NULL,
					  suffix = "_functional") { 
    ## Description: volcano plot with TE subfamily activities on the x axis, activity coefficient significance on the y axis
        ## para: -res: result file returned by getSignifActivities.R
                #-alpha: BH-adjusted threshold for significance (black circles)
                #-n_label: number of subfamilies to label on the plot, starting from the most significant subfamily
                #-main: title of the plot
                #-n_subplots: when plotting several plots on the same figure (e.g par(mar=c(2, 2))), we only reset par() after the last subplot is drawn
                #-idx: belongs to the interval [1, n_subplots] and triggers the reset of par() when it equals n_subplots
    
    #stopifnot(alpha <= 1 & alpha >= 0)
    stopifnot(idx <= n_subplots)
    stopifnot(n_subplots >= 0)
    
    # saving current plot margins
    old_par = par(no.readonly = TRUE)
    
    if(idx == n_subplots) {
        on.exit(par(old_par))
    }
    # changing margin sizes
    par(mar=c(4.1, 4.8, 5, 2.5))
    
    # when doing subplots, no need for space for the x and y labels
    if(n_subplots > 1) {
        par(mar=c(3.5, 3.5, 5, 2.5))
    }
    
    # for the selected subfams, we only plot three values:
    # 1) the activity of the functional fraction (from res_functional, subfam_functional) in red
    # 2) the activity of the whole subfamily (from res) in black
    # 3) the activity of the non functional fraction of the subfamily (from res_functional, subfam) in blue
    
    subfams_to_plot_functional = paste0(subfams_to_plot, suffix)

    act_functional = res_functional$coefs[subfams_to_plot_functional, ]
    act_non_functional = res_functional$coefs[subfams_to_plot, ]
    rownames(act_non_functional) = paste0(rownames(act_non_functional), paste0('non', suffix))
    
    act_standard = res$coefs[subfams_to_plot, ]
    
    # gathering all activities in a single table
    all_act = rbind(act_functional, act_non_functional, act_standard)    
    
    # defining the length of the largest label (for margin)
    largest_label_length = max(sapply(rownames(act_non_functional), nchar))
    
    
    # correcting the plot margins for the error bars and the labels to fit in
    std_offset = 1.96*all_act$`Std. Error`
    std_offset_margin = max(std_offset)
    
    # empty plot with limits
    x_min = min(all_act$Estimate)
    x_max = max(all_act$Estimate)
    
    x_margin = max(std_offset_margin, largest_label_length/200)
    
    y_min = 0
    y_max = max(na.omit((-log10(all_act[, 'p_adj']))))
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
    
    # plotting the 0.95 confidence interval for the estimated coefficient
    #arrows(x0=coefs_signif$Estimate - std_offset, y0=-log10(coefs_signif$p_adj), x1=coefs_signif$Estimate + std_offset, y1=-log10(coefs_signif$p_adj), 
    #       code=3, angle=90, length=0.05,
    #       col="grey", lwd=1.5)
    
    # plotting the 0.95 confidence interval for the estimated coefficient
    # does it really mean something to plot the CI for different regression ?
    
    #arrows(x0=all_act$Estimate - std_offset, y0=-log10(all_act$p_adj), x1=all_act$Estimate + std_offset, y1=-log10(all_act$p_adj), 
    #       code=3, angle=90, length=0.05,
    #       col="grey", lwd=1.5)
    
    # drawing segments between functional and standard activities
    segments(x0 = act_functional$Estimate, y0 = -log10(act_functional$p_adj), 
             x1 = act_standard$Estimate, -log10(act_standard$p_adj))
    
    # drawing segments between standard and non functional
    segments(x0 = act_standard$Estimate, y0 = -log10(act_standard$p_adj), 
             x1 = act_non_functional$Estimate, -log10(act_non_functional$p_adj))
    
    # functional_activities
    points(act_functional$Estimate, -log10(act_functional$p_adj), col = 'red', cex=1.3)
    
    # non_functional_activities
    points(act_non_functional$Estimate, -log10(act_non_functional$p_adj), col = 'blue')

    # standard activities
    points(act_standard$Estimate, -log10(act_standard$p_adj), col = 'black', cex=1.3)
    
    # adding the line of significance
    abline(h = -log10(0.05), lty = 'dotted')
    
    # labeling, only for subfams_to_legend 
    if(is.null(subfams_to_label)) {
        text(x = act_functional$Estimate, y=-log10(act_functional$p_adj), labels = gsub(suffix, '', rownames(act_functional)), adj = c(-0.1, -0.4), cex = 1.4) 
    } else {
        
        subfam_selection = rownames(act_functional)[rownames(act_functional) %in% paste0(subfams_to_label, suffix)]
        act_functional_subfam_selection = act_functional[subfam_selection, ]
        text(x = act_functional_subfam_selection$Estimate, y=-log10(act_functional_subfam_selection$p_adj), labels = gsub(suffix, '', rownames(act_functional_subfam_selection)), adj = c(-0.1, -0.4), cex = 1.8) 
    }
        
}
