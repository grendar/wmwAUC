#' Plot Method for wmw_test Objects
#'
#' Creates empirical ROC curve plot with test results (p-value, eAUC with confidence 
#' interval) displayed in subtitle. If `ci_method = 'boot'` was used in `wmw_test()`, 
#' the plot includes confidence bands for the ROC curve constructed using the same 
#' bootstrap resamples used for the AUC confidence interval.
#' 
#' When `special_case = TRUE` was used in `wmw_test()`, an additional boxplot with 
#' swarmplot overlay is created, showing the eAUC as effect size estimate with 
#' confidence interval in the subtitle (demonstrating the dual interpretation of 
#' eAUC in the location-shift case).
#' 
#' @param x Object of class 'wmw_test' returned by `wmw_test()`
#' @param combine_plots Logical, whether to return combined plot using patchwork 
#'   (TRUE) or list of individual plots (FALSE). Only relevant when special_case = TRUE
#' @param ... Additional arguments (not currently used)
#'
#' @return No return value, called for side effects. Creates a plot visualizing  
#' the Wilcoxon-Mann-Whitney test results including distributions,  
#' test statistic, and confidence information.
#'
#'
#' @importFrom stats as.formula
#' 
#' @export
plot.wmw_test <- function(x, combine_plots = TRUE, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

   
  ##############################################################################
  #
  # ROC
  #
  ##############################################################################
  #
  roc_plot = plot_roc(x$roc_object)
  #
  # Add AUC annotation
  auc_text <- paste0(", eAUC = ", round(x$auc, 3))
  if (!is.null(x$auc.conf.int)) {
    auc_text <- paste0('WMW p-value ',
                       ifelse(x$p.value < 0.001, "< 0.001", paste0('= ', formatC(x$p.value, format = "f", digits = 3))),
                auc_text, " (95% CI: [", round(x$auc.conf.int[1], 3),
                       ", ", round(x$auc.conf.int[2], 3), "])")
  }
  roc_plot <- roc_plot + ggplot2::labs(title = 'eROC', subtitle = auc_text)
  #
  
  ##############################################################################
  #
  # special case
  #
  ##############################################################################
  #
  if (x$special_case == TRUE) {

    ##############################################################################
    #
    # Boxplot with beeswarm
    #
    ##############################################################################
    #
    da = data.frame(response = c(x$x_vals, x$y_vals),
                    group = x$groups)
    da$group = as.factor(da$group)
    da$group = relevel(da$group, ref = x$group_ref_level)
    #
    # quadruplot
    fml = as.formula('response ~ group')
    qp = quadruplot(fml, da, 
                    ref_level = x$group_ref_level,
                    test = 'ks', # 
                    combine_plots = FALSE, 
                    show_legend = TRUE)
    box_plot = qp$ggb
    #
    # subtitle
    #wmw_test_text <- paste0('Wilcoxon-Mann-Whitney test p-value ',
    #                        ifelse(x$p.value < 0.001, "< 0.001", paste0('= ', formatC(x$p.value, format = "f", digits = 3))), '\n')
    #
    wmw_test_text = paste0('WMW p-value ',
                           ifelse(x$p.value < 0.001, "< 0.001", paste0('= ', formatC(x$p.value, format = "f", digits = 3))))
    auc_text <- paste0(". Effect size: eAUC = ", round(x$auc, 3))
    if (!is.null(x$auc.conf.int)) {
      auc_text <- paste0(auc_text, " (95% CI: [", round(x$auc.conf.int[1], 3),
                         ", ", round(x$auc.conf.int[2], 3), "])")
    }
    subtitle_text = paste0(wmw_test_text, auc_text)
    #
    box_plot <- box_plot + ggplot2::labs(subtitle = subtitle_text)


    ############################################################################
    #
    # Return based on combine_plots preference
    #
    if (combine_plots) {
      # Combine all three plots in 1x2 layout
      combined_plot <- roc_plot + box_plot
      return(combined_plot)
    } else {
      return(list(roc_plot = roc_plot, box_plot = box_plot))
    }
    #
  } else if (x$special_case == FALSE) {
    #
    # Return based on combine_plots preference
    if (combine_plots) {
      return(roc_plot)
    } else {
      return(list(roc = roc_plot))
    }
    #
  }  
  #
}


