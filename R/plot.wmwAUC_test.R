#' Plot Method for wmwAUC_test Objects
#'
#' @description Creates a combined plot ("diplot"): a boxplot with beeswarm
#' overlay of the raw group values, alongside the empirical ROC curve (eROC)
#' with confidence band (when available) and AUC confidence interval.
#' Optionally adds a third panel with overlaid group density curves.
#'
#' @param x Object of class 'wmwAUC_test' returned by \code{\link{wmwAUC_test}}.
#' @param show_density Logical; if TRUE, adds a third panel with overlaid
#'   kernel density estimates for the two groups. Default FALSE.
#' @param combine_plots Logical; if TRUE (default), returns a single combined
#'   plot via \pkg{patchwork}; if FALSE, returns a named list of the
#'   individual \pkg{ggplot2} objects instead.
#' @param ... Additional arguments (not currently used).
#'
#' @return If \code{combine_plots = TRUE}, a combined \pkg{patchwork} object.
#' If \code{combine_plots = FALSE}, a named list with elements \code{box_plot},
#' \code{roc_plot}, and (if \code{show_density = TRUE}) \code{density_plot}.
#'
#' @export
plot.wmwAUC_test <- function(x, show_density = FALSE, combine_plots = TRUE, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Package 'ggbeeswarm' is required for plotting. Please install it.")
  }
  
  ##############################################################################
  #
  # data frame for box/swarm and (optional) density panels
  #
  ##############################################################################
  case_label <- x$group_levels[x$group_levels != x$group_ref_level]
  ref_label  <- x$group_ref_level
  
  da <- data.frame(
    response = c(x$x_vals, x$y_vals),
    group = factor(
      c(rep(case_label, length(x$x_vals)), rep(ref_label, length(x$y_vals))),
      levels = c(ref_label, case_label)
    )
  )
  
  ##############################################################################
  #
  # shared subtitle text: p-value, eAUC (+CI), and pseudomedian (+CI) if requested
  #
  ##############################################################################
  p_text <- ifelse(x$p_value < 0.001, "< 0.001",
                   paste0("= ", formatC(x$p_value, format = "f", digits = 3)))
  auc_text <- paste0("wmwAUC p-value ", p_text, ", eAUC = ", round(x$auc, 3))
  if (!is.null(x$auc_conf_int)) {
    auc_text <- paste0(auc_text, " (95% CI: [", round(x$auc_conf_int[1], 3),
                       ", ", round(x$auc_conf_int[2], 3), "])")
  }
  if (isTRUE(x$pseudomedian_requested) && !is.null(x$pseudomedian)) {
    pm_text <- paste0("; pseudomedian = ", round(x$pseudomedian, 3))
    if (!is.null(x$pseudomedian_conf_int)) {
      pm_text <- paste0(pm_text, " (", round(100 * x$pseudomedian_conf_level), "% CI: [",
                        round(x$pseudomedian_conf_int[1], 3), ", ",
                        round(x$pseudomedian_conf_int[2], 3), "])")
    }
    auc_text <- paste0(auc_text, pm_text)
  }
  
  ##############################################################################
  #
  # box + beeswarm panel
  #
  ##############################################################################
  box_plot <- ggplot2::ggplot(da, ggplot2::aes(x = group, y = response)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    ggbeeswarm::geom_beeswarm(ggplot2::aes(color = group), alpha = 0.7, size = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "Group values", subtitle = auc_text, x = NULL, y = NULL)
  #
  lancet_colors <- ggsci::pal_lancet()(9)
  box_plot <- box_plot + ggplot2::scale_color_manual(values = lancet_colors)
  #
  
  ##############################################################################
  #
  # ROC panel: reuses plot_roc(), which reads auc/auc_ci directly off
  # x$roc_object (see plot_roc.R)
  #
  ##############################################################################
  roc_plot <- plot_roc(x$roc_object) + ggplot2::labs(title = "eROC")
  
  panels <- list(box_plot = box_plot, roc_plot = roc_plot)
  
  ##############################################################################
  #
  # optional density panel
  #
  ##############################################################################
  if (show_density) {
    density_plot <- ggplot2::ggplot(da, ggplot2::aes(x = response, fill = group, color = group)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Group densities", x = NULL, y = "density")
    #
    density_plot <- density_plot + 
      ggplot2::scale_color_manual(values = lancet_colors) +
      ggplot2::scale_fill_manual(values = lancet_colors) +
      ggplot2::theme(legend.position = "none")
    #
    panels <- list(box_plot = box_plot, density_plot = density_plot, roc_plot = roc_plot)
    #
    # panels$density_plot <- density_plot
    #
  }
  
  ##############################################################################
  #
  # combine or return as a list
  #
  ##############################################################################
  if (combine_plots) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' is required for combine_plots = TRUE. ",
           "Please install it, or set combine_plots = FALSE.")
    }
    combined <- Reduce(`+`, panels)
    return(combined)
  } else {
    return(panels)
  }
}