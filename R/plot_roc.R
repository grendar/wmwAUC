#' ROC plot with confidence band -- internal function
#'
#' @param x Object of class `roc_ci` returned by `roc_with_ci()`
#' @param ... not used
#' 
#' @keywords internal
#' @export
plot_roc = function(x, ...) {
  #
  #
  roc_df <- x$roc_df
  roc_band <- x$roc_band
  auc <- x$roc_object$auc
  auc_ci <- x$roc_object$auc_ci
  #
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = roc_df,
                       ggplot2::aes(x = fpr, y = tpr), linewidth = 1.0, color = "blue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "ROC Curve",
      subtitle = sprintf(
        "AUC = %.3f (%.3f, %.3f)",
        auc,
        auc_ci[1],
        auc_ci[2]
      ),
      x = "False Positive Rate",
      y = "True Positive Rate"
    )

  if (!is.null(roc_band)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = roc_band,
        ggplot2::aes(x = fpr, ymin = tpr_lo, ymax = tpr_hi),
        alpha = 0.2,
        fill = "blue"
      )
  }
  #
  return(p)
  #
}

