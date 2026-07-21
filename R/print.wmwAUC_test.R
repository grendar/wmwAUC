#' Print Method for wmwAUC_test Objects
#'
#' Prints summary of wmwAUC test results.
#'
#' @param x Object of class 'wmwAUC_test' returned by `wmwAUC_test()`
#' @param digits Integer, number of digits to display for numeric results (default: 3)
#' @param ... Additional arguments (not currently used)
#'
#' @return Invisibly returns the input object \code{x} (of class "wmwAUC_test").
#' Called primarily for side effects to print a formatted summary
#' of the wmwAUC test results to the console.
#'
#' @export
print.wmwAUC_test <- function(x, digits = 3, ...) {
  
  cat("\n")
  cat("        wmwAUC Test of No Group Discrimination (H0: AUC = 1/2)\n\n")
  
  cat("data: ", x$data_name, " (n1 = ", as.numeric(x$n)[1], ", n2 = ",  as.numeric(x$n)[2],
      ")\n", sep = "")
  i_ref_level <- which(x$group_levels == x$group_ref_level)
  i_nonref_level <- which(x$group_levels != x$group_ref_level)
  cat('groups: ', x$group_levels[i_nonref_level], ' vs ',
      x$group_levels[i_ref_level], ' (reference)\n', sep = '')
  cat("U = ", x$U_statistic, ", eAUC = ", formatC(x$auc, format = "f", digits = digits),
      ", p-value = ", formatC(x$p_value, format = "f", digits = 6),
      ", method = ", x$pvalue_method,
      "\n", sep = "")
  cat("alternative hypothesis for AUC:", x$alternative, "\n")
  #
  cat("95 percent confidence interval for AUC (",  x$ci_method, "): \n", sep = '')
  cat(" ", formatC(x$auc_conf_int[1], format = "f", digits = digits), " ",
      formatC(x$auc_conf_int[2], format = "f", digits = digits),
      "\n", sep = "")
  
  
  if (isTRUE(x$pseudomedian_requested)) {
    #
    if (!is.null(x$pseudomedian)) {
      cat("\nAUC-equalizing shift (pseudomedian):\n")
      cat(" ", formatC(x$pseudomedian, format = "f", digits = digits),
          " [solves P(X < Y + delta) = 0.5]\n", sep = "")
    }
    
    if (!is.null(x$pseudomedian_conf_int)) {
      cat("95 percent confidence interval for the pseudomedian:\n")
      cat(" ", formatC(x$pseudomedian_conf_int[1], format = "f", digits = digits), " ",
          formatC(x$pseudomedian_conf_int[2], format = "f", digits = digits),
          "\n", sep = "")
    }
    #
  }
  
  cat("\n")
  invisible(x)
}