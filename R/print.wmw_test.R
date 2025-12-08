#' Print Method for wmw_test Objects
#' 
#' Prints summary of Wilcoxon-Mann-Whitney discrimination test results.
#'
#' @param x Object of class 'wmw_test' returned by `wmw_test()`
#' @param digits Integer, number of digits to display for numeric results (default: 4)
#' @param ... Additional arguments (not currently used)
#'
#' @return Invisibly returns the input object \code{x} (of class "wmw_test"). 
#' Called primarily for side effects to print a formatted summary 
#' of the Wilcoxon-Mann-Whitney test results to the console.
#'
#' @export
print.wmw_test <- function(x, digits = 3, ...) {

  cat("\n")
  cat("        Wilcoxon-Mann-Whitney Test of No Group Discrimination\n\n")
  
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
  
  
  if (x$special_case == TRUE) {
    #
    # Show appropriate alternative hypothesis for location mode
    if (x$special_case == TRUE) {
      cat("\nLocation-shift analysis (under F1(x) = F2(x - delta)):\n")
      cat("alternative hypothesis for location:", x$alternative, "\n")
    }
    
    if (!is.null(x$pseudomedian)) {
      cat("Hodges-Lehmann median of all pairwise distances:\n")
      cat(" ", formatC(x$pseudomedian, format = "f", digits = digits),
          " [location effect size: eAUC = ", formatC(x$auc, format = "f", digits = digits),
          "]\n", sep = "")
    }
    
    # Location-specific results 
    if (!is.null(x$pseudomedian_conf_int)) {
      cat("95 percent confidence interval for median of all pairwise distances:\n")
      cat(" ", formatC(x$pseudomedian_conf_int[1], format = "f", digits = digits), " ",
          formatC(x$pseudomedian_conf_int[2], format = "f", digits = digits),
          "\n", sep = "")
    }
    #
  }

  cat("\n")
  invisible(x)
}
