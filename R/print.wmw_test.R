#' Print Method for wmw_test Objects
#' 
#' Prints summary of Wilcoxon-Mann-Whitney discrimination test results.
#'
#' @param x Object of class 'wmw_test' returned by `wmw_test()`
#' @param digits Integer, number of digits to display for numeric results (default: 4)
#' @param ... Additional arguments (not currently used)
#'
#' @export
print.wmw_test <- function(x, digits = 3, ...) {

  cat("\n")
  cat("        Wilcoxon-Mann-Whitney Test of No Group Discrimination\n\n")
  
  cat("data: ", x$data.name, " (n1 = ", as.numeric(x$n)[1], ", n2 = ",  as.numeric(x$n)[2],
      ")\n", sep = "")
  i_ref_level <- which(x$group_levels == x$group_ref_level)
  i_nonref_level <- which(x$group_levels != x$group_ref_level)
  cat('groups: ', x$group_levels[i_nonref_level], ' vs ', 
      x$group_levels[i_ref_level], ' (reference)\n', sep = '')
  cat("W = ", x$statistic, ", p-value = ",
      formatC(x$p.value, format = "f", digits = 6), "\n", sep = "")
  cat("alternative hypothesis for AUC:", x$alternative, "\n")
  #
  cat("95 percent confidence interval for AUC (",  x$ci_method, "): \n", sep = '')
  cat(" ", formatC(x$auc.conf.int[1], format = "f", digits = digits), " ",
      formatC(x$auc.conf.int[2], format = "f", digits = digits),
      "\n", sep = "")
  #
  cat("empirical AUC (eAUC):\n")
  
  
  if (x$special_case == FALSE) {
    #
    cat(" ", formatC(x$auc, format = "f", digits = digits),
        "\n", sep = "")
    #
  } else if (x$special_case == TRUE) {
    #
    cat(" ", formatC(x$auc, format = "f", digits = digits),
        " [discrimination effect size]\n", sep = "")
    #
    #
    # Show appropriate alternative hypothesis for location mode
    if (x$special_case == TRUE && !is.null(x$null.value)) {
      cat("\nLocation-shift analysis (under F1(x) = F2(x - delta)):\n")
      cat("alternative hypothesis for location:", x$alternative, "\n")
    }
    
    # Location-specific results 
    if (!is.null(x$conf.int)) {
      cat("95 percent confidence interval for median of all pairwise distances:\n")
      cat(" ", formatC(x$conf.int[1], format = "f", digits = digits), " ",
          formatC(x$conf.int[2], format = "f", digits = digits),
          "\n", sep = "")
    }
    if (!is.null(x$estimate)) {
      cat("Hodges-Lehmann median of all pairwise distances:\n")
      cat(" ", formatC(x$estimate, format = "f", digits = digits),
          " [location effect size: eAUC = ", formatC(x$auc, format = "f", digits = digits),
          "]\n", sep = "")
    }
    #
  }

  cat("\n")
  invisible(x)
}
