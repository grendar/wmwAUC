#' Test of equality of distributions from twosamples library applied to median-centered data  
#'
#' Applies a specified test from `twosamples` library to median-centered data.
#'
#' @param x vector
#' @param y vector
#' @param test one of c("ks", "kuiper", "cvm", "ad", "wass", "dts")
#'
#' @references
#' For more details see the [Two Sample Test Package Website](https://twosampletest.com/)
#' 
#' Dowd, C. (2020). A new ECDF two-sample test statistic. arXiv preprint arXiv:2007.01360.
#'
#' @keywords internal
#' @export
test_shift_equivalence <- function(x, y, test = "ks") {
  #
  # Center by medians
  x_centered <- x - median(x, na.rm = T)
  y_centered <- y - median(y, na.rm = T)

  # Set seed for reproducible bootstrap results
  set.seed(123L)

  # Select the appropriate test - use twosamples for all tests for consistency
  assumption_test <- switch(test,
                            "ks" = twosamples::ks_test(x_centered, y_centered),
                            'kuiper' = twosamples::kuiper_test(x_centered, y_centered),
                            "cvm" = twosamples::cvm_test(x_centered, y_centered),
                            "ad" = twosamples::ad_test(x_centered, y_centered),
                            "wass" = twosamples::wass_test(x_centered, y_centered),
                            "dts" = twosamples::dts_test(x_centered, y_centered))

  #
  # Convert to list format compatible with the wmw_test() (extract p-value)
  test_pvalue <- as.numeric(assumption_test[2])  
  assumption_test <- list(
    statistic = as.numeric(assumption_test[1]),
    p.value = test_pvalue,
    method = paste(test, "test"),
    full.result = assumption_test  # Store full result for reference
  )
  #
  return(assumption_test)
  #
}








