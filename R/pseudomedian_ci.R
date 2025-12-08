#' Confidence Interval for Hodges-Lehmann Pseudomedian via Test Inversion
#' 
#' Computes confidence interval for the pseudomedian under \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5}
#' by test inversion.
#' 
#' @param x numeric vector, first sample
#' @param y numeric vector, second sample  
#' @param conf.level confidence level (default 0.95)
#' @param pvalue_method character, either 'EU' or 'BC'
#' @param n_grid number of grid points for search (default 1000)
#' @return list with conf.int, estimate and conf.level
#' 
#'
#' @importFrom stats mad
#' @importFrom stats pnorm
#' @importFrom stats var
#' 
#' @export
pseudomedian_ci <- function(x, y, conf.level = 0.95, 
                            pvalue_method = 'EU',
                            n_grid = 1000) {
  
  # Input validation
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors")
  }
  if (length(x) == 0 || length(y) == 0) {
    stop("x and y must have length > 0")
  }
  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be between 0 and 1")
  }
  
  alpha <- 1 - conf.level
  
  # Compute all pairwise differences to get search range
  all_diffs <- outer(x, y, "-")
  pseudomedian_est <- median(all_diffs)
  
  # Set search range around the estimate
  range_width <- 2 * mad(all_diffs, constant = 1.4826) # robust scale
  if (range_width == 0) range_width <- 1 # fallback for constant differences
  
  delta_min <- pseudomedian_est - 3 * range_width
  delta_max <- pseudomedian_est + 3 * range_width
  
  # Grid of candidate pseudomedian values
  delta_grid <- seq(delta_min, delta_max, length.out = n_grid)
  
  # Test H0: AUC = 0.5 for each shifted sample 
  # 
  n1 = length(x)
  n2 = length(y)
  n_kernel_evals <- n1 * n2
  total_computations <- n_kernel_evals * n_grid
  if (total_computations > 1e6) {
    warning(sprintf(
      "Large computational load detected: %dx%d = %s kernel evaluations x %d grid points = %s total computations. Consider reducing n_grid for faster execution.",
      n1, n2, 
      format(n_kernel_evals, big.mark = ","),
      n_grid,
      format(total_computations, big.mark = ",")
    ), call. = FALSE, immediate. = TRUE)
  }  

  if (pvalue_method == 'BC') { # using wmw_pvalue()
   #    
    p_values <- sapply(delta_grid, function(delta) {
      y_shifted <- y + delta
      wmw_pvalue(x, y_shifted, alternative = "two.sided")
    })
    #
  } else { # EU -> using wmw_pvalue_ties()
    #
    p_values <- sapply(delta_grid, function(delta) {
      y_shifted <- y + delta
      suppressWarnings(wmw_pvalue_ties(x, y_shifted, alternative = "two.sided"))
    })
    #
  }
  
  # Find confidence interval endpoints
  accepted <- delta_grid[p_values >= alpha]
  
  if (length(accepted) == 0) {
    warning("No accepted values found. Try increasing n_grid or check data.")
    ci <- c(NA, NA)
  } else {
    ci <- c(min(accepted), max(accepted))
  }
  
  list(
    conf.int = ci,
    estimate = pseudomedian_est,
    conf.level = conf.level
  )
}

