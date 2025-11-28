#' Confidence Interval for Hodges-Lehmann Pseudomedian via Test Inversion
#' 
#' Computes confidence interval for the pseudomedian under \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5}
#' by test inversion.
#' 
#' @param x numeric vector, first sample
#' @param y numeric vector, second sample  
#' @param conf.level confidence level (default 0.95)
#' @param n_grid number of grid points for search (default 1000)
#' @return list with conf.int, estimate and conf.level
#' 
#'
#' @importFrom stats mad
#' @importFrom stats pnorm
#' @importFrom stats var
#' 
#' @export
pseudomedian_ci <- function(x, y, conf.level = 0.95, n_grid = 1000) {
  
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
  
  # Test H0: AUC = 0.5 for each shifted sample using wmw_pvalue()
  p_values <- sapply(delta_grid, function(delta) {
    y_shifted <- y + delta
    wmw_pvalue(x, y_shifted, alternative = "two.sided")
  })
  
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

#' Test H0: AUC = 0.5 (wilcox.test-style edge case handling)
#' 
#' Test for AUC = 0.5 following wilcox.test edge case philosophy
#' 
#' @param x first sample
#' @param y second sample
#' @return p-value for testing H0: AUC = 0.5
test_auc_half <- function(x, y) {
  
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  
  # Input validation (like wilcox.test)
  if (n1 < 1 || n2 < 1) {
    stop("not enough observations")
  }
  
  # Compute empirical AUC using wilcox.test() statistic (consistent with wmwAUC)
  wt_stat <- wilcox.test(x, y, exact = FALSE)$statistic
  auc_emp <- as.numeric(wt_stat) / (n1 * n2)
  
  # Placement values: G(X_i) and F(Y_j)
  G_at_X <- sapply(x, function(xi) mean(y <= xi))
  F_at_Y <- sapply(y, function(yj) mean(x <= yj))
  
  # Simple variance estimates (no bias correction)
  var_G_X <- var(G_at_X)
  var_F_Y <- var(F_at_Y) 
  
  # Sample proportion
  lambda_n <- n1 / n
  
  # Asymptotic variance under H0: AUC = 0.5
  sigma2 <- var_G_X / lambda_n + var_F_Y / (1 - lambda_n)
  
  # Edge case handling inspired by wilcox.test:
  # Use mathematical approach, not conservative p=1.0
  if (sigma2 <= 1e-12) {
    # When empirical variance is essentially zero,
    # use the general asymptotic variance formula with empirical estimates
    
    # Under general H0: AUC = 0.5, asymptotic variance is:
    # σ² = Var(G(X))/λ + Var(F(Y))/(1-λ)
    # When both are degenerate, fall back to simple test based on AUC distance
    
    z_stat <- abs(auc_emp - 0.5) * sqrt(12 * n1 * n2 / (n1 + n2))
    return(2 * (1 - pnorm(abs(z_stat))))
  }
  
  # Normal case: use empirical variance estimate
  z_stat <- sqrt(n) * (auc_emp - 0.5) / sqrt(sigma2)
  
  # Two-sided p-value
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  
  # Never return exactly 1.0 unless mathematically justified
  # (wilcox.test only returns 1.0 for trivial cases like n=1 vs n=1)
  if (p_value > 0.999999 && n1 > 1 && n2 > 1) {
    p_value <- 0.999999
  }
  
  return(p_value)
}

