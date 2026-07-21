#' P-value for the Bias-Corrected (BC) Test of AUC
#'
#' @description Tests \eqn{H_0\colon \mathrm{AUC} = A_0}{H0: AUC = A0} vs
#' the specified alternative, using a bias-corrected finite-sample variance
#' estimator with the mid-rank kernel.
#'
#' @param x Numeric vector of cases (group 1) values.
#' @param y Numeric vector of reference/control (group 2) values.
#' @param alternative Character: \code{"two.sided"}, \code{"greater"}, or
#'   \code{"less"}.
#' @param A0 Numeric null value of \eqn{\mathrm{AUC} = P(X < Y)}{AUC = P(X<Y)}.
#'   Defaults to 0.5.
#' @param min_n_warn_threshold Integer; if \code{min(length(x), length(y))}
#'   is below this threshold, a warning is issued that power may be very low
#'   at this sample size. Default 10.
#'
#' @return Numeric p-value.
#'
#' @details BC estimates \eqn{\mathrm{Var}(\hat A)}{Var(Ahat)} by correcting
#' each placement-variance component for its \eqn{O(1/n)} upward bias, using
#' a plug-in estimate of the bias subtracted from the naive placement
#' variance; each corrected component is floored independently at a small
#' \eqn{\epsilon > 0}{epsilon > 0} if it would otherwise go negative. The
#' mid-rank kernel \eqn{h(x,y) = 1\{x<y\} + \frac{1}{2} 1\{x=y\}}{h(x,y) =
#' 1(x<y) + 0.5*1(x=y)} is used throughout, for both the point estimate and
#' the variance components.
#'
#' Uses one-tier approach with \eqn{\hat\sigma^2_{\mathrm{adj}}}.
#'
#' BC is a conservative test: observed size stays below nominal across a wide
#' range of sample sizes, heteroskedasticity, and tie proportions, at a real
#' cost in power for small or imbalanced samples. See
#' \code{min_n_warn_threshold} and its warning text; the EU method
#' (\code{\link{wmwAUC_pvalue_EU}}) is recommended when \code{min(n1, n2)} is
#' small.
#'
#' \code{x} is taken to represent cases and \code{y} the reference/control
#' group, matching the convention of \code{wilcox.test()}. Internally, the
#' test statistic and variance components are computed in the
#' \eqn{P(X<Y)}{P(X<Y)} framework.
#'
#' @importFrom stats pt
#' @export
wmwAUC_pvalue_BC <- function(x, y, alternative = "two.sided", A0 = 0.5,
                             min_n_warn_threshold = 10) {
  n1 <- length(x); n2 <- length(y); n <- n1 + n2
  if (n1 < 3 || n2 < 3) stop("Sample sizes must be at least 3")
  if (!alternative %in% c("two.sided", "greater", "less")) {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  # See NEWS.md for the empirical basis of min_n_warn_threshold and the EU
  # recommendation below.
  if (min(n1, n2) < min_n_warn_threshold) {
    warning(sprintf(
      "min(n1,n2) = %d is small relative to the other group. Simulation shows ",
      min(n1, n2)),
      "BC's power can be close to zero at min(n1,n2) as large as 8, even in ",
      "balanced, homoskedastic, tie-free data (the easiest case for this test). ",
      "A non-significant result at this sample size is likely uninformative, ",
      "not evidence for the null. The EU method (wmwAUC_pvalue_EU()) ",
      "retained real power in simulations at this same ",
      "sample size and is recommended instead.")
  }
  
  out <- compute_t_stat_and_df_BC(x, y, A0 = A0)
  t_stat <- out$t_stat; df <- out$df
  
  if (alternative == "two.sided") {
    p_value <- 2 * pt(-abs(t_stat), df = df)
  } else if (alternative == "greater") {
    p_value <- pt(-t_stat, df = df)
  } else {
    p_value <- pt(t_stat, df = df)
  }
  p_value
}

compute_t_stat_and_df_BC <- function(x, y, A0 = 0.5, eps = 1e-8) {
  n1 <- length(x); n2 <- length(y); n <- n1 + n2
  lambda_n <- n1 / n
  
  # mid-rank placement values, consistent with the mid-rank point estimate below
  G_x_vals <- sapply(x, function(xi) mean((y < xi) + 0.5 * (y == xi)))
  F_y_vals <- sapply(y, function(yj) mean((x < yj) + 0.5 * (x == yj)))
  
  # mid-rank AUC point estimate, P(X < Y) convention: x = cases, y = reference
  auc_hat <- mean(outer(x, y, function(a, b) (a < b) + 0.5 * (a == b)))
  
  # Denominator is n (not n-1): E[(Ghat(x)-(1-A0))^2] centers at the fixed,
  # known value (1-A0), not an estimated sample mean, so Bessel's correction
  # is not justified here. 
  var_G_X <- sum((G_x_vals - (1 - A0))^2) / n1
  var_F_Y <- sum((F_y_vals - A0)^2) / n2
  
  # Bias correction
  omega_1_hat <- sum(G_x_vals * (1 - G_x_vals)) / n1
  omega_2_hat <- sum(F_y_vals * (1 - F_y_vals)) / n2
  
  var_G_X <- max(var_G_X - omega_1_hat / n2, eps)
  var_F_Y <- max(var_F_Y - omega_2_hat / n1, eps)
  
  w1 <- var_G_X / lambda_n
  w2 <- var_F_Y / (1 - lambda_n)
  sigma_sq_adj <- w1 + w2
  
  # 
  t_stat <- sqrt(n) * (auc_hat - A0) / sqrt(sigma_sq_adj)
  
  df_numerator <- sigma_sq_adj^2
  df_term1 <- w1^2 / (n1)
  df_term2 <- w2^2 / (n2)
  df <- df_numerator / (df_term1 + df_term2)
  
  list(t_stat = t_stat, df = df, auc_hat = auc_hat)
}

