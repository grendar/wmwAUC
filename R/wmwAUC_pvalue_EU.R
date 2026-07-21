#' P-value for the Exact Unbiased (EU) Test of AUC
#'
#' @description Tests \eqn{H_0\colon \mathrm{AUC} = A_0}{H0: AUC = A0} vs
#' the specified alternative, using an exact finite-sample unbiased variance
#' estimator with the mid-rank kernel. 
#'
#' @param x Numeric vector of cases (group 1) values.
#' @param y Numeric vector of reference/control (group 2) values.
#' @param alternative Character: \code{"two.sided"}, \code{"greater"}, or
#'   \code{"less"}.
#' @param A0 Numeric null value of \eqn{\mathrm{AUC} = P(X < Y)}{AUC = P(X<Y)}.
#'   Defaults to 0.5. Only supported when \code{length(x) + length(y) >= 20}
#'   (see Details); an error is raised otherwise.
#' @param max_exact Integer; the permutation branch (used when
#'   \code{length(x) + length(y) < 20}) enumerates all permutations exactly
#'   when their count is at most \code{max_exact}, and falls back to Monte
#'   Carlo sampling above that. Default 10000.
#' @param n_perm Integer; number of Monte Carlo permutation replicates used
#'   when exact enumeration is not feasible. Default 2000.
#'
#' @return Numeric p-value.
#'
#' @details Uses two-tier approach: studentized permutation 
#' for \code{length(x) + length(y) < 20} and the exact finite-sample unbiased
#' estimator for \code{length(x) + length(y) >= 20}.
#' 
#' For \code{length(x) + length(y) < 20}, a studentized permutation test is
#' used: the same t-statistic is recomputed on each permuted split of
#' the pooled data, and the p-value is the proportion of permuted statistics
#' at least as extreme as the observed one. This permutation scheme relies on
#' group-relabeling exchangeability, which preserves \eqn{H_0}{H0} only when
#' \eqn{A_0 = 0.5}{A0 = 0.5}; general \code{A0} is therefore not currently
#' supported in this small-sample regime and will raise an error.
#'
#' For \code{length(x) + length(y) >= 20}, EU estimates
#' \eqn{\mathrm{Var}(\hat A)}{Var(Ahat)} by the exact finite-sample unbiased
#' combination derived from the Hoeffding decomposition of the mid-rank
#' kernel, with Welch--Satterthwaite degrees of freedom.
#'
#'
#' \code{x} is taken to represent cases and \code{y} the reference/control
#' group, matching the convention of \code{wilcox.test()}. Internally, the
#' test statistic and variance components are computed in the
#' \eqn{P(X<Y)}{P(X<Y)} framework.
#'
#' @export
wmwAUC_pvalue_EU <- function(x, y, alternative = "two.sided", A0 = 0.5, max_exact = 10000, n_perm = 2000) {
  n1 <- length(x); n2 <- length(y); n <- n1 + n2
  if (n1 < 3 || n2 < 3) stop("Sample sizes must be at least 3")
  if (!alternative %in% c("two.sided", "greater", "less")) {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  if (n < 20) {
    if (A0 != 0.5) {
      stop("General A0 is not supported for n1+n2<20: the permutation branch ",
           "relies on group-relabeling exchangeability, which only preserves ",
           "H0 when A0=0.5 (relabeling maps A0 to 1-A0 otherwise). This is a ",
           "known limitation deferred to a future version, not a validated ",
           "restriction workaround -- do not bypass it.")
    }
    return(studentized_permutation_test_EU(x, y, alternative, A0 = A0, max_exact = max_exact, n_perm = n_perm))
  }
  out <- compute_stat_ties_EU(x, y, A0 = A0)
  t_stat <- out$t_stat; df <- out$df
  if (alternative == "two.sided") {
    p_value <- 2 * stats::pt(-abs(t_stat), df = df)
  } else if (alternative == "greater") {
    p_value <- stats::pt(-t_stat, df = df)
  } else {
    p_value <- stats::pt(t_stat, df = df)
  }
  p_value
}


compute_stat_ties_EU <- function(x, y, A0 = 0.5) {
  n1 <- length(x); n2 <- length(y); n <- n1 + n2; M <- n1 * n2
  
  h <- outer(x, y, function(a, b) (a < b) + 0.5 * (a == b))  # mid-rank, P(X<Y) kernel directly
  Ahat <- mean(h)
  
  v_hat <- sum((h - Ahat)^2) / (M - 1)
  zeta1_hat2 <- if (n1 > 1) var(rowMeans(h)) else 0
  zeta2_hat2 <- if (n2 > 1) var(colMeans(h)) else 0
  
  a <- -(M - 1) / (M * (n1 - 1) * (n2 - 1))
  b <- n2^2 / (M * (n2 - 1))
  c <- n1^2 / (M * (n1 - 1))
  Var_hat_unbiased <- a * v_hat + b * zeta1_hat2 + c * zeta2_hat2
  
  if (is.na(Var_hat_unbiased) || Var_hat_unbiased <= 0) {
    Var_hat_unbiased <- (zeta1_hat2 + zeta2_hat2) / 2 / n  # unvalidated fallback, unchanged
  }
  
  t_stat <- (Ahat - A0) / sqrt(Var_hat_unbiased)
  
  if (M >= 4 && n1 >= 3 && n2 >= 3) {
    comp1 <- b * zeta1_hat2; comp2 <- c * zeta2_hat2; comp3 <- a * v_hat
    df_num <- Var_hat_unbiased^2
    df_den <- comp1^2 / (n1 - 1) + comp2^2 / (n2 - 1) + comp3^2 / (M - 1)
    df <- if (df_den > 0) max(1, min(df_num / df_den, n - 2)) else n - 2
  } else {
    df <- max(1, n - 2)
  }
  
  list(t_stat = t_stat, df = df, Ahat = Ahat)
}

# Studentized permutation test for small samples. Reuses
# compute_stat_ties_EU()'s t-statistic on each permuted split, so the
# reference distribution is that of the studentized statistic rather than
# the raw AUC deviation -- see NEWS.md for why this matters and its
# validated operating characteristics.
#
# This permutation scheme (relabeling X/Y under group exchange) is only
# exact for A0=0.5 -- relabeling maps A0 to 1-A0, so it does not preserve H0
# for other A0. Do not call this with A0 != 0.5; the top-level
# wmwAUC_pvalue_EU() guards against that below.
studentized_permutation_test_EU <- function(x, y, alternative, A0 = 0.5, max_exact = 10000, n_perm = 2000) {
  n1 <- length(x); n2 <- length(y); n_total <- n1 + n2
  T_obs <- compute_stat_ties_EU(x, y, A0 = A0)$t_stat
  pooled <- c(x, y)
  
  perm_stat <- function(indices) {
    perm_x <- pooled[indices]; perm_y <- pooled[-indices]
    out <- compute_stat_ties_EU(perm_x, perm_y, A0 = A0)
    out$t_stat
  }
  
  if (choose(n_total, n1) <= max_exact) {
    all_combs <- utils::combn(n_total, n1)
    T_perm <- apply(all_combs, 2, perm_stat)
  } else {
    T_perm <- replicate(n_perm, perm_stat(sample(n_total, n1)))
  }
  T_perm <- T_perm[is.finite(T_perm)]  # guard against degenerate permuted variance
  
  if (alternative == "two.sided") {
    mean(abs(T_perm) >= abs(T_obs))
  } else if (alternative == "greater") {
    mean(T_perm >= T_obs)
  } else {
    mean(T_perm <= T_obs)
  }
}

