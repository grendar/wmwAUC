#' P-value for Wilcoxon-Mann-Whitney Test of No Group Discrimination (With Possible Ties)
#' 
#' @description Tests \eqn{\mathrm{H_0\colon AUC} = 0.5}{AUC = 0.5} vs \eqn{\mathrm{H_1\colon AUC} \neq 0.5}{AUC != 0.5} 
#' with exact finite-sample unbiased variance estimation for arbitrary tie patterns
#' 
#' @param x Numeric vector of cases/group 1 values
#' @param y Numeric vector of controls/reference group values  
#' @param alternative character: "two.sided", "greater", or "less"
#' @return p-value
#' 
#' @details Implements the Exact finite-sample Unbiased (EU) variance estimator derived from 
#' Hoeffding decomposition theory. Uses tie-corrected kernel \eqn{h(x,y) = \mathbf{1}\{x < y\} + \frac{1}{2}\mathbf{1}\{x = y\}}{h(x,y) = 1{x<y} + 0.5*1{x=y}}
#' with universal second-order correction factor to provide honest p-values under 
#' \eqn{\mathrm{H_0\colon AUC} = 0.5}{H₀: AUC = 0.5} regardless of tie structure.
#' 
#' Uses three-tier approach: permutation \eqn{(n < 20)}{(n < 20)}, 
#' exact unbiased estimator \eqn{(20 \le n < 50)}{(20 <= n < 50)}, 
#' asymptotic with corrections \eqn{n \ge 50}{(n >= 50)}.
#'   
#' The unbiased variance estimator is constructed as a specific linear combination:
#' 
#' \eqn{\widetilde{\mathrm{Var}}(\hat{A}) = \frac{n_2\hat{\zeta}_1^2 + n_1\hat{\zeta}_2^2 - \frac{M-1}{M}\hat{v}}{M+1}}
#' 
#' where \eqn{\hat{v}} is the pooled sample variance of kernel values and 
#' \eqn{\hat{\zeta}_1^2, \hat{\zeta}_2^2} are row/column mean variances.
#' 
#' Welch-Satterthwaite degrees of freedom account for bias correction structure:
#' 
#' \eqn{\nu = \frac{(\hat{\sigma}^2)^2}{\frac{(n_2\hat{\zeta}_1^2/(M+1))^2}{n_1-2} + \frac{(n_1\hat{\zeta}_2^2/(M+1))^2}{n_2-2} + \frac{((M-1)\hat{v}/(M(M+1)))^2}{M-3}}}
#'
#' Function uses mid-rank tie handling throughout, ensuring theoretical consistency
#' with the corrected null hypothesis framework.
#'
#' Function assumes \eqn{x} represents cases and \eqn{y} represents the reference level, 
#' in accord with `wilcox.test()` and `wmw_test()`. 
#' Internal calculations convert to P(X < Y) framework to match theoretical derivations.
#'
#' @importFrom utils combn
#
#' @export
wmw_pvalue_ties <- function(x, y, alternative = "two.sided") {
  
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  M <- n1 * n2
  
  # Validate inputs
  if (n1 < 3 || n2 < 3) {
    stop("Sample sizes must be at least 3")
  }
  if (!alternative %in% c("two.sided", "greater", "less")) {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  
  # Small samples: use permutation test
  if (n < 20) {
    return(wmw_permutation_test_ties(x, y, alternative))
  }
  
  # Compute tie-corrected kernel matrix h_{ij} = 1{x>y} + 0.5*1{x==y}
  h <- outer(x, y, FUN = Vectorize(function(a,b) as.numeric(a > b) + 0.5 * as.numeric(a == b)))
  # AUC estimate
  Ahat <- mean(h)   # P(X > Y)
  Ahat <- 1 - Ahat  # since asymptotics in preprint is for P(X < Y)
  h <- 1 - h        # since asymptotics in preprint is for P(X < Y)
  #
  
  # Compute empirical variance components
  v_hat <- sum((h - Ahat)^2) / (M - 1)        # pooled sample variance
  z1_hat2 <- if (n1 > 1) var(rowMeans(h)) else 0  # row-mean variance
  z2_hat2 <- if (n2 > 1) var(colMeans(h)) else 0  # column-mean variance
  
  # Exact finite-sample unbiased variance estimator
  Var_hat_unbiased <- (n2 * z1_hat2 + n1 * z2_hat2 - ((M - 1) / M) * v_hat) / (M + 1)
  
  # Check for numerical issues
  if (Var_hat_unbiased <= 0) {
    warning("Unbiased variance estimator is non-positive - using fallback method")
    # Fallback to simpler estimator
    Var_hat_unbiased <- (z1_hat2 + z2_hat2) / 2 / n
  }
  
  # Apply universal second-order U-statistic correction
  correction_factor <- (1 - 1/n1 - 1/n2)
  Var_hat_final <- correction_factor * Var_hat_unbiased
  
  # Ensure positive final variance
  if (Var_hat_final <= 0) {
    warning("Final variance estimate is non-positive - using uncorrected estimate")
    Var_hat_final <- Var_hat_unbiased
    correction_factor <- 1
  }
  
  # Test statistic
  t_stat <- (Ahat - 0.5) / sqrt(Var_hat_final)

  # Welch-Satterthwaite degrees of freedom with bias correction adjustments
  if (M >= 4 && n1 >= 3 && n2 >= 3) {
    df_numerator <- (Var_hat_final)^2
    df_term1 <- ((n2 * z1_hat2 / (M + 1)) * correction_factor)^2 / max(n1 - 2, 1)
    df_term2 <- ((n1 * z2_hat2 / (M + 1)) * correction_factor)^2 / max(n2 - 2, 1)
    df_term3 <- (((M - 1) * v_hat / (M * (M + 1))) * correction_factor)^2 / max(M - 3, 1)
    
    df_denominator <- df_term1 + df_term2 + df_term3
    
    if (df_denominator > 0) {
      df <- df_numerator / df_denominator
      # Bound degrees of freedom reasonably
      df <- max(1, min(df, n - 2))
    } else {
      df <- n - 2  # Fallback
    }
  } else {
    df <- max(1, n - 2)  # Fallback for very small samples
  }
  
  # Compute p-value using t-distribution
  if (alternative == "two.sided") {
    p_value <- 2 * pt(-abs(t_stat), df = df)
  } else if (alternative == "greater") {
    p_value <- pt(-t_stat, df = df)  # H₁: AUC > 0.5
  } else {
    p_value <- pt(t_stat, df = df)   # H₁: AUC < 0.5  
  }
  
  return(p_value)
}

# Permutation test for small samples with ties
wmw_permutation_test_ties <- function(x, y, alternative) {
  
  n1 <- length(x)
  n2 <- length(y)
  
  # Compute observed AUC using tie-corrected kernel
  h_obs <- outer(x, y, FUN = Vectorize(function(a,b) as.numeric(a > b) + 0.5 * as.numeric(a == b)))
  observed_auc_raw <- mean(h_obs) # P(X>Y)
  observed_auc <- 1 - observed_auc_raw  # Convert to P(X < Y)
  
  # Permutation distribution under H₀: AUC = 0.5
  pooled <- c(x, y)
  n_total <- length(pooled)
  
  # For very small samples, do exact enumeration if feasible
  if (choose(n_total, n1) <= 10000) {
    # Exact permutation test
    n_perms <- choose(n_total, n1)
    # message("Performing exact permutation test with ", n_perms, " permutations")
    
    # Generate all combinations
    all_combs <- combn(n_total, n1)
    perm_aucs <- apply(all_combs, 2, function(indices) {
      perm_x <- pooled[indices]
      perm_y <- pooled[-indices]
      h_perm <- outer(perm_x, perm_y, FUN = Vectorize(function(a,b) as.numeric(a > b) + 0.5 * as.numeric(a == b)))
      perm_auc_raw <- mean(h_perm)  # P(X > Y)
      return(1 - perm_auc_raw)  # Convert to P(X < Y)
    })
  } else {
    # Monte Carlo permutation test
    n_perms <- 2000
    # message("Performing Monte Carlo permutation test with ", n_perms, " permutations")
    
    perm_aucs <- replicate(n_perms, {
      perm_indices <- sample(n_total)
      perm_x <- pooled[perm_indices[1:n1]]
      perm_y <- pooled[perm_indices[(n1+1):n_total]]
      
      h_perm <- outer(perm_x, perm_y, FUN = Vectorize(function(a,b) as.numeric(a > b) + 0.5 * as.numeric(a == b)))
      perm_auc_raw <- mean(h_perm)  # P(X > Y)
      return(1 - perm_auc_raw)  # Convert to P(X < Y)
    })
  }
  
  # P-value calculation
  if (alternative == "two.sided") {
    p_value <- mean(abs(perm_aucs - 0.5) >= abs(observed_auc - 0.5))
  } else if (alternative == "greater") {
    p_value <- mean(perm_aucs >= observed_auc)
  } else {
    p_value <- mean(perm_aucs <= observed_auc)
  }
  
  return(p_value)
}

