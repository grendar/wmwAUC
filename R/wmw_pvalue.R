#' P-value for Wilcoxon-Mann-Whitney Test of No Group Discrimination (Continuous Variables)
#' 
#' @description Tests \eqn{\mathrm{H_0\colon AUC} = 0.5}{AUC = 0.5} vs \eqn{\mathrm{H_1\colon AUC} \neq 0.5}{AUC != 0.5} 
#' with proper finite-sample corrections
#' 
#' @param x numeric vector for first group
#' @param y numeric vector for second group  
#' @param alternative character: "two.sided", "greater", or "less"
#' @return p-value
#' 
#' @details Implements the bias-corrected variance estimator with second-order
#' U-statistic correction to provide honest p-values under \eqn{\mathrm{H_0\colon AUC} = 0.5}{H₀: AUC = 0.5}.
#' Uses three-tier approach: permutation \eqn{(n < 20)}{(n < 20)}, 
#' bias-corrected \eqn{(20 \le n \lt 50)}{(20 <= n < 50)}, 
#' asymptotic with correction \eqn{n \ge 50}{(n >= 50)}.
#'   
#' For medium samples, the naive variance estimators \eqn{\widehat{\mathrm{Var}}(G(X))}{Var̂(G(X))} 
#' and \eqn{\widehat{\mathrm{Var}}(F(Y))}{Var̂(F(Y))} are 
#' corrected by subtracting O(1/n) bias terms of the form 
#' \eqn{(n_1 n_2)^{-1} \sum_i \hat{G}(X_i)(1 - \hat{G}(X_i))}
#' to prevent variance underestimation that would inflate Type I error rates.
#'
#' @export
wmw_pvalue <- function(x, y, alternative = "two.sided") {
  
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  lambda_n <- n1 / n
  
  # Validate inputs
  if (n1 < 3 || n2 < 3) {
    stop("Sample sizes must be at least 3")
  }
  if (!alternative %in% c("two.sided", "greater", "less")) {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  
  # Compute empirical AUC
  wt <- wilcox.test(x, y)
  W <- wt$statistic
  auc_hat <- as.numeric(W) / (n1 * n2)
  
  # 
  # Small samples: use permutation test
  if (n < 20) {
    return(wmw_permutation_test(x, y, alternative))
  }
  
  # Compute placement values
  # Ĝ(Xᵢ) = n₂⁻¹ Σⱼ 1{Yⱼ ≤ Xᵢ}
  G_x_vals <- sapply(x, function(xi) mean(y <= xi))
  
  # F̂(Yⱼ) = n₁⁻¹ Σᵢ 1{Xᵢ ≤ Yⱼ}  
  F_y_vals <- sapply(y, function(yj) mean(x <= yj))
  
  # Compute variance estimates
  # Var̂(G(X)) = (n₁−1)⁻¹ Σᵢ (Ĝ(Xᵢ) − 0.5)²
  var_G_X <- sum((G_x_vals - 0.5)^2) / (n1 - 1)
  
  # Var̂(F(Y)) = (n₂−1)⁻¹ Σⱼ (F̂(Yⱼ) − 0.5)²
  var_F_Y <- sum((F_y_vals - 0.5)^2) / (n2 - 1)
  
  # Apply bias correction for moderate samples
  if (n < 50) {
    # Subtract extra variance from estimating empirical CDFs
    # ω̂₁ = (1/n₁) Σᵢ Ĝ(Xᵢ)(1−Ĝ(Xᵢ))
    omega_1_hat <- sum(G_x_vals * (1 - G_x_vals)) / n1
    
    # ω̂₂ = (1/n₂) Σⱼ F̂(Yⱼ)(1−F̂(Yⱼ))  
    omega_2_hat <- sum(F_y_vals * (1 - F_y_vals)) / n2
    
    # Bias-corrected estimates (subtract the extra variance)
    var_G_X <- var_G_X - omega_1_hat / n2
    var_F_Y <- var_F_Y - omega_2_hat / n1
    
    # Check for over-correction
    if (var_G_X <= 0 || var_F_Y <= 0) {
      warning("Bias correction too aggressive - using uncorrected estimates")
      var_G_X <- sum((G_x_vals - 0.5)^2) / (n1 - 1)
      var_F_Y <- sum((F_y_vals - 0.5)^2) / (n2 - 1)
    }
  }
  
  # Welch-type variance combination  
  sigma_sq_adj <- var_G_X / lambda_n + var_F_Y / (1 - lambda_n)
  
  # Apply second-order U-statistic correction
  # σ̂²_{final} = (1 − 1/n₁ − 1/n₂) · σ̂²_{adj}
  correction_factor <- (1 - 1/n1 - 1/n2)
  sigma_sq_final <- correction_factor * sigma_sq_adj
  
  # Test statistic
  t_stat <- sqrt(n) * (auc_hat - 0.5) / sqrt(sigma_sq_final)
  
  # Satterthwaite degrees of freedom
  # df = (σ̂²)² / [(Var̂(G(X))/λₙ)²/(n₁−1) + (Var̂(F(Y))/(1−λₙ))²/(n₂−1)]
  df_numerator <- sigma_sq_final^2
  df_term1 <- (var_G_X / lambda_n * correction_factor)^2 / (n1 - 1)
  df_term2 <- (var_F_Y / (1 - lambda_n) * correction_factor)^2 / (n2 - 1)
  df <- df_numerator / (df_term1 + df_term2)
  
  # Compute p-value
  if (alternative == "two.sided") {
    p_value <- 2 * pt(-abs(t_stat), df = df)
  } else if (alternative == "greater") {
    p_value <- pt(-t_stat, df = df)  # H₁: AUC > 0.5
  } else {
    p_value <- pt(t_stat, df = df)   # H₁: AUC < 0.5  
  }
  
  # Determine method used
  if (n < 20) {
    method <- "Honest WMW test (permutation)"
  } else if (n < 50) {
    method <- "Honest WMW test (bias-corrected)"
  } else {
    method <- "Honest WMW test (asymptotic)"
  }
  
  return(p_value)
  #  
  
}



# Placeholder for permutation test (for n < 20)
wmw_permutation_test <- function(x, y, alternative) {
  # Simple permutation test implementation
  observed_auc <- mean(outer(x, y, ">")) + 0.5 * mean(outer(x, y, "=="))
  
  # Permutation distribution under H₀: AUC = 0.5
  pooled <- c(x, y)
  n1 <- length(x)
  n2 <- length(y)
  
  perm_aucs <- replicate(2000, {
    perm_indices <- sample(length(pooled))
    perm_x <- pooled[perm_indices[1:n1]]
    perm_y <- pooled[perm_indices[(n1+1):(n1+n2)]]
    mean(outer(perm_x, perm_y, ">")) + 0.5 * mean(outer(perm_x, perm_y, "=="))
  })
  
  # P-value calculation
  if (alternative == "two.sided") {
    p_value <- mean(abs(perm_aucs - 0.5) >= abs(observed_auc - 0.5))
  } else if (alternative == "greater") {
    p_value <- mean(perm_aucs >= observed_auc)
  } else {
    p_value <- mean(perm_aucs <= observed_auc)
  }
  
  return(list(
    statistic = c("AUC" = observed_auc),
    parameter = NULL,
    p.value = p_value,
    auc = observed_auc,
    method = "Honest WMW test (permutation)",
    alternative = alternative,
    sigma.sq = NA,
    correction.factor = NA
  ))
}

