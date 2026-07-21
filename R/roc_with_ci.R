#' ROC related computations -- internal function
#'
#' @param probs Vector of class probabilities or values of continuous predictor
#' @param labels Vector, factor with two levels
#' @param positive Character giving the level that corresponds to 'case'
#' @param auc Numeric value of AUC
#' @param ci_method Character from c("none", "delong", "bootstrap")
#' @param n_boot Numeric value giving the number of bootstrap replicates (default: 1000)
#' @param alpha  Level of significance (default: 0.05)
#'
#' @details \code{ci_method = 'bootstrap'} makes this function's result
#' stochastic (resampled AUC and ROC confidence band). For reproducible
#' results, call \code{set.seed()} immediately beforehand.
#'
#' @return List with components:
#'   \item{roc_df}{data frame for plotting ROC curve}
#'   \item{roc_band}{data frame for plotting confidence band of ROC}
#'   \item{auc}{auc}
#'   \item{auc_ci}{confint for auc}
#'
#' @importFrom stats approx median qnorm quantile relevel wilcox.test
#' @importFrom utils head tail
#'
#' @keywords internal
#' @export
roc_with_ci <- function(
    probs,
    labels,
    positive,
    auc,
    ci_method = c("none", "delong", "bootstrap"),
    n_boot = 1000,
    alpha = 0.05
) {
  ci_method <- match.arg(ci_method)
  
  # Convert labels to numeric 0/1 (1 = positive)
  y <- ifelse(labels == positive, 1, 0)
  
  # Sort by decreasing score
  ord <- order(probs, decreasing = TRUE)
  probs_ord <- probs[ord]
  y_ord <- y[ord]
  
  ##############################################################################
  # ROC computation
  ##############################################################################
  #
  compute_roc <- function(y_sorted) {
    P <- sum(y_sorted == 1)
    N <- sum(y_sorted == 0)
    
    TPR <- cumsum(y_sorted == 1) / P
    FPR <- cumsum(y_sorted == 0) / N
    
    data.frame(fpr = c(0, FPR), tpr = c(0, TPR))
  }
  
  roc_df <- compute_roc(y_ord)
  
  # ci for AUC
  auc_ci <- c(NA, NA)
  roc_band <- NULL
  
  ##############################################################################
  # DeLong et al. (1988) structural-components CI
  ##############################################################################
  #
  if (ci_method == "delong") {
    x_case <- probs[y == 1]
    y_ref  <- probs[y == 0]
    m <- length(x_case)
    n <- length(y_ref)
    
    # mid-rank kernel, consistent with the auc = P(case > reference)
    # convention used throughout this package
    psi_mat <- outer(x_case, y_ref, function(a, b) (a > b) + 0.5 * (a == b))
    
    V10 <- rowMeans(psi_mat)  # one structural component per case
    V01 <- colMeans(psi_mat)  # one structural component per reference obs
    
    S10 <- var(V10)  # var() default: divide by (m-1)
    S01 <- var(V01)  # divide by (n-1)
    
    se_auc <- sqrt(S10 / m + S01 / n)
    
    z <- qnorm(1 - alpha / 2)
    auc_ci <- c(auc - z * se_auc, auc + z * se_auc)
  }
  
  ##############################################################################
  # Bootstrap CI for AUC + ROC band
  ##############################################################################
  #
  if (ci_method == "bootstrap") {
    auc_boot <- numeric(n_boot)
    
    # Common FPR grid for ROC confidence band
    fpr_grid <- seq(0, 1, length.out = 200)
    tpr_mat <- matrix(NA, nrow = n_boot, ncol = length(fpr_grid))
    
    for (b in seq_len(n_boot)) {
      idx <- sample(seq_along(y), replace = TRUE)
      
      probs_b <- probs[idx]
      y_b <- y[idx]
      
      ord_b <- order(probs_b, decreasing = TRUE)
      roc_b <- compute_roc(y_b[ord_b])
      
      # Interpolate ROC on common grid
      tpr_mat[b, ] <- suppressWarnings(approx(roc_b$fpr, roc_b$tpr,
                                              xout = fpr_grid, rule = 2)$y)
      
      # Bootstrap AUC
      auc_boot[b] <- sum(diff(roc_b$fpr) *
                           (head(roc_b$tpr, -1) + tail(roc_b$tpr, -1)) / 2)
    }
    
    # AUC CI from bootstrap percentiles
    auc_ci <- quantile(auc_boot, c(alpha/2, 1 - alpha/2))
    
    # Pointwise ROC confidence band
    roc_band <- data.frame(
      fpr = fpr_grid,
      tpr_lo = apply(tpr_mat, 2, quantile, alpha / 2),
      tpr_hi = apply(tpr_mat, 2, quantile, 1 - alpha / 2)
    )
  }
  #
  ##############################################################################
  #
  ##############################################################################
  #
  roc_ci_list <- list(roc_df = roc_df, roc_band = roc_band, auc = auc, auc_ci = auc_ci)
  #
  class(roc_ci_list) <- "roc_ci"
  return(roc_ci_list) 
  
}