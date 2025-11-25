#' ROC related computations -- internal function
#'
#' @param probs Vector of class probabilities or values of continuous predictor
#' @param labels Vector, factor with two levels
#' @param positive Character giving the level that corresponds to 'case'
#' @param auc Numeric value of AUC
#' @param ci_method Character from c("none", "hanley", "bootstrap")
#' @param n_boot Numeric value giving the number of bootstrap replicates (default: 1000)
#' @param alpha  Level of significance (default: 0.05)
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
    ci_method = c("none", "hanley", "bootstrap"),
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
  # Hanley & McNeil (analytical) CI
  ##############################################################################
  #
  if (ci_method == "hanley") {
    P <- sum(y == 1)
    N <- sum(y == 0)

    Q1 <- auc / (2 - auc)
    Q2 <- 2 * auc^2 / (1 + auc)

    se_auc <- sqrt((auc * (1 - auc) +
                      (P - 1) * (Q1 - auc^2) +
                      (N - 1) * (Q2 - auc^2)) / (P * N))

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
      tpr_mat[b, ] <- approx(roc_b$fpr, roc_b$tpr,
                             xout = fpr_grid, rule = 2)$y

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
