#' wmwAUC Test of No Group Discrimination
#'
#' Performs the wmwAUC test of \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5} based on
#' the Wilcoxon-Mann-Whitney statistic.
#' 
#'
#' @param formula Formula of the form `response ~ group`
#' @param data Data frame containing continuous response variable and grouping factor
#' @param ref_level Character, reference level of grouping factor (if NULL, uses first level)
#' @param pvalue_method Character, method ('EU', 'BC') used for computing p-values (default 'EU')
#' @param alternative Character, alternative hypothesis is c("two.sided", "greater", "less")
#' @param ci_method Character, confidence interval method for eAUC: c('delong', 'boot', 'none')
#' @param conf_level Numeric, confidence level for intervals (default 0.95)
#' @param pseudomedian Logical; if TRUE, additionally reports the AUC-equalizing
#'   shift (pseudomedian) and its confidence interval (default FALSE)
#' @param n_grid Numeric, number of grid points for search in `wmwAUC_pseudomedian_ci()` (default 1000)
#' @param ... Additional arguments passed to `roc_with_ci()`
#'
#' @details
#' The function tests the null hypothesis \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5}
#' against \eqn{\mathrm{H_1\colon AUC} \neq 0.5}{H1: AUC != 0.5},
#' where AUC represents the Area Under the ROC Curve.
#'
#' The Exact Unbiased ('EU') method is used by default for computing p-values.
#' Bias-Corrected ('BC') method is available through `pvalue_method = 'BC'`
#' and is markedly conservative at small or imbalanced sample sizes; EU is
#' recommended unless BC's specific properties are wanted (see
#' \code{\link{wmwAUC_pvalue_BC}}).
#'
#' Following the convention of `wilcox.test()` AUC equals the probability \eqn{P(X > Y)} 
#' that a randomly selected observation from the first group exceeds a randomly
#' selected observation from the second group.
#' For `response ~ group`, observations from the non-reference group constitute \eqn{X},
#' while observations from the reference group (specified by `ref_level`) constitute \eqn{Y}.
#' Thus AUC = P(non-reference > reference). If `ref_level` is not specified, the first
#' factor level is used as reference. The \eqn{U}-statistic and the resulting empirical AUC (eAUC)
#' are calculated consistently with this group assignment.
#'
#' The test statistic is eAUC, which estimates the true AUC.
#' The empirical ROC curve (eROC) is constructed by varying the classification
#' threshold across all observed values and computing sensitivity and 1-specificity
#' at each threshold.
#'
#' When `pseudomedian = TRUE`, the function additionally reports the AUC-equalizing
#' shift \eqn{\delta}, defined as the value solving
#' \eqn{P(X < Y + \delta) = 0.5}{P(X < Y + delta) = 0.5}; see
#' \code{\link{wmwAUC_pseudomedian_ci}} for details. 
#'
#' Confidence intervals for the true AUC are computed using either the
#' DeLong et al. (1988) structural-components method based on asymptotic
#' normality, or bootstrap resampling. 
#' If bootstrap resampling is selected, it is also used for constructing the
#' confidence band for the ROC curve.
#'
#' This function can call two independent sources of randomness: bootstrap
#' resampling (\code{ci_method = 'boot'}), and, when \code{pseudomedian =
#' TRUE} with a small sample (\code{n1 + n2 < 20}), \code{\link{wmwAUC_pvalue_EU}}'s
#' Monte Carlo permutation fallback, called once per grid point inside
#' \code{\link{wmwAUC_pseudomedian_ci}}. For reproducible results, call
#' \code{set.seed()} immediately before \code{wmwAUC_test()} rather than
#' relying on the ambient RNG state; a single seed covers both sources, since
#' they draw from the same stream in a fixed order.
#'
#' @return Object of class 'wmwAUC_test' containing:
#'   \item{pseudomedian_requested}{Logical indicating whether the pseudomedian was computed}
#'   \item{n}{Named vector with components n1, n2 giving sample sizes for each group}
#'   \item{U_statistic}{U statistic}
#'   \item{p_value}{P-value for testing H0: AUC = 0.5}
#'   \item{alternative}{Alternative hypothesis specification}
#'   \item{pvalue_method}{Character string describing the test method}
#'   \item{data_name}{Character string giving the name of the data}
#'   \item{pseudomedian}{AUC-equalizing shift estimate (when pseudomedian = TRUE)}
#'   \item{pseudomedian_conf_int}{Confidence interval for AUC-equalizing shift (when pseudomedian = TRUE)}
#'   \item{pseudomedian_conf_level}{Confidence level for the pseudomedian interval (when pseudomedian = TRUE)}
#'   \item{ci_method}{Method used to compute confidence interval for AUC}
#'   \item{roc_object}{ROC analysis object returned by `roc_with_ci()`}
#'   \item{auc}{Empirical AUC (eAUC), the standardized U statistic}
#'   \item{auc_conf_int}{Confidence interval for true AUC using DeLong et al. or bootstrap method}
#'   \item{x_vals}{Numeric vector of observations from non-reference group}
#'   \item{y_vals}{Numeric vector of observations from reference group}
#'   \item{groups}{Character vector of group labels from original data}
#'   \item{group_levels}{Character vector of factor levels for grouping variable}
#'   \item{group_ref_level}{Character string indicating which level corresponds to reference group}
#'
#' @examples
#' library('wmwAUC')
#' \donttest{
#' library('gemR')
#' data(MS)
#' da <- MS
#' # preparing data frame
#' class(da$proteins) <- setdiff(class(da$proteins), "AsIs")
#' df <- as.data.frame(da$proteins)
#' df$MS <- da$MS
#' # wmwAUC test
#' wmd <- wmwAUC_test(P19099 ~ MS, data = df, ref_level = 'no')
#' wmd
#' plot(wmd)
#' # compute pseudomedian
#' set.seed(123L)
#' wmd_pm <- wmwAUC_test(P19099 ~ MS, data = df, ref_level = 'no', pseudomedian = TRUE)
#' wmd_pm
#' # compute confint for AUC by bootstrap
#' set.seed(123L)
#' wmd_ci_boot <- wmwAUC_test(P19099 ~ MS, data = df, ref_level = 'no', ci_method = 'boot')
#' wmd_ci_boot
#' plot(wmd_ci_boot)
#' # BC method
#' wmd_bc <- wmwAUC_test(P19099 ~ MS, data = df, ref_level = 'no', pvalue_method = 'BC')
#' wmd_bc
#' }
#'
#' @references
#'
#' Grendar, M. (2025). Wilcoxon-Mann-Whitney test of no group discrimination.
#' arXiv:2511.20308. (Full bibliography, including all methods and sources
#' cited throughout this package, is given there.)
#'
#' @seealso
#' \code{\link{print.wmwAUC_test}} for formatted output of `wmwAUC_test()`.
#' \code{\link{plot.wmwAUC_test}} for plot of output of `wmwAUC_test()`.
#' \code{\link{wmwAUC_pvalue_BC}} for details on computing p-values using the 'BC' method.
#' \code{\link{wmwAUC_pvalue_EU}} for details on computing p-values using the 'EU' method.
#' \code{\link{wmwAUC_pseudomedian_ci}} for details on computing confidence intervals for the pseudomedian.
#' \code{\link{wilcox.test}} for the Wilcoxon-Mann-Whitney test in base R.
#' \code{\link[rankFD]{rank.two.samples}} in package \pkg{rankFD} for an
#' implementation of the Brunner-Munzel test.
#'
#' @export
wmwAUC_test <- function(formula,
                        data,
                        ref_level = NULL,
                        pvalue_method = c("EU", "BC"),
                        alternative = c("two.sided", "greater", "less"),
                        ci_method = 'delong',
                        conf_level = 0.95,
                        pseudomedian = FALSE,
                        n_grid = 1000,
                        ...) {
  
  alternative <- match.arg(alternative)
  pvalue_method <- match.arg(pvalue_method)
  
  #
  # check
  if (is.null(data)) {
    stop("For formula interface, 'data' must be provided.")
  }
  if (!is.data.frame(data)) {
    stop("For formula interface, 'data' must be a data frame.")
  }
  
  vars <- all.vars(formula)
  if (length(vars) != 2) {
    stop("Formula must be of the form 'response ~ group'.")
  }
  
  response_var <- vars[1]
  group_var <- vars[2]
  
  if (!(response_var %in% names(data))) {
    stop(paste("Variable", response_var, "not found in data."))
  }
  if (!(group_var %in% names(data))) {
    stop(paste("Variable", group_var, "not found in data."))
  }
  
  ############################################################
  #
  # data frame columns
  #
  values <- data[[response_var]]
  groups <- data[[group_var]]
  
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }
  
  group_levels <- levels(groups)
  if (length(group_levels) != 2) {
    stop("Group variable must have exactly 2 levels.")
  }
  
  # Store original group levels and reference level
  group_levels_original <- group_levels
  #
  if (is.null(ref_level) == TRUE) {
    #
    group_ref_level <- group_levels[1]  # First level is reference; ie 'controls'
    #
  } else {
    #
    group_ref_level <- ref_level
    #
  }
  
  
  ##############################################################################
  #
  # wilcox.test() requires x to be 'cases' to get correct AUC
  x_vals <- values[groups != group_ref_level]    # cases
  y_vals <- values[groups ==  group_ref_level]   # controls
  #
  data_name <- paste(response_var, "by", group_var)
  #
  
  
  ##############################################################################
  #
  # Remove missing values
  x_vals <- x_vals[!is.na(x_vals)]
  y_vals <- y_vals[!is.na(y_vals)]
  
  if (length(x_vals) == 0 || length(y_vals) == 0) {
    stop("Groups cannot be empty after removing missing values.")
  }
  
  n1 <- length(x_vals)
  n2 <- length(y_vals)
  
  
  ##############################################################################
  #
  # Perform main Wilcoxon test
  #
  # To get eAUC
  pooled <- c(x_vals, y_vals)
  # compute eAUC
  ranks <- rank(pooled)
  rank_sum <- sum(ranks[1:n1])
  # Convert to Mann-Whitney U
  U_statistic <- rank_sum - n1 * (n1 + 1) / 2
  Ahat <- U_statistic/(n1*n2)
  #
  # p-value
  #
  # NOTE: wmwAUC_pvalue_EU()/wmwAUC_pvalue_BC() internally test H0 in the
  # P(X<Y) convention, whereas this function's Ahat and the "alternative"
  # argument use the P(X>Y) = P(case > reference) convention documented
  # above. Passing (y_vals, x_vals), rather than (x_vals, y_vals), makes
  # the internal P(X<Y) equal P(reference < case) = P(case > reference),
  # matching this function's convention exactly -- so "greater"/"less" mean
  # what the documentation says without needing to also flip the
  # alternative string.
  if (pvalue_method == 'EU') {
    #
    p_value <- wmwAUC_pvalue_EU(y_vals, x_vals, alternative = alternative, A0 = 0.5)
    #
  } else {
    # compute pval
    p_value <- wmwAUC_pvalue_BC(y_vals, x_vals, alternative = alternative, A0 = 0.5)
    #
  }
  #
  # To get correct confint for the pseudomedian
  # by inverting the test
  if (pseudomedian == TRUE) {
    #
    if(n1 < 20 || n2 < 20) {
      warning("Small sample sizes (n1=", n1, ", n2=", n2, ") detected.\n",
              "Computing exact pseudomedian CI requires permutation tests\n",
              "for each grid point (~100 points x 2000 permutations).\n",
              "This may take several minutes.",
              call. = FALSE, immediate. = TRUE)
    }
    #
    pseudo_median <-  wmwAUC_pseudomedian_ci(x_vals, y_vals, conf.level = conf_level,
                                             pvalue_method = pvalue_method, n_grid = n_grid)
    pseudomedian_estimate <- pseudo_median$estimate
    pseudomedian_conf_int <- pseudo_median$conf.int
    pseudomedian_conf_level <- pseudo_median$conf.level
    #
  } else {
    #
    pseudomedian_estimate <-  NULL
    pseudomedian_conf_int <- NULL
    pseudomedian_conf_level <- NULL
    #
  }
  
  
  ##############################################################################
  #
  # Calculate AUC and confidence interval
  #
  # AUC confidence interval
  rc = roc_with_ci(probs = values,
                   labels = groups,
                   positive = group_levels[group_levels != group_ref_level],
                   auc = Ahat,
                   ci_method = ci_method,
                   alpha = 1 - conf_level, ...)
  auc_ci = rc$auc_ci
  
  ##############################################################################
  #
  # Create result object with group level information
  result <- list(
    #
    pseudomedian_requested = pseudomedian,
    n = c(n1 = n1, n2 = n2),
    U_statistic = U_statistic,
    p_value = p_value,
    alternative = alternative,
    pvalue_method = pvalue_method,
    data_name = data_name,
    pseudomedian = pseudomedian_estimate,
    pseudomedian_conf_int = pseudomedian_conf_int,
    pseudomedian_conf_level = pseudomedian_conf_level,
    # eROC, eAUC
    ci_method = ci_method,
    roc_object = rc,
    auc = Ahat, # eAUC
    auc_conf_int = auc_ci,
    #
    x_vals = x_vals,
    y_vals = y_vals,
    #
    groups = groups,
    group_levels = group_levels_original,
    group_ref_level = group_ref_level
  )
  
  class(result) <- "wmwAUC_test"
  return(result)
}