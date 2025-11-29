#' Wilcoxon-Mann-Whitney Test of No Group Discrimination  (Continuous Variables)
#'
#' Performs distribution-free Wilcoxon-Mann-Whitney test for AUC-detectable 
#' group discrimination, testing \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5}  
#' against \eqn{\mathrm{H_1\colon AUC} \neq 0.5}{H1: AUC != 0.5}.
#' Under location-shift assumption, equivalently tests zero location difference.
#' 
#' Constructs confidence intervals for the pseudomedian via test inversion.
#' Under location-shift assumptions (\eqn{G(x) = F(x - \delta)}), the pseudomedian
#' represents the location difference between groups.
#' 
#' Derived under \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5}.
#'
#' @note Current implementation assumes continuous data without ties.
#'
#'
#' @param formula Formula of the form `response ~ group`
#' @param data Data frame containing continuous response variable and grouping factor
#' @param ref_level Character, reference level of grouping factor (if NULL, uses first level)
#' @param special_case Logical, location-shift assumption (default FALSE)
#' @param alternative Character, alternative hypothesis is c("two.sided", "greater", "less")
#' @param ci_method Character, confidence interval method for eAUC: c('hanley', 'boot', 'none')
#' @param conf_level Numeric, confidence level for intervals (default 0.95)
#' @param exact Logical, passed to `wilcox.test()`  for exact p-values (default NULL)
#' @param correct Logical, passed to `wilcox.test()` for continuity correction (default TRUE)
#' @param tol_root Numeric, passed to `wilcox.test()` for root tolerance (default 1e-4)
#' @param digits_rank Numeric, passed to `wilcox.test()` for ranking precision (default Inf)
#' @param ... Additional arguments passed to `roc_with_ci()`
#'
#' @details
#' The function tests the null hypothesis \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5} 
#' against \eqn{\mathrm{H_0\colon AUC} \neq 0.5}{H1: AUC != 0.5}, 
#' where AUC represents the Area Under the ROC Curve and equals the probability 
#' \eqn{P(X > Y)} that a randomly selected observation from the first group exceeds a randomly 
#' selected observation from the second group.
#' 
#' For `response ~ group`, observations from the non-reference group constitute X, 
#' while observations from the reference group (specified by `ref_level`) constitute Y.
#' Thus AUC = P(non-reference > reference). If `ref_level` is not specified, the first 
#' factor level is used as reference. The W statistic from \code{\link{wilcox.test}} and the 
#' resulting empirical AUC (eAUC) are calculated consistently with this group assignment.
#'
#' The test statistic is eAUC, which estimates the true AUC. 
#' The empirical ROC curve (eROC) is constructed by varying the classification 
#' threshold across all observed values and computing sensitivity and 1-specificity 
#' at each threshold.
#' 
#' When special_case = TRUE, the function additionally reports location-shift 
#' parameters under the assumption that \eqn{F_1(x) = F_2(x - \delta)}{F1(x) = F2(x - delta)}. 
#' Under this assumption, the discrimination test \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5} is mathematically 
#' equivalent to testing H0: \eqn{\delta = 0}{delta = 0} (zero location shift). 
#' In this special case, eAUC takes the dual role of both test statistic and effect size 
#' for the location difference.
#' 
#' Confidence intervals for the true AUC are computed using either the Hanley and 
#' McNeil (1982) method based on asymptotic normality, or bootstrap resampling.
#' If bootstrap resampling is selected, it is also used for constructing the 
#' confidence band for the ROC curve.
#' 
#' **Statistical Methodology:**
#' Unlike standard implementations that assume the erroneously broad null hypothesis 
#' \eqn{\mathrm{H_0\colon F = G}}{H0: F = G}, 
#' this function derives p-values under the correct null hypothesis 
#' \eqn{\mathrm{H_0\colon AUC} = 0.5}{H0: AUC = 0.5} 
#' that WMW actually tests. P-values are computed using asymptotic distribution 
#' theory with sample-size dependent bias corrections to maintain proper Type I 
#' error control under heteroskedasticity. Confidence intervals for the pseudomedian
#' are obtained by inverting the test.  
#' 
#'
#' @return Object of class 'wmw_test' containing:
#'   \item{special_case}{Logical indicating whether special case (location-shift) analysis was performed}
#'   \item{n}{Named vector with components n1, n2 giving sample sizes for each group}
#'   \item{statistic}{W statistic returned by `wilcox.test`}
#'   \item{p.value}{P-value for testing H0: AUC = 0.5}
#'   \item{alternative}{Alternative hypothesis specification used in `wilcox.test`}
#'   \item{method}{Character string describing the test method}
#'   \item{data.name}{Character string giving the name of the data}
#'   \item{estimate}{Hodges-Lehmann median difference estimate (when special_case = TRUE)}
#'   \item{conf.int}{Confidence interval for the location shift (when special_case = TRUE)}
#'   \item{conf.level}{Confidence level for the confidence interval for HL estimator (when special_case = TRUE)}
#'   \item{null.value}{Null hypothesis value (mu = 0 for location shift)}
#'   \item{ci_method}{Method used to compute confidence interval for AUC}
#'   \item{roc_object}{ROC analysis object returned by `roc_with_ci` function}
#'   \item{auc}{Empirical AUC (eAUC) calculated from Wilcoxon W statistic}
#'   \item{auc.conf.int}{Confidence interval for true AUC using Hanley-McNeil or bootstrap method}
#'   \item{x_vals}{Numeric vector of observations from non-reference group}
#'   \item{y_vals}{Numeric vector of observations from reference group}
#'   \item{groups}{Character vector of group labels from original data}
#'   \item{group_levels}{Character vector of factor levels for grouping variable}
#'   \item{group_ref_level}{Character string indicating which level corresponds to reference group}
#'
#' @examples
#'
#' #############################################################################
#' #
#' library(wmwAUC)  
#' #
#' #############################################################################
#'
#' #############################################################################
#' #
#' # Simulation 1: H0: F=G is erroneously too broad
#' #
#' #############################################################################
#' #
#' # This simulation takes several minutes to complete
#' if (FALSE)  {
#' N = 10000
#' n = 1000
#' set.seed(123L)
#' pval_wt = pval_wmw = eauc = numeric(N)
#' for (i in 1:N) {
#' #
#'   x = rnorm(n, sd = 0.1)
#'   y = rnorm(n, sd = 3)
#'   # wilcox.test() of H0: F = G
#'   wt = wilcox.test(x, y)
#'   pval_wt[i] = wt$p.value
#'   # wmw_test() of H0: AUC = 0.5
#'   pval_wmw[i] = wmw_pvalue(x, y)
#'   # eAUC
#'   eauc[i] = wt$statistic/(n*n)
#'   #
#' }
#' }
#' data(simulation1)  # List eauc, pval_wt, pval_wmw
#' # Empirical AUC centered at 0.5 despite F != G
#' hist(simulation1$eauc)
#' # Traditional p-values (incorrectly testing F = G)
#' hist(simulation1$pval_wt) 
#' # Correct p-values (correctly testing AUC = 0.5)  
#' hist(simulation1$pval_wmw)
#'  
#' #############################################################################
#' #
#' # Simulation 2: H1: F stoch. dominates G is erroneously too narrow
#' #               WMW is consistent for broader H1: AUC != 0.5   
#' #
#' #############################################################################
#' #
#' # This simulation takes several minutes to complete
#' if (FALSE) {
#' N = 10000
#' n = 1000
#' set.seed(123L)
#' pval_wt = pval_wmw = eauc = numeric(N)
#' for (i in 1:N) {
#' #
#'   # gaussians with different location and scale
#'   # does not satisfy stochastic dominance
#'   x = rnorm(n, 0, sd = 0.1)
#'   y = rnorm(n, 0.5, sd = 3)
#'   # wilcox.test H0: F = G vs H1: (F stochastically dominates G) OR (G stochastically dominates F)
#'   wt = wilcox.test(x, y)
#'   pval_wt[i] = wt$p.value
#'   # wmw_test H0: AUC = 0.5 vs H1: AUC neq 0.5
#'   pval_wmw[i] = wmw_pvalue(x, y)
#'   # eAUC
#'   eauc[i] = wt$statistic/(n*n)
#' #
#' }
#' }
#'
#' data(simulation2)  # List of eauc, pval_wt, pval_wmw
#' # WMW detects broader alternatives than traditional stochastic dominance
#' hist(simulation2$eauc)
#' hist(simulation2$pval_wmw)
#' hist(simulation2$pval_wt)
#'
#' #############################################################################
#' #
#' # Simulation 3: confidence interval for pseudomedian derived under H0: AUC = 0.5   
#' #               MC study of N = 500 replicas
#' #               x ~ rnorm(300, 0,1)
#' #               y ~ rlaplace(300, 0,1)
#' #
#' #############################################################################
#' #
#' # This simulation takes long time to complete
#' if (FALSE) {
#' 
#' if (!requireNamespace("VGAM", quietly = TRUE)) {
#' install.packages("VGAM")
#' }
#' library('VGAM')
#
#' N <- 500
#' n_test <- 300
#'
#' set.seed(123L)
#' wmw_ci = wt_ci = list(N)
#' eauc = pseudomed = numeric(N)
#' for (i in 1:N) {
#'  #
#'  x_test <- rnorm(n_test, 0, 1)
#'  y_test <- VGAM::rlaplace(n_test, 0, 1)
#'
#'  wmw_test <- pseudomedian_ci(x_test, y_test, conf.level = 0.95)
#'  wmw_ci[[i]] = wmw_test$conf.int
#'  wt_test <- wilcox.test(x_test, y_test, conf.int = TRUE)
#'  wt_ci[[i]] = wt_test$conf.int
#'  eauc[i] = wt_test$statistic/(n_test*n_test)
#'  pseudomed[i] = as.numeric(wt_test$estimate)
#'  #  
#' }
#'  
#' wmw_ci = do.call(rbind, wmw_ci)
#' wt_ci = do.call(rbind, wt_ci)
#' colMeans(wmw_ci)
#' colMeans(wt_ci)
#' mean(eauc)
#' 
#' # coverage
#' length(which((wmw_ci[,1] < 0) & (wmw_ci[,2] > 0)))
#' length(which((wt_ci[,1] < 0) & (wt_ci[,2] > 0)))
#'
#' mean(pseudomed)
#' }
#'
#' data(simulation3)  # List of wmw_ci, wt_ci, eauc, pseudomedian
#' # 
#' # Average across MC of confidence intervals obtained under H0: AUC=0.5
#' colMeans(simulation3$wmw_ci)
#' # Average across MC of confidence intervals from wilcox.test()
#' colMeans(simulation3$wt_ci)
#' # 
#' # Average across MC of eAUC
#' mean(simulation3$eauc)
#' 
#' # Coverage
#' length(which((simulation3$wmw_ci[,1] < 0) & (simulation3$wmw_ci[,2] > 0)))
#' length(which((simulation3$wt_ci[,1] < 0) & (simulation3$wt_ci[,2] > 0)))
#' 
#' # Mean pseudomedian
#' mean(simulation3$pseudomed)
#'
#'
#' #############################################################################
#' #
#' # Ex 1
#' #
#' #############################################################################
#' #
#' if (!requireNamespace("gemR", quietly = TRUE)) {
#' install.packages("gemR")
#' }
#' library('gemR')
#'
#' data(MS)
#' da <- MS
#'
#' # preparing data frame
#' class(da$proteins) <- setdiff(class(da$proteins), "AsIs")
#' df <- as.data.frame(da$proteins)
#' df$MS <- da$MS
#'
#' # WMW test 
#' wmd <- wmw_test(P19099 ~ MS, data = df, ref_level = 'no')
#' wmd
#' plot(wmd)
#'
#' # EDA to assess location shift assumption validity
#' qp <- quadruplot(P19099 ~ MS, data = df, ref_level = 'no')
#' qp
#' # => location shift assumption is not valid
#'
#' #############################################################################
#' #
#' # Ex 2
#' #
#' #############################################################################
#' #
#' data(Ex2)
#' da <- Ex2
#'
#' # WMW test
#' wmd <- wmw_test(y ~ group, data = da, ref_level = 'control')
#' wmd
#' plot(wmd)
#'
#' # Check location-shift assumption with EDA
#' qp <- quadruplot(y ~ group, data = da, ref_level = 'control')
#' qp
#' # => location-shift assumption not tenable
#'
#' # Note that medians are essentially the same:
#' median(da$y[da$group == 'case'])
#' # 0.495
#' median(da$y[da$group == 'control'])
#' # 0.493
#'
#' # Erroneous use of location-shift special case would falsely 
#' # conclude significant median difference despite identical medians
#' wml <- wmw_test(y ~ group, data = da, special_case = TRUE,
#'                 ref_level = 'control')
#' wml
#'
#'
#' #############################################################################
#' #
#' # Ex 3
#' #
#' #############################################################################
#' #
#' if (!requireNamespace("gss", quietly = TRUE)) {
#' install.packages("gss")
#' }
#' library('gss')
#'
#' data(wesdr)
#' da = wesdr
#' da$ret = as.factor(da$ret)
#'
#' #############################################################################
#' # WARNING: data with ties
#' #          current version of wmwAUC does not take ties into account
#' #############################################################################
#' # WMW 
#' wmd <- wmw_test(bmi ~ ret, data = da, ref_level = '0')
#' wmd
#' plot(wmd)
#'
#'
#' # EDA to assess location shift assumption validity
#' qp <- quadruplot(bmi ~ ret, data = da, ref_level = '0')
#' qp
#' # => location shift assumption is tenable
#' 
#' # Special case of WMW test
#' suppressWarnings({ # ties in data
#' wml <- wmw_test(bmi ~ ret, data = da, ref_level = '0', 
#'                 ci_method = 'boot', special_case = TRUE)
#' })                 
#' wml
#' plot(wml)
#'
#' 
#' 
#'
#' @references
#' 
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods. Biometrics Bulletin, 
#' 1(6), 80-83.
#' 
#' Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random 
#' variables is stochastically larger than the other. The Annals of Mathematical 
#' Statistics, 18(1), 50-60.
#' 
#' Van Dantzig, D. (1951). On the consistency and the power of Wilcoxon's two sample
#' test. Proceedings KNAW, Series A, 54(1), 1-8.
#' 
#' Birnbaum, Z. W. (1956). On a use of the Mann-Whitney statistic. 
#' In Proceedings of the Third Berkeley Symposium on Mathematical Statistics and 
#' Probability, Volume 1: Contributions to the Theory of Statistics (Vol. 3, pp. 13-18). 
#' University of California Press.
#' 
#' Bamber, D. (1975). The area above the ordinal dominance graph 
#' and the area below the receiver operating characteristic graph. 
#' Journal of mathematical psychology, 12(4), 387-415.
#' 
#' Lehmann, E. L., & Abrera, H. B. D. (1975). Nonparametrics. Statistical methods 
#' based on ranks. San Francisco, CA, Holden-Day.
#' 
#' Hanley, J. A., & McNeil, B. J. (1982). The meaning and use of the area under 
#' a receiver operating characteristic (ROC) curve. Radiology, 143(1), 29-36.
#'
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. 
#' Psychological bulletin, 114(3), 494.
#' 
#' Arcones, M. A., Kvam, P. H., & Samaniego, F. J. (2002). Nonparametric estimation 
#' of a distribution subject to a stochastic precedence constraint. Journal of the 
#' American Statistical Association, 97(457), 170-182.
#' 
#' Pepe, M. S. (2003). The statistical evaluation of medical tests for classification 
#' and prediction. Oxford university press.
#' 
#' Conroy, R. M. (2012). What hypotheses do “nonparametric” two-group tests actually test?. 
#' The Stata Journal, 12(2), 182-190.
#' 
#' del Barrio, E., Cuesta-Albertos, J. A., & Matrán, C. (2025). Invariant measures 
#' of disagreement with stochastic dominance. The American Statistician, 1-13.
#' 
#' Grendar, M. (2025). Wilcoxon-Mann-Whitney test of no group discrimination. arXiv:2511.20308.
#'
#' @seealso
#' \code{\link{print.wmw_test}} for formated output of `wmw_test()`.
#' \code{\link{plot.wmw_test}} for plot of output of  `wmw_test()`.
#' \code{\link{wmw_pvalue}} for details on computing p-values
#' \code{\link{pseudomedian_ci}} for details on computing confidence intervals for pseudomedian
#' \code{\link{quadruplot}} for exploratory data analysis plots that assist in evaluating location-shift assumption validity.
#' \code{\link{wilcox.test}} for the underlying Wilcoxon-Mann-Whitney test.
#'
#' @export
wmw_test <- function(formula,
                     data,
                     ref_level = NULL,
                     special_case = FALSE,
                     alternative = c("two.sided", "greater", "less"),
                     ci_method = 'hanley',
                     conf_level = 0.95,
                     exact = NULL,
                     correct = TRUE,
                     tol_root = 1e-4,
                     digits_rank = Inf,
                     ...) {

  alternative <- match.arg(alternative)

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
  if (is.null(ref_level) == T ) {
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
  wtest <- stats::wilcox.test(x_vals, y_vals, alternative = alternative,
                       mu = 0,
                       conf.int = TRUE, conf.level = conf_level,
                       exact = exact, correct = correct, tol.root = tol_root,
                       digits.rank = digits_rank)
  # To get the correct p-value
  pval = wmw_pvalue(x_vals, y_vals, alternative = alternative)
  wtest$p.value = pval
  
  # To get correct confint for HL pseudomedian
  # by inverting the test
  p_median =  pseudomedian_ci(x_vals, y_vals, conf.level = conf_level, n_grid = 1000) 
  wtest$conf.int = p_median$conf.int
  wtest$estimate = p_median$estimate
  wtest$conf.level = p_median$conf.level


  ##############################################################################
  #
  # Calculate AUC and confidence interval
  #
  W <- wtest$statistic
  auc <- as.numeric(W) / (n1 * n2)
  #
  # AUC confidence interval
  rc = roc_with_ci(probs = values,
           labels = groups,
           positive = group_levels[group_levels != group_ref_level],
           auc = auc,
           ci_method = ci_method,
           alpha = 1 - conf_level, ...)
  auc_ci = rc$auc_ci

  ##############################################################################
  #
  # Create result object with group level information
  result <- list(
    # 
    special_case = special_case,
    n = c(n1 = n1, n2 = n2),
    statistic = wtest$statistic,
    p.value = wtest$p.value,
    alternative = wtest$alternative,
    method = wtest$method,
    data.name = data_name,
    estimate = wtest$estimate,
    conf.level = wtest$conf.level,
    conf.int = wtest$conf.int,
    null.value = wtest$null.value,
    # eROC, eAUC
    ci_method = ci_method,
    roc_object = rc,
    auc = auc,
    auc.conf.int = auc_ci,
    #
    x_vals = x_vals,
    y_vals = y_vals,
    #
    groups = groups,
    group_levels = group_levels_original,
    group_ref_level = group_ref_level
  )

  class(result) <- "wmw_test"
  return(result) 
}


