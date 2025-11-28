#' @name simulation1
#' @title Synthetic data
#' @docType data
#' @usage data(simulation1)
#' @format A list containing simulation results (N=10000, n=1000):
#' \describe{
#'   \item{eauc}{Empirical AUC values}
#'   \item{pval_wt}{Traditional wilcox.test p-values}
#'   \item{pval_wmw}{WMW p-values under H0: AUC = 0.5}
#' }
NULL
#'
#' @name simulation2
#' @title Synthetic data
#' @docType data
#' @usage data(simulation2)
#' @format A list containing simulation results (N=10000, n=1000):
#' \describe{
#'   \item{eauc}{Empirical AUC values}
#'   \item{pval_wt}{Traditional wilcox.test p-values}
#'   \item{pval_wmw}{WMW p-values under H0: AUC = 0.5}
#' }
NULL
#'
#'
#' @name simulation3
#' @title Synthetic data
#' @docType data
#' @usage data(simulation3)
#' @format A list containing simulation results (N=500, n=300):
#' \describe{
#'   \item{wmw_ci}{95% confidence intervals obtained by pseudomedian_ci()}
#'   \item{wt_ci}{95% confidence intervals obtained by wilcox.test()}
#'   \item{eauc}{Values of eAUC}
#'   \item{pseudomedian}{Values of the pseudomedian}
#' }
NULL
#'
#' @name Ex2
#' @title Synthetic data
#' @description A data frame with numeric `y` and factor `group`
#' @docType data
#' @usage data(Ex2)
#' @format A data frame with 200 observations on 2 variables.
NULL
