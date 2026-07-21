#' Deprecated Function Names (wmwAUC 1.0.0)
#'
#' @description These functions have been renamed in wmwAUC 1.0.0, alongside
#' argument and behavior changes in some cases. Calling them raises an
#' informative error pointing to the current function and to NEWS.md, rather
#' than attempting to forward the call -- for functions with argument-name or
#' argument-value changes (see Details), a silent passthrough could otherwise
#' produce a different result than the old function used to, without any
#' indication that anything had changed.
#'
#' @param ... Not used; calling any of these functions always raises an error.
#'
#' @details
#' \itemize{
#'   \item \code{wmw_test()} -> \code{\link{wmwAUC_test}}. The \code{special_case}
#'   argument was renamed \code{pseudomedian}, and the \code{ci_method = 'hanley'}
#'   option was removed (not renamed), and replaced by the DeLong et al. (1988) method.
#'   \item \code{wmw_pvalue()} -> \code{\link{wmwAUC_pvalue_BC}}.
#'   \item \code{wmw_pvalue_ties()} -> \code{\link{wmwAUC_pvalue_EU}}.
#'   \item \code{pseudomedian_ci()} -> \code{\link{wmwAUC_pseudomedian_ci}}.
#' }
#'
#' @name wmwAUC-deprecated
NULL

#' @rdname wmwAUC-deprecated
#' @export
wmw_test <- function(...) {
  stop("wmw_test() has been renamed to wmwAUC_test() in wmwAUC 1.0.0, with ",
       "argument changes (special_case -> pseudomedian) and a removed ",
       "ci_method option ('hanley', replaced by 'delong'). Please update your ",
       "code to call wmwAUC_test() directly; see NEWS.md for details.",
       call. = FALSE)
}

#' @rdname wmwAUC-deprecated
#' @export
wmw_pvalue <- function(...) {
  stop("wmw_pvalue() has been renamed to wmwAUC_pvalue_BC() in wmwAUC 1.0.0. ",
       "Please update your code to call wmwAUC_pvalue_BC() directly; see ",
       "NEWS.md for details.",
       call. = FALSE)
}

#' @rdname wmwAUC-deprecated
#' @export
wmw_pvalue_ties <- function(...) {
  stop("wmw_pvalue_ties() has been renamed to wmwAUC_pvalue_EU() in wmwAUC ",
       "1.0.0. Please update your code to call wmwAUC_pvalue_EU() directly; ",
       "see NEWS.md for details.",
       call. = FALSE)
}

#' @rdname wmwAUC-deprecated
#' @export
pseudomedian_ci <- function(...) {
  stop("pseudomedian_ci() has been renamed to wmwAUC_pseudomedian_ci() in ",
       "wmwAUC 1.0.0. Please update your code to call ",
       "wmwAUC_pseudomedian_ci() directly; see NEWS.md for details.",
       call. = FALSE)
}