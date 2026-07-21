#' Confidence Interval for the AUC-Equalizing Shift (Hodges-Lehmann Pseudomedian) via Test Inversion
#' 
#' Computes a confidence interval for the pseudomedian by inverting the test
#' \eqn{\mathrm{H_0\colon AUC}(X, Y+\delta) = 0.5}{H0: AUC(X, Y+delta) = 0.5}.
#' 
#' @param x numeric vector, first sample
#' @param y numeric vector, second sample  
#' @param conf.level confidence level (default 0.95)
#' @param pvalue_method character, either 'EU' or 'BC'
#' @param n_grid number of grid points for search (default 1000)
#' @return list with conf.int, estimate and conf.level
#' 
#' @details
#' The pseudomedian \eqn{\delta} is defined generally, for any pair of
#' distributions \eqn{F, G}, as the shift solving
#' \eqn{P(X < Y + \delta) = 0.5}{P(X < Y + delta) = 0.5} --- i.e. the value of
#' \eqn{\delta} that equalizes the AUC between \eqn{X} and the shifted
#' \eqn{Y + \delta}. This definition does not require a location-shift model
#' relating \eqn{F} and \eqn{G}. Under location-shift
#' (\eqn{F(t) = G(t - \Delta)}{F(t) = G(t - Delta)}), this \eqn{\delta}
#' coincides with the classical Hodges-Lehmann pseudomedian and equals the
#' location difference \eqn{\Delta}; this is a special case, not the
#' definition. Outside location-shift, \eqn{\delta} remains well-defined as
#' the AUC-equalizing shift, but should not be read as "the median
#' difference" between the two groups.
#'
#' The confidence interval is obtained by grid search over candidate
#' \eqn{\delta} values, retaining those for which the AUC=0.5 test applied
#' to \eqn{(X, Y+\delta)} is not rejected at level \code{1 - conf.level}.
#' The search range is set heuristically from \eqn{\pm 3} times the robust
#' MAD scale of the pairwise differences around the point estimate; if the
#' true interval extends beyond this range, or if the acceptance region is
#' not contiguous (which can occur in small, imbalanced samples), a warning
#' is issued. See Warnings.
#' 
#' EU and BC pseudomedian CIs are computed on independent calls, and may occasionally 
#' be reported as numerically identical at coarse n_grid even though the underlying
#' p-value curves genuinely differ — both methods' crossing points can land in the 
#' same discrete grid cell.
#'
#' At small \code{min(length(x), length(y))}, BC's near-zero power (see
#' \code{\link{wmwAUC_pvalue_BC}}) can make its pseudomedian interval very
#' wide and largely uninformative rather than merely conservative, since the
#' grid search rarely rejects anywhere in the searched range. EU is
#' recommended in that regime; see the warning issued below.
#'
#' When \code{pvalue_method = 'EU'} and \code{length(x) + length(y) < 20},
#' each grid point's p-value may come from Monte Carlo permutation (see
#' \code{\link{wmwAUC_pvalue_EU}}), making this function's result stochastic
#' at small sample sizes. For reproducible results, call \code{set.seed()}
#' immediately beforehand.
#'
#' @section Warnings:
#' Two diagnostic warnings are issued when the grid search may be
#' unreliable:
#' \itemize{
#'   \item \strong{Non-contiguous acceptance region}: the interval is
#'   reported as \code{[min(accepted), max(accepted)]}, which assumes the
#'   set of non-rejected \eqn{\delta} values is a single contiguous block.
#'   If it is not, the reported interval may be wider than the true
#'   acceptance region and should be inspected directly.
#'   \item \strong{Search range too narrow}: if the p-value has not dropped
#'   below \eqn{\alpha}{alpha} at either edge of the search grid, the true
#'   endpoint may lie outside the searched range; consider widening the
#'   range or increasing \code{n_grid}.
#' }
#'
#' @importFrom stats mad
#' @importFrom stats pnorm
#' @importFrom stats var
#' 
#' @export
wmwAUC_pseudomedian_ci <- function(x, y, conf.level = 0.95, 
                                   pvalue_method = c('EU', 'BC'),
                                   n_grid = 1000) {
  
  # Input validation
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors")
  }
  if (length(x) == 0 || length(y) == 0) {
    stop("x and y must have length > 0")
  }
  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be between 0 and 1")
  }
  pvalue_method <- match.arg(pvalue_method)
  
  alpha <- 1 - conf.level
  
  n1 = length(x)
  n2 = length(y)
  
  # Single upfront warning (not one per grid point -- see suppressWarnings()
  # below) when BC is used at small min(n1,n2): BC's near-zero power there
  # (see wmwAUC_pvalue_BC()) causes the grid search to rarely reject
  # anywhere in the searched range, producing a very wide, largely
  # uninformative pseudomedian CI rather than merely a conservative one.
  if (pvalue_method == 'BC' && min(n1, n2) < 10) {
    warning("min(n1,n2) = ", min(n1, n2), " is small. Simulation found BC's ",
            "power near zero at min(n1,n2) as large as 6, which causes its ",
            "pseudomedian confidence interval to become very wide and largely ",
            "uninformative (in one check, mean width ~16 for a true shift of ",
            "1, at min(n1,n2)=6) rather than merely conservative. The EU ",
            "method is recommended instead at this sample size.", call. = FALSE)
  }
  
  # Compute all pairwise differences to get search range
  all_diffs <- outer(x, y, "-")
  pseudomedian_est <- median(all_diffs)
  
  # Set search range around the estimate
  range_width <- 2 * mad(all_diffs, constant = 1.4826) # robust scale
  if (range_width == 0) range_width <- 1 # fallback for constant differences
  
  delta_min <- pseudomedian_est - 3 * range_width
  delta_max <- pseudomedian_est + 3 * range_width
  
  # Grid of candidate pseudomedian values
  delta_grid <- seq(delta_min, delta_max, length.out = n_grid)
  
  # Test H0: AUC = 0.5 for each shifted sample 
  # 
  n_kernel_evals <- n1 * n2
  total_computations <- n_kernel_evals * n_grid
  if (total_computations > 1e6) {
    warning(sprintf(
      "Large computational load detected: %dx%d = %s kernel evaluations x %d grid points = %s total computations. Consider reducing n_grid for faster execution.",
      n1, n2, 
      format(n_kernel_evals, big.mark = ","),
      n_grid,
      format(total_computations, big.mark = ",")
    ), call. = FALSE, immediate. = TRUE)
  }  
  
  if (pvalue_method == 'BC') { 
    #
    # suppressWarnings(): without this, wmwAUC_pvalue_BC()'s own small-n
    # power warning fires once per grid point (up to n_grid times per call,
    # e.g. thousands of near-identical warnings) -- the single upfront
    # warning above covers this instead, matching the EU branch below
    # which already suppressed its own per-grid-point warnings.
    p_values <- sapply(delta_grid, function(delta) {
      y_shifted <- y + delta
      suppressWarnings(wmwAUC_pvalue_BC(x, y_shifted, alternative = "two.sided"))
    })
    #
  } else { 
    #
    p_values <- sapply(delta_grid, function(delta) {
      y_shifted <- y + delta
      suppressWarnings(wmwAUC_pvalue_EU(x, y_shifted, alternative = "two.sided"))
    })
    #
  }
  
  # Find confidence interval endpoints
  accepted_idx <- which(p_values >= alpha)
  accepted <- delta_grid[accepted_idx]
  
  if (length(accepted) == 0) {
    warning("No accepted values found. Try increasing n_grid or check data.")
    ci <- c(NA, NA)
  } else {
    ci <- c(min(accepted), max(accepted))
    
    # Diagnostic: acceptance region should be a single contiguous block.
    # A gap indicates the [min, max] interval may overstate the true
    # (possibly disjoint) acceptance region.
    if (length(accepted_idx) > 1 && !all(diff(accepted_idx) == 1)) {
      warning("Acceptance region is not contiguous; reported CI may be unreliable. Consider inspecting p_values across delta_grid directly.", call. = FALSE)
    }
    
    # Diagnostic: search range should bracket the full CI, i.e. the
    # p-value should have dropped below alpha before reaching either edge
    # of delta_grid. If not, the true endpoint may lie outside the range.
    if (p_values[1] >= alpha || p_values[n_grid] >= alpha) {
      warning("Search range may not bracket the full confidence interval; consider widening the range (currently +/- 3*MAD) or increasing n_grid.", call. = FALSE)
    }
  }
  
  list(
    conf.int = ci,
    estimate = pseudomedian_est,
    conf.level = conf.level
  )
}
