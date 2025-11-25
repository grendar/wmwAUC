#' Adds simultaneous confidence band to ECDF using sfsmisc
#'
#' @param p `ggcdf` ggplot returned by `quadruplot()`
#' @param data data frame used in `quadruplot()`
#' @param response_col character giving the name of response variable
#' @param group_col character giving the name of group factor
#' @param ref_level character giving the reference level of group factor
#' @param alpha size of test (0.05) used to provide confidence level
#'
#' @keywords internal
#' @export
add_simultaneous_bands_sfsmisc <- function(p, data, response_col, group_col, ref_level = NULL, alpha = 0.05) {


  if (is.null(ref_level)) {
    ref_level <- levels(factor(data[[group_col]]))[1]
  }
  
  # APPLY THE SAME MEDIAN-CENTERING AS QUADRUPLOT
  data_centered <- data
  groups <- unique(data[[group_col]])
  for (g in groups) {
    idx <- data[[group_col]] == g
    data_centered[idx, response_col] <- data[idx, response_col] - median(data[idx, response_col], na.rm = TRUE)
  }
  
  # Calculate bands using MEDIAN-CENTERED data and sfsmisc
  all_bands <- data.frame()
  
  for (g in groups) {
    group_data <- data_centered[data_centered[[group_col]] == g, response_col]
    group_data <- group_data[!is.na(group_data)]
    
    # Calculate simultaneous bands using sfsmisc::approx.ksD
    bands <- calc_simultaneous_ecdf_bands_sfsmisc(group_data, alpha)
    bands$group <- g
    bands$is_ref <- (g == ref_level)
    
    all_bands <- rbind(all_bands, bands)
  }
  
  # Define colors to match ECDF lines
  colors <- c("TRUE" = "blue", "FALSE" = "red")
  
  # Add simultaneous confidence bands to the plot
  p <- p + 
    ggplot2::geom_ribbon(data = all_bands,
                ggplot2::aes(x = x, ymin = lower, ymax = upper, fill = is_ref),
                alpha = 0.15, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = colors, guide = "none")
  
  return(p)
}





###############################################################################
#'
#'
#' Confidence bands for ECDF using sfsmisc::KSd 
#'
#' @param x numeric vector
#' @param alpha size of test; hence confidence level is 1 - alpha
#' 
#' @importFrom stats ecdf pt
#'
#' @keywords internal
#' @export
calc_simultaneous_ecdf_bands_sfsmisc <- function(x, alpha = 0.05) {
  
  if (!requireNamespace("sfsmisc", quietly = TRUE)) {
    stop("Package 'sfsmisc' is required for ECDF bands. Please install it.")
  }
  
  # Note: KSd() is hardcoded for alpha = 0.05 (95% confidence)
  if (alpha != 0.05) {
    warning("sfsmisc::KSd() only supports alpha = 0.05 (95% confidence). Using alpha = 0.05.")
    alpha <- 0.05
  }
  
  n <- length(x)
  
  # This exactly mirrors what ecdf.ksCI() does internally:
  ec <- ecdf(x)
  xx <- get("x", envir = environment(ec))  # x coordinates (sorted unique values)
  yy <- get("y", envir = environment(ec))  # ECDF values at those points
  
  # Use sfsmisc's authoritative KSd calculation (95% only)
  D <- sfsmisc::KSd(n)
  
  # Calculate simultaneous confidence bands
  upper <- pmin(yy + D, 1)
  lower <- pmax(yy - D, 0)
  
  return(data.frame(
    x = xx,
    ecdf = yy,
    upper = upper,
    lower = lower
  ))
}

