#' Four EDA Plots for Visual Assessment of Location-Shift Assumption
#'
#' Creates four diagnostic plots to visually assess whether the location-shift 
#' assumption \eqn{F_1(x) = F_2(x - \delta)}{F1(x) = F2(x - delta)} holds: 
#' (1) boxplot with swarmplot overlay, 
#' (2) density plot comparison, (3) wormplot of median-centered residuals, and 
#' (4) empirical CDF comparison with confidence band for median-centered data.
#'
#'
#' @param formula Formula of the form `response ~ group`
#' @param data Data frame containg response, group
#' @param ref_level Character, reference level of the grouping factor. If NULL (default), 
#'   uses first factor level
#' @param test Character, statistical test for shift-equivalence assumption. 
#'   Tests for distributional equality applied to median-centered data: 
#'   "ks" (Kolmogorov-Smirnov) (default), "kuiper" (Kuiper), "cvm" (Cramér-von Mises), 
#'   "ad" (Anderson-Darling), "wass" (Wasserstein), "dts" (DTS test).
#' @param ylab Character, label for y-axis. If NULL (default), uses variable name
#' @param color_palette Character, color palette to use. One of "viridis", "plasma", 
#'   "inferno", "magma", or "cividis"
#' @param combine_plots Logical, whether to return combined plot using patchwork 
#'   (TRUE) or list of individual plots (FALSE)
#' @param distribution Character, theoretical distribution for Q-Q plot comparison. 
#'   Default is "norm" for normal distribution
#' @param show_colors Logical, whether to use colors (TRUE) or grayscale (FALSE)
#' @param show_legend Logical, whether to display legend in plots (default TRUE)
#'
#' @return If combine_plots = TRUE, returns a combined ggplot object created by 
#'   patchwork. If FALSE, returns a list of four ggplot objects named 'boxplot', 
#'   'density', 'wormplot', and 'ecdf'.
#'
#' @details
#' The location-shift assumption is assessed by applying a test of H0: equality 
#' of distributions to median-centered data. One of the tests from the `twosamples` 
#' package can be used. The empirical CDF plot includes 95% confidence bands for 
#' the difference between distributions, computed using the `sfsmisc::KSd` function 
#' based on the Kolmogorov-Smirnov distribution. These bands help assess whether 
#' observed differences between median-centered distributions exceed what would be 
#' expected under the location-shift assumption.
#' 
#' 
#' @note Uses \pkg{twosamples} for distribution comparison and 
#' \code{KSd} from \pkg{sfsmisc} for exact confidence bands.
#'   
#'
#' @examples
#'
#' library(wmwAUC)
#' 
#' data(Ex2)
#' da <- Ex2
#' qp = quadruplot(y ~ group, data = da, ref_level = 'control')
#' qp
#'
#'
#' @references 
#' 
#'   O'Dowd, C. (2025). Statistical Code Examples. 
#'   \url{https://codowd.com/code} (accessed November 28, 2025).
#'   
#'    Maechler M (2024). _sfsmisc: Utilities from 'Seminar fuer Statistik' ETH
#'    Zurich_. R package version 1.1-20,
#'    <https://CRAN.R-project.org/package=sfsmisc>.
#'
#' @export
quadruplot <- function(formula, data, 
                       ref_level = NULL, 
                       test = 'ks', 
                       ylab = NULL,
                       color_palette = "lancet",
                       combine_plots = TRUE,
                       distribution = "norm",
                       show_colors = TRUE,
                       show_legend = TRUE) {

  test <- match.arg(test)
  
  # Input validation
  vars <- all.vars(formula)
  if (length(vars) != 2) {
    stop("Formula must be of the form 'response ~ group'.")
  }
  
  what <- vars[1]
  by <- vars[2]
  
  if (!(what %in% names(data))) {
    stop(paste("Variable", what, "not found in data."))
  }
  if (!(by %in% names(data))) {
    stop(paste("Variable", by, "not found in data."))
  }
  

  # Check required packages
  required_pkgs <- c("ggplot2", "rlang", "ggbeeswarm", "qqplotr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Required packages missing: ", paste(missing_pkgs, collapse = ", "))
  }

  if (combine_plots && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required when combine_plots = TRUE")
  }

  # Convert grouping variable to factor if needed
  data[[by]] <- as.factor(data[[by]])
  data[[by]] <- relevel(data[[by]], ref = ref_level)

  # Set up color scales (both color and fill need to match)
  #
  # Extract Lancet palette colors (default number of colors)
  lancet_colors <- ggsci::pal_lancet()(9)
  #
  color_scale_color <- switch(color_palette,
                              "lancet" = if (requireNamespace("ggsci", quietly = TRUE)) ggplot2::scale_color_manual(values = lancet_colors) else ggplot2::scale_color_discrete(),
                              "viridis" = ggplot2::scale_color_viridis_d(),
                              "brewer" = ggplot2::scale_color_brewer(type = "qual"),
                              "default" = ggplot2::scale_color_discrete(),
                              ggplot2::scale_color_discrete()  # fallback
  )

  if (show_colors) {
    color_scale_fill <- switch(color_palette,
                               "lancet" = if (requireNamespace("ggsci", quietly = TRUE)) ggplot2::scale_fill_manual(values = lancet_colors) else ggplot2::scale_fill_discrete(),
                               "viridis" = ggplot2::scale_fill_viridis_d(),
                               "brewer" = ggplot2::scale_fill_brewer(type = "qual"),
                               "default" = ggplot2::scale_fill_discrete(),
                               ggplot2::scale_fill_discrete()  # fallback
    )
  }

  # Base ggplot for consistency
  ggp_base <- ggplot2::ggplot(data) +
    ggplot2::theme_classic()

  ##############################################################################
  #
  # Boxplot w beeswarm
  #
  ##############################################################################
  #
  # 1. Boxplot with beeswarm
  if (show_colors) {
    ggb <- ggp_base +
      ggplot2::aes(x = "", y = !!rlang::sym(what), color = !!rlang::sym(by)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom() +
      ggplot2::facet_wrap(rlang::sym(by)) +
      ggplot2::labs(y = ylab) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     legend.position = if(show_legend) "right" else "none") +
      color_scale_color

  } else {
    ggb <- ggp_base +
      ggplot2::aes(x = "", y = !!rlang::sym(what)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom() +
      ggplot2::facet_wrap(rlang::sym(by)) +
      ggplot2::labs(y = ylab) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  ggb <- ggb + ggplot2::labs(title = 'Boxplot with Swarmplot') +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())


  ##############################################################################
  #
  # Density plot
  #
  ##############################################################################
  #
  # 2. Density plot
  if (show_colors) {
    ggd <- ggp_base +
      ggplot2::aes(x = !!rlang::sym(what), color = !!rlang::sym(by), fill = !!rlang::sym(by)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::facet_wrap(rlang::sym(by)) +
      ggplot2::labs(y = 'Density', x = '') +
      ggplot2::theme(legend.position = if(show_legend) "right" else "none") +
      color_scale_color +
      color_scale_fill
  } else {
    ggd <- ggp_base +
      ggplot2::aes(x = !!rlang::sym(what)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::facet_wrap(rlang::sym(by)) +
      ggplot2::labs(y = 'Density')
  }
  ggd <- ggd + ggplot2::labs(title = 'Density Plot')  +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_blank())


  ##############################################################################
  #
  # Wormplot of median centered data (both groups in one plot)
  #
  ##############################################################################
  #

  # Helper function to create median-centered wormplot without confidence bands
  create_centered_wormplot <- function(data, what, by, distribution, show_colors, color_scale_color, show_legend) {
    # Temporarily suppress warnings for qqplotr
    old_warn <- getOption("warn")
    options(warn = -1)
    on.exit(options(warn = old_warn))

    # Center both groups by their medians
    data_centered <- data
    groups <- unique(data[[by]])
    for (g in groups) {
      idx <- data[[by]] == g
      data_centered[idx, what] <- data[idx, what] - median(data[idx, what], na.rm = TRUE)
    }

    if (show_colors) {
      p <- ggplot2::ggplot(data_centered, ggplot2::aes(sample = !!rlang::sym(what), color = !!rlang::sym(by))) +
        qqplotr::stat_qq_line(distribution = distribution, detrend = TRUE) +
        qqplotr::stat_qq_point(distribution = distribution, detrend = TRUE) +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = if(show_legend) "bottom" else "none") +
        color_scale_color #+
        #ggplot2::labs(subtitle = "Similar patterns suggest location shift")
    } else {
      p <- ggplot2::ggplot(data_centered, ggplot2::aes(sample = !!rlang::sym(what))) +
        qqplotr::stat_qq_line(distribution = distribution, detrend = TRUE) +
        qqplotr::stat_qq_point(distribution = distribution, detrend = TRUE) +
        ggplot2::theme_classic() +
        ggplot2::facet_wrap(rlang::sym(by)) +
        ggplot2::labs(subtitle = "Median-Centered Data")
    }
    return(p)
  }

  # 3. Median-centered wormplot (location shift assessment)
  ggw_centered <- create_centered_wormplot(data, what, by, distribution, show_colors, color_scale_color, show_legend) +
    ggplot2::labs(x = "Theoretical Quantiles", y = "Detrended Sample Quantiles",
                  title = 'Wormplots of Median-Centered Data')# "Location Shift Assessment")


  ##############################################################################
  #
  # 4. Empirical CDFs of median centered (both groups in one plot)
  #
  ##############################################################################
  #
  # Center both groups by their medians
  data_centered <- data
  groups <- unique(data[[by]])
  for (g in groups) {
    idx <- da[[by]] == g
    data_centered[idx, what] <- data[idx, what] - median(data[idx, what], na.rm = TRUE)
  }
  #
  ggp_base_centered <- ggplot2::ggplot(data_centered) +
    ggplot2::theme_classic()
  #

  if (show_colors) {
    ggcdf <- ggp_base_centered +
      ggplot2::aes(x = !!rlang::sym(what), color = !!rlang::sym(by)) +
      ggplot2::stat_ecdf(linewidth = 1) +
      ggplot2::labs(
        title = "Empirical CDFs of Median-Centered Data",
        x = ylab,
        y = "Cumulative Probability"
      ) +
      ggplot2::theme(legend.position = if(show_legend) "bottom" else "none") +
      color_scale_color
  } else {
    ggcdf <- ggp_base +
      ggplot2::aes(x = !!rlang::sym(what)) +
      ggplot2::stat_ecdf(linewidth = 1) +
      ggplot2::facet_wrap(rlang::sym(by)) +
      ggplot2::labs(
        title = "Empirical CDFs of Median-Centered Data",
        x = ylab,
        y = "Cumulative Probability"
      )
  }
  
  ##############################################################################
  #
  # add confidence band
  #
  ##############################################################################
  #
  # add confidence band using sfsmisc::KSd 
  #
  calc_simultaneous_ecdf_bands_sfsmisc <- function(x, alpha = 0.05) {
    
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
  #
  
  ##############################################################################
  #
  # add_simultaneous_bands function using sfsmisc
  #
  ##############################################################################
  #
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
  
  # call it 
  ggcdf <- add_simultaneous_bands_sfsmisc(p = ggcdf, 
                                 data = data, 
                                 response_col = what,
                                 group_col = by,
                                 ref_level = ref_level,
                                 alpha = 0.05)
  

  ##############################################################################
  #
  # add to plot output from test_location_shift_assumption()
  #
  ##############################################################################
  #
  values <- data[[what]]  #ta[[response_var]]
  groups <- data[[by]]  #ta[[group_var]]
  #
  x_vals <- values[groups != ref_level]    # cases
  y_vals <- values[groups ==  ref_level]   # controls
  #
  test_loc = test_shift_equivalence(x_vals, y_vals, test = test)
  #
  ##########################################################
  #
  # Add assumption test annotation to the CDF plot
  test_name <- switch(test,
                      "ks" = "Kolmogorov-Smirnov test",
                      "cvm" = "Cramer-von Mises test",
                      "ad" = "Anderson-Darling test",
                      "wass" = "Wasserstein test",
                      "dts" = "DTS test"
  )
  test_loc$assumption.violated <- test_loc$p.value < 0.05
  assumption_test_text <- paste0(test_name, " p-value ",
                                 ifelse(test_loc$p.value < 0.001, "< 0.001", paste0('= ', formatC(test_loc$p.value, format = "f", digits = 3))))
  if (test_loc$assumption.violated) {
    assumption_test_text <- paste0(assumption_test_text, " (assumption violated)")
  }
  #
  ggcdf <- ggcdf + ggplot2::labs(subtitle = assumption_test_text) #+
    #ggplot2::theme(legend.position = "none")



  ##############################################################################
  #
  #
  #
  ##############################################################################
  #
  # Combine plots if requested
  if (combine_plots) {
    combined <- (ggb + ggd) / (ggw_centered + ggcdf)
    return(combined)
  } else {
    return(list(ggb = ggb, ggd = ggd, ggw_centered = ggw_centered, ggcdf = ggcdf))
  }
}
