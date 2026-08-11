# Plot Method for wmwAUC_test Objects

Creates a combined plot ("diplot"): a boxplot with beeswarm overlay of
the raw group values, alongside the empirical ROC curve (eROC) with
confidence band (when available) and AUC confidence interval. Optionally
adds a third panel with overlaid group density curves.

## Usage

``` r
# S3 method for class 'wmwAUC_test'
plot(x, show_density = FALSE, combine_plots = TRUE, ...)
```

## Arguments

- x:

  Object of class 'wmwAUC_test' returned by
  [`wmwAUC_test`](https://grendar.github.io/wmwAUC/reference/wmwAUC_test.md).

- show_density:

  Logical; if TRUE, adds a third panel with overlaid kernel density
  estimates for the two groups. Default FALSE.

- combine_plots:

  Logical; if TRUE (default), returns a single combined plot via
  patchwork; if FALSE, returns a named list of the individual ggplot2
  objects instead.

- ...:

  Additional arguments (not currently used).

## Value

If `combine_plots = TRUE`, a combined patchwork object. If
`combine_plots = FALSE`, a named list with elements `box_plot`,
`roc_plot`, and (if `show_density = TRUE`) `density_plot`.
