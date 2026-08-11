# ROC plot with confidence band – internal function

ROC plot with confidence band – internal function

## Usage

``` r
plot_roc(x, ...)
```

## Arguments

- x:

  Object of class `roc_ci` returned by
  [`roc_with_ci()`](https://grendar.github.io/wmwAUC/reference/roc_with_ci.md)

- ...:

  not used

## Value

No return value, called for side effects. Creates an ROC curve plot
showing the receiver operating characteristic with AUC information and
confidence intervals if available.
