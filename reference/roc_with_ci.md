# ROC related computations – internal function

ROC related computations – internal function

## Usage

``` r
roc_with_ci(
  probs,
  labels,
  positive,
  auc,
  ci_method = c("none", "delong", "bootstrap"),
  n_boot = 1000,
  alpha = 0.05
)
```

## Arguments

- probs:

  Vector of class probabilities or values of continuous predictor

- labels:

  Vector, factor with two levels

- positive:

  Character giving the level that corresponds to 'case'

- auc:

  Numeric value of AUC

- ci_method:

  Character from c("none", "delong", "bootstrap")

- n_boot:

  Numeric value giving the number of bootstrap replicates (default:
  1000)

- alpha:

  Level of significance (default: 0.05)

## Value

List with components:

- roc_df:

  data frame for plotting ROC curve

- roc_band:

  data frame for plotting confidence band of ROC

- auc:

  auc

- auc_ci:

  confint for auc

## Details

`ci_method = 'bootstrap'` makes this function's result stochastic
(resampled AUC and ROC confidence band). For reproducible results, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) immediately
beforehand.
