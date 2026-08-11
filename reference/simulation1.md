# Synthetic data

Synthetic data

## Usage

``` r
data(simulation1)
```

## Format

A list containing simulation results (N=10000, n=1000):

- eauc:

  Empirical AUC values

- pval_wt:

  Traditional wilcox.test p-values

- pval_wmwAUC:

  wmwAUC p-values under H0: AUC = 0.5
