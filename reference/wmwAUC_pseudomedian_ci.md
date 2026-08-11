# Confidence Interval for the AUC-Equalizing Shift (Hodges-Lehmann Pseudomedian) via Test Inversion

Computes a confidence interval for the pseudomedian by inverting the
test \\\mathrm{H_0\colon AUC}(X, Y+\delta) = 0.5\\.

## Usage

``` r
wmwAUC_pseudomedian_ci(
  x,
  y,
  conf.level = 0.95,
  pvalue_method = c("EU", "BC"),
  n_grid = 1000
)
```

## Arguments

- x:

  numeric vector, first sample

- y:

  numeric vector, second sample

- conf.level:

  confidence level (default 0.95)

- pvalue_method:

  character, either 'EU' or 'BC'

- n_grid:

  number of grid points for search (default 1000)

## Value

list with conf.int, estimate and conf.level

## Details

The pseudomedian \\\delta\\ is defined generally, for any pair of
distributions \\F, G\\, as the shift solving \\P(X \< Y + \delta) =
0.5\\ — i.e. the value of \\\delta\\ that equalizes the AUC between
\\X\\ and the shifted \\Y + \delta\\. This definition does not require a
location-shift model relating \\F\\ and \\G\\. Under location-shift
(\\F(t) = G(t - \Delta)\\), this \\\delta\\ coincides with the classical
Hodges-Lehmann pseudomedian and equals the location difference
\\\Delta\\; this is a special case, not the definition. Outside
location-shift, \\\delta\\ remains well-defined as the AUC-equalizing
shift, but should not be read as "the median difference" between the two
groups.

The confidence interval is obtained by grid search over candidate
\\\delta\\ values, retaining those for which the AUC=0.5 test applied to
\\(X, Y+\delta)\\ is not rejected at level `1 - conf.level`. The search
range is set heuristically from \\\pm 3\\ times the robust MAD scale of
the pairwise differences around the point estimate; if the true interval
extends beyond this range, or if the acceptance region is not contiguous
(which can occur in small, imbalanced samples), a warning is issued. See
Warnings.

EU and BC pseudomedian CIs are computed on independent calls, and may
occasionally be reported as numerically identical at coarse n_grid even
though the underlying p-value curves genuinely differ — both methods'
crossing points can land in the same discrete grid cell.

At small `min(length(x), length(y))`, BC's near-zero power (see
[`wmwAUC_pvalue_BC`](https://grendar.github.io/wmwAUC/reference/wmwAUC_pvalue_BC.md))
can make its pseudomedian interval very wide and largely uninformative
rather than merely conservative, since the grid search rarely rejects
anywhere in the searched range. EU is recommended in that regime; see
the warning issued below.

When `pvalue_method = 'EU'` and `length(x) + length(y) < 20`, each grid
point's p-value may come from Monte Carlo permutation (see
[`wmwAUC_pvalue_EU`](https://grendar.github.io/wmwAUC/reference/wmwAUC_pvalue_EU.md)),
making this function's result stochastic at small sample sizes. For
reproducible results, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) immediately
beforehand.

## Warnings

Two diagnostic warnings are issued when the grid search may be
unreliable:

- **Non-contiguous acceptance region**: the interval is reported as
  `[min(accepted), max(accepted)]`, which assumes the set of
  non-rejected \\\delta\\ values is a single contiguous block. If it is
  not, the reported interval may be wider than the true acceptance
  region and should be inspected directly.

- **Search range too narrow**: if the p-value has not dropped below
  \\\alpha\\ at either edge of the search grid, the true endpoint may
  lie outside the searched range; consider widening the range or
  increasing `n_grid`.
