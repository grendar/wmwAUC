# P-value for the Bias-Corrected (BC) Test of AUC

Tests \\H_0\colon \mathrm{AUC} = A_0\\ vs the specified alternative,
using a bias-corrected finite-sample variance estimator with the
mid-rank kernel.

## Usage

``` r
wmwAUC_pvalue_BC(
  x,
  y,
  alternative = "two.sided",
  A0 = 0.5,
  min_n_warn_threshold = 10
)
```

## Arguments

- x:

  Numeric vector of cases (group 1) values.

- y:

  Numeric vector of reference/control (group 2) values.

- alternative:

  Character: `"two.sided"`, `"greater"`, or `"less"`.

- A0:

  Numeric null value of \\\mathrm{AUC} = P(X \< Y)\\. Defaults to 0.5.

- min_n_warn_threshold:

  Integer; if `min(length(x), length(y))` is below this threshold, a
  warning is issued that power may be very low at this sample size.
  Default 10.

## Value

Numeric p-value.

## Details

BC estimates \\\mathrm{Var}(\hat A)\\ by correcting each
placement-variance component for its \\O(1/n)\\ upward bias, using a
plug-in estimate of the bias subtracted from the naive placement
variance; each corrected component is floored independently at a small
\\\epsilon \> 0\\ if it would otherwise go negative. The mid-rank kernel
\\h(x,y) = 1\\x\<y\\ + \frac{1}{2} 1\\x=y\\\\ is used throughout, for
both the point estimate and the variance components.

Uses one-tier approach with \\\hat\sigma^2\_{\mathrm{adj}}\\.

BC is a conservative test: observed size stays below nominal across a
wide range of sample sizes, heteroskedasticity, and tie proportions, at
a real cost in power for small or imbalanced samples. See
`min_n_warn_threshold` and its warning text; the EU method
([`wmwAUC_pvalue_EU`](https://grendar.github.io/wmwAUC/reference/wmwAUC_pvalue_EU.md))
is recommended when `min(n1, n2)` is small.

`x` is taken to represent cases and `y` the reference/control group,
matching the convention of
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html). Internally,
the test statistic and variance components are computed in the
\\P(X\<Y)\\ framework.
