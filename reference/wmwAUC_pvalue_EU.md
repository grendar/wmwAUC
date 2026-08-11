# P-value for the Exact Unbiased (EU) Test of AUC

Tests \\H_0\colon \mathrm{AUC} = A_0\\ vs the specified alternative,
using an exact finite-sample unbiased variance estimator with the
mid-rank kernel.

## Usage

``` r
wmwAUC_pvalue_EU(
  x,
  y,
  alternative = "two.sided",
  A0 = 0.5,
  max_exact = 10000,
  n_perm = 2000
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
  Only supported when `length(x) + length(y) >= 20` (see Details); an
  error is raised otherwise.

- max_exact:

  Integer; the permutation branch (used when
  `length(x) + length(y) < 20`) enumerates all permutations exactly when
  their count is at most `max_exact`, and falls back to Monte Carlo
  sampling above that. Default 10000.

- n_perm:

  Integer; number of Monte Carlo permutation replicates used when exact
  enumeration is not feasible. Default 2000.

## Value

Numeric p-value.

## Details

Uses two-tier approach: studentized permutation for
`length(x) + length(y) < 20` and the exact finite-sample unbiased
estimator for `length(x) + length(y) >= 20`.

For `length(x) + length(y) < 20`, a studentized permutation test is
used: the same t-statistic is recomputed on each permuted split of the
pooled data, and the p-value is the proportion of permuted statistics at
least as extreme as the observed one. This permutation scheme relies on
group-relabeling exchangeability, which preserves \\H_0\\ only when
\\A_0 = 0.5\\; general `A0` is therefore not currently supported in this
small-sample regime and will raise an error.

For `length(x) + length(y) >= 20`, EU estimates \\\mathrm{Var}(\hat A)\\
by the exact finite-sample unbiased combination derived from the
Hoeffding decomposition of the mid-rank kernel, with Welch–Satterthwaite
degrees of freedom.

`x` is taken to represent cases and `y` the reference/control group,
matching the convention of
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html). Internally,
the test statistic and variance components are computed in the
\\P(X\<Y)\\ framework.
