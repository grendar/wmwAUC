# wmwAUC (development version)

* Reorganized README structure.
* Added reference to Grendár (2026), arXiv:2608.11177, in README.

# wmwAUC 1.0.0

Major revision correcting several errors identified in 0.2.0. Function names
changed (see Renaming below).

## Bug fixes

### In `wmw_pvalue_ties()` (renamed to `wmwAUC_pvalue_EU`)

`Var_hat_unbiased` (the exact finite-sample unbiased estimator from Theorem 1) 
is used directly, with no further rescaling. In other words,
correction_factor <- (1 - 1/n1 - 1/n2) is a bug; removed.

Errors in exact expectations of `hat zeta1^2` and `hat zeta2^2`,
which propagate to estimator of `Var(Ahat)` are corrected in the v3 of 
the arXiv preprint (arXiv:2511.20308).
The corrected formulas are implemented in `wmwAUC_pvalue_EU()`.

The three-tier approach is replaced by two-tier, with the exact unbiased estimates
used for `{length(x) + length(y) >= 20}`, and studentized permutation 
(instead of the  naive permutation) test  for `{length(x) + length(y) < 20}`.

Uses natural degrees of freedom (n1-1, n2-1, M-1) for the three 
Satterthwaite components.


### In `wmw_pvalue()` (renamed to `wmwAUC_pvalue_BC`) 

The three-tier approach is replaced by one-tier, with the bias-adjusted sigma-squared
estimator:

#### BC variance estimator: denominator fix and removal of the higher-order correction
 
Two related corrections to BC's placement-variance components:
 
- **Denominator fix.** `var_G_X`/`var_F_Y` (and the corresponding manuscript
  quantities `zeta1hat^2`/`zeta2hat^2`) now use `n1`/`n2` as the denominator,
  not `n1-1`/`n2-1`. Bessel's `-1` correction is only justified when centering
  at an *estimated* mean; here each term centers at the fixed, known value
  `1-A0`/`A0`, so it does not apply. With the corrected denominator,
  `E[zeta1hat^2] = zeta1^2 + omega1/n2` holds *exactly* (no asymptotic
  remainder).
  
- **Higher-order correction removed.** The `(1 - 1/n1 - 1/n2)` shrinkage
  applied on top of `sigma2_adj` is no longer applied. It was compensating
  for the bias the denominator fix now removes at the source; once that bias
  is gone, applying the old shrinkage on top overcorrects, underestimating
  the true variance by 11-23% (worst at small n) rather than correcting it.
  `sigma2_adj` alone is now the final variance estimator BC uses in its test
  statistic; no `sigma2_final` remains.



### In `wmw_test()` (renamed to `wmwAUC_test`) 

Bug in treating one-sided alternative, fixed.

`wmwAUC_test()`'s displayed `auc` and its `alternative` argument use the
`P(case > reference)` convention (matching `wilcox.test()`), but
`wmwAUC_pvalue_EU()`/`wmwAUC_pvalue_BC()` internally test in the `P(X < Y)`
convention. Calling them directly as `wmwAUC_pvalue_EU(x_vals, y_vals, ...)`
would silently invert the meaning of `alternative = "greater"`/`"less"`
relative to what `wmwAUC_test()` documents (two-sided p-values are
unaffected, since the test is symmetric there). Fixed by swapping the
argument order at the call site --
`wmwAUC_pvalue_EU(y_vals, x_vals, alternative = alternative)` -- rather than
altering the `alternative` string itself, so that the internal `P(X < Y)`
becomes `P(reference < case) = P(case > reference)`, matching
`wmwAUC_test()`'s own convention exactly. Verified by direct construction of
a case where cases are stochastically larger than reference: with the fix,
`alternative = "greater"` correctly gives a small p-value and `"less"` a
large one; without it, both come out backwards. `wmwAUC_pseudomedian_ci()`
is not affected -- it has no `alternative` argument and always tests
two-sided internally.


### Degrees-of-freedom caveat (both methods)

Both methods' Welch-Satterthwaite degrees of freedom are heuristics, not
formulas with a proven accuracy order:
  
- **EU**: the standardized pivot showed real undercoverage at small or
imbalanced n (up to ~4-6 Monte Carlo standard errors in some
configurations), traced to correlation between the estimated degrees of
freedom and the variance estimate itself, not to a shape mismatch in the
Satterthwaite approximation.

- **BC**: the standardized pivot is over-concentrated relative to its
claimed t-distribution even under homoscedasticity (e.g. variance of the
probability-integral-transformed pivot ~0.078 against the nominal 0.083),
and this miscalibration worsens monotonically with the heteroscedasticity
ratio.


## New features 

### In `wmwAUC_pvalue_BC()`

The function `wmwAUC_pvalue_BC()` allows for ties in data.
Though Sect 4 of the arXiv preprint assumes continuous data,
simulation confirmed that BC’s conservative behavior persists under ties 
(size remained below nominal across tie proportions from 10% to 99.8%, 
at n1 = n2 = 15, homoskedastic), though this extension is validated 
empirically rather than re-derived theoretically.
`wmwAUC 1.0.0` implements in `wmwAUC_pvalue_BC()` mid-rank placement values 
consistently in every place.

Each of the two bias-corrected variance components (`var_G_X`, `var_F_Y`) 
is now floored independently at a small epsilon if it would go negative. 
A prior both-or-nothing version reverted *both* components to their naive 
(biased) values if *either* one alone would go negative.


## In `wmwAUC_test()`

The `ci_method = hanley` is replaced by  DeLong et al. (1988) confidence interval formula; `ci_method = delong`.

## In `plot.wmwAUC_test()`

By default produces ggplot comprising boxplot overlaid with swarmplot of data,
and the empirical ROC curve plot. Allows for adding density plot.


## General A0 support

`A0` is now a parameter (default 0.5) rather than hardcoded, in both
`wmwAUC_pvalue_BC()` and the asymptotic branch of `wmwAUC_pvalue_EU()`
(`n1+n2 >= 20`). BC's higher-order correction was validated by simulation
across A0 in {0.14, 0.50, 0.86}, both balanced and imbalanced n, with no
evidence of A0-specific miscalibration (the ratio of estimated to true
variance stayed consistently in 0.90-0.93 throughout). EU's asymptotic
branch was not restricted to A0=0.5 in its derivation to begin with.

General `A0` is **not** supported in EU's permutation branch
(`n1+n2 < 20`): group-relabeling exchangeability only preserves H0 when
A0=0.5 (relabeling maps A0 to 1-A0 otherwise), and `wmwAUC_pvalue_EU()` will
error rather than silently return an unvalidated result if this is
attempted. A general-A0 small-sample method (e.g. a bootstrap resampled
under H0 rather than relying on relabeling symmetry) is deferred to a later
version.

Full general-A0 support (including the small-sample branch, and exposure
through `wmwAUC_test()`) is deferred to a later release; 1.0.0's scope is
limited to correcting the errors above.



## Abandoned

Unlike the previous version, `wmwAUC 1.0.0` does not promote using the `wmwAUC`
test for testing equality of medians. Consequently, functions
`calc_simultaneous_ecdf_bands_sfsmisc()`, `add_simultaneous_bands_sfsmisc()`, `test_shift_equivalence()` and `quadruplot()`, are abandoned.


## Renaming

- `wmw_test()` -> `wmwAUC_test()`
- `wmw_pvalue()` -> `wmwAUC_pvalue_BC()`
- `wmw_pvalue_ties()` -> `wmwAUC_pvalue_EU()`
- `pseudomedian_ci()` -> `wmwAUC_pseudomedian_ci()`

The old names are not retained as functional aliases. Calling any of them
raises an informative error naming the new function, rather than forwarding
the call -- `wmw_test()` in particular has argument changes (`special_case`
renamed to `pseudomedian`; `ci_method = 'hanley'` removed, not renamed) that
a silent passthrough could get wrong without any indication that anything
had changed.

