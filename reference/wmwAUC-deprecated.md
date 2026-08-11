# Deprecated Function Names (wmwAUC 1.0.0)

These functions have been renamed in wmwAUC 1.0.0, alongside argument
and behavior changes in some cases. Calling them raises an informative
error pointing to the current function and to NEWS.md, rather than
attempting to forward the call – for functions with argument-name or
argument-value changes (see Details), a silent passthrough could
otherwise produce a different result than the old function used to,
without any indication that anything had changed.

## Usage

``` r
wmw_test(...)

wmw_pvalue(...)

wmw_pvalue_ties(...)

pseudomedian_ci(...)
```

## Arguments

- ...:

  Not used; calling any of these functions always raises an error.

## Details

- `wmw_test()` -\>
  [`wmwAUC_test`](https://grendar.github.io/wmwAUC/reference/wmwAUC_test.md).
  The `special_case` argument was renamed `pseudomedian`, and the
  `ci_method = 'hanley'` option was removed (not renamed), and replaced
  by the DeLong et al. (1988) method.

- `wmw_pvalue()` -\>
  [`wmwAUC_pvalue_BC`](https://grendar.github.io/wmwAUC/reference/wmwAUC_pvalue_BC.md).

- `wmw_pvalue_ties()` -\>
  [`wmwAUC_pvalue_EU`](https://grendar.github.io/wmwAUC/reference/wmwAUC_pvalue_EU.md).

- `pseudomedian_ci()` -\>
  [`wmwAUC_pseudomedian_ci`](https://grendar.github.io/wmwAUC/reference/wmwAUC_pseudomedian_ci.md).
