## Summary

This is a major update (0.2.0 -> 1.0.0) correcting several statistical
errors identified in the previous version, with accompanying breaking API
changes. See NEWS.md for full details. Highlights:

- Corrected variance formula and degrees of freedom for the Exact Unbiased (EU) method.
- Corrected variance formula of the Bias Correction (BC) method. 
- Fixed a directionality bug in wmwAUC_test()'s one-sided alternative
  handling.
- Both methods now handle tied data consistently via a mid-rank kernel.
- Replaced the Hanley-McNeil AUC confidence interval with the DeLong et al. (1988) method.
- Renamed core functions (wmw_test -> wmwAUC_test, wmw_pvalue ->
  wmwAUC_pvalue_BC, wmw_pvalue_ties -> wmwAUC_pvalue_EU, pseudomedian_ci ->
  wmwAUC_pseudomedian_ci). Old names now raise an informative error rather
  than being silently forwarded, since several involve argument-name or
  argument-value changes that could otherwise produce silently different
  results.
- Removed quadruplot(), test_shift_equivalence(), and two sfsmisc-based
  helper functions with no remaining callers.

## Test environments

* CentOS / R 4.4.0

## R CMD check results

0 errors | 0 warnings | 1 notes (unable to verify current time)


## Downstream dependencies

None. 
