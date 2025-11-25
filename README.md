Wilcoxon-Mann-Whitney Test of No Group Discrimination (Continuous
Variables)
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# wmwAUC

<!-- badges: start -->
<!-- badges: end -->

The goal of wmwAUC is to provide the inferences for the
Wilcoxon-Mann-Whitney test of $\mathrm{H0: AUC = 0.5}$.

## Installation

You can install the development version of wmwAUC like so:

``` r
devtools::install_github('grendar/wmwAUC')
```

## Simulation 1

Consider the setting of two zero-mean different-scale gaussians. Then
the traditional $\mathrm{H0: F = G}$ of WMW test is false and
$\mathrm{H1: F \neq G}$ holds.

The Monte Carlo simulation demonstrates that the normalized test
statistic $U/(n_1n_2)$ which is just eAUC, concentrates asymptotically
on 0.5 - the value expected under a true null hypothesis.

If WMW tested distributional equality, the test statistic should not
concentrate on its null value when distributions clearly differ.

Also note that under $\mathrm{H1: F \neq G}$, p-values should
concentrate near zero, yet the observed distribution is nearly uniform
with a slightly elevated first bins, consistent with testing a true null
hypothesis ($\mathrm{AUC = 0.5}$) using miscalibrated variance
estimation.

``` r
#############################################################################
#
# Simulation 1: H0: F=G is erroneously too broad
#
#############################################################################
#
# This simulation takes several minutes to complete
# N = 10000
# n = 1000
# set.seed(123L)
# pval_wt = pval_wmw = eauc = numeric(N)
# for (i in 1:N) {
#
#  x = rnorm(n, sd = 0.1)
#   y = rnorm(n, sd = 3)
#   # wilcox.test() of H0: F = G
#   wt = wilcox.test(x, y)
#   pval_wt[i] = wt$p.value
#   # wmw_test() of H0: AUC = 0.5
#   pval_wmw[i] = wmw_pvalue(x, y)
#   # eAUC
#   eauc[i] = wt$statistic/(n*n)
#   #
# }
data(simulation1)  # List eauc, pval_wt, pval_wmw
#
```

<img src="man/figures/README-hist_sim1_1-1.png" width="50%" style="display: block; margin: auto;" />
Empirical AUC centered at 0.5 despite F != G

<img src="man/figures/README-hist_sim1_2-1.png" width="50%" style="display: block; margin: auto;" />
Traditional p-values should under H1 concentrate near 0.

<img src="man/figures/README-hist_sim1_3-1.png" width="50%" style="display: block; margin: auto;" />
Correct p-values (correctly testing AUC = 0.5)

## Simulation 2

``` r
#############################################################################
#
# Simulation 2: H1: F stoch. dominates G is erroneously too narrow
#               WMW is consistent for broader H1: AUC != 0.5   
#
#############################################################################
#
#
# This simulation takes several minutes to complete
# N = 10000
# n = 1000
# set.seed(123L)
# pval_wt = pval_wmw = eauc = numeric(N)
# for (i in 1:N) {
# #
#   # gaussians with different location and scale
#   # does not satisfy stochastic dominance
#   x = rnorm(n, 0, sd = 0.1)
#   y = rnorm(n, 0.5, sd = 3)
#   # wilcox.test H0: F = G vs H1: (F stochastically dominates G) OR (G stochastically dominates F)
#   wt = wilcox.test(x, y)
#   pval_wt[i] = wt$p.value
#   # wmw_test H0: AUC = 0.5 vs H1: AUC neq 0.5
#   pval_wmw[i] = wmw_pvalue(x, y)
#   # eAUC
#   eauc[i] = wt$statistic/(n*n)
# #
# }
data(simulation2)  # List of eauc, pval_wt, pval_wmw
# WMW detects broader alternatives than traditional stochastic dominance
```

<img src="man/figures/README-hist_sim2-1.png" width="50%" style="display: block; margin: auto;" /><img src="man/figures/README-hist_sim2-2.png" width="50%" style="display: block; margin: auto;" /><img src="man/figures/README-hist_sim2-3.png" width="50%" style="display: block; margin: auto;" />

## Example 1

Real data analyzed by WMW test of no group discrimination

``` r
data(gemR::MS)
da <- MS

# preparing data frame
class(da$proteins) <- setdiff(class(da$proteins), "AsIs")
df <- as.data.frame(da$proteins)
df$MS <- da$MS
```

### Test of no group discrimination

``` r
wmd <- wmw_test(P19099 ~ MS, data = df, ref_level = 'no')
wmd
#> 
#>         Wilcoxon-Mann-Whitney Test of No Group Discrimination
#> 
#> data: P19099 by MS (n1 = 37, n2 = 64)
#> groups: yes vs no (reference)
#> W = 1726, p-value = 0.001131
#> alternative hypothesis for AUC: two.sided 
#> 95 percent confidence interval for AUC (hanley): 
#>  0.623 0.835
#> empirical AUC (eAUC):
#>  0.729
```

<img src="man/figures/README-plot_ex1-1.png" width="50%" style="display: block; margin: auto;" />

## Example 2

Synthetic data illustrating the special case of location shift
assumption.

``` r
data(Ex2)
da <- Ex2

# WMW test
wmd <- wmw_test(y ~ group, data = da, ref_level = 'control')
wmd
#> 
#>         Wilcoxon-Mann-Whitney Test of No Group Discrimination
#> 
#> data: y by group (n1 = 100, n2 = 100)
#> groups: case vs control (reference)
#> W = 3705, p-value = 0.003197
#> alternative hypothesis for AUC: two.sided 
#> 95 percent confidence interval for AUC (hanley): 
#>  0.294 0.447
#> empirical AUC (eAUC):
#>  0.370
```

<img src="man/figures/README-ROC_example2-1.png" width="50%" style="display: block; margin: auto;" />

### Check location-shift assumption with EDA

<img src="man/figures/README-quadruplot_Ex2-1.png" width="100%" style="display: block; margin: auto;" />
location-shift assumption not tenable

Note that medians are essentially the same:

``` r
median(da$y[da$group == 'case'])
#> [1] 0.4949383
median(da$y[da$group == 'control'])
#> [1] 0.4926145
```

Erroneous use of location-shift special case would falsely conclude
significant median difference despite identical medians

    #> 
    #>         Wilcoxon-Mann-Whitney Test of No Group Discrimination
    #> 
    #> data: y by group (n1 = 100, n2 = 100)
    #> groups: case vs control (reference)
    #> W = 3705, p-value = 0.003197
    #> alternative hypothesis for AUC: two.sided 
    #> 95 percent confidence interval for AUC (hanley): 
    #>  0.294 0.447
    #> empirical AUC (eAUC):
    #>  0.370 [discrimination effect size]
    #> 
    #> Location-shift analysis (under F1(x) = F2(x - delta)):
    #> alternative hypothesis for location: two.sided 
    #> 95 percent confidence interval for median of all pairwise distances:
    #>  -0.094 -0.016
    #> Hodges-Lehmann median of all pairwise distances:
    #>  -0.048 [location effect size: eAUC = 0.370]

## Example 3

WMW applied to another real-life data set.

``` r
data(wesdr)
da = wesdr
da$ret = as.factor(da$ret)
# WMW 
wmd <- wmw_test(bmi ~ ret, data = da, ref_level = '0')
wmd
#> 
#>         Wilcoxon-Mann-Whitney Test of No Group Discrimination
#> 
#> data: bmi by ret (n1 = 278, n2 = 391)
#> groups: 1 vs 0 (reference)
#> W = 59417.5, p-value = 0.037970
#> alternative hypothesis for AUC: two.sided 
#> 95 percent confidence interval for AUC (hanley): 
#>  0.502 0.591
#> empirical AUC (eAUC):
#>  0.547
```

<img src="man/figures/README-roc_Ex3-1.png" width="50%" style="display: block; margin: auto;" />

### EDA to assess location shift assumption validity

<img src="man/figures/README-quadruplot_Ex3-1.png" width="100%" style="display: block; margin: auto;" />
hence, location shift assumption is tenable

### Special case of WMW test

``` r
suppressWarnings({ # ties in data
wml <- wmw_test(bmi ~ ret, data = da, ref_level = '0', 
                 ci_method = 'boot', special_case = TRUE)
})                 
wml
#> 
#>         Wilcoxon-Mann-Whitney Test of No Group Discrimination
#> 
#> data: bmi by ret (n1 = 278, n2 = 391)
#> groups: 1 vs 0 (reference)
#> W = 59417.5, p-value = 0.037970
#> alternative hypothesis for AUC: two.sided 
#> 95 percent confidence interval for AUC (boot): 
#>  0.499 0.591
#> empirical AUC (eAUC):
#>  0.547 [discrimination effect size]
#> 
#> Location-shift analysis (under F1(x) = F2(x - delta)):
#> alternative hypothesis for location: two.sided 
#> 95 percent confidence interval for median of all pairwise distances:
#>  0.000 1.100
#> Hodges-Lehmann median of all pairwise distances:
#>  0.600 [location effect size: eAUC = 0.547]
```

Plot
<img src="man/figures/README-plot_ex3-1.png" width="100%" style="display: block; margin: auto;" />

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
