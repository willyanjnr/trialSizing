# Comparing the methods

``` r

library(trialSizing)
```

## Why compare

The plot-size methods disagree by construction, and the published
articles report that disagreement as a finding rather than a nuisance.
The MCM optimum is typically the smallest, the LRP intermediate and the
QRP the largest – an ordering seen across many crops.
[`compare_methods()`](https://willyanjnr.github.io/trialSizing/reference/compare_methods.md)
runs them on one set of data and tabulates what each recommends, so the
spread is visible instead of inferred from separate calls.

The real question it answers is not *which method gives the largest
number* but *whether the methods disagree by more than their own
uncertainty*. With `bootstrap = TRUE` the answer comes with intervals,
and the spread between methods is often smaller than the imprecision
within any one of them.

## From a CV table

Given a data frame with `x` and `cv` columns, the three CV-based methods
run:

``` r

grid1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
                                    grep("^col", names(uniformity_trial))])
cv_tab <- calc_cv_shapes(grid1)
#> CV by shape for 1 trial(s).

cmp <- compare_methods(data.frame(x = cv_tab$x, cv = cv_tab$cv), step = 0.05)
#> Comparing 3 method(s) on 1 trial(s).
cmp
#> Plot-size methods compared
#> Trials: 1  | source: CV table  | weighted: FALSE 
#> 
#>    trial method     Xo   CVxo    R2  RMSE
#>  Trial 1    MCM  4.495 12.742 0.975 0.725
#>  Trial 1    LRP  9.150  6.960 0.893 1.513
#>  Trial 1    QRP 11.950  7.028 0.918 1.324
#> 
#> Xo ranges from 4.495 to 11.950 (a factor of 2.7).
```

The footer reports the range of $`X_o`$ across methods as a factor: how
many times larger the largest recommendation is than the smallest.

## From the raw grid

Given a grid of basic experimental units instead of a CV table,
[`compare_methods()`](https://willyanjnr.github.io/trialSizing/reference/compare_methods.md)
builds the CV table on the way with \[calc_cv_shapes()\], and the
Paranaíba method joins the comparison – it works on the basic units
directly and so is only available when the grid is supplied:

``` r

compare_methods(grid1, step = 0.05)
#> Comparing 4 method(s) on 1 trial(s).
#> Plot-size methods compared
#> Trials: 1  | source: raw grid  | weighted: FALSE 
#> 
#>    trial    method     Xo   CVxo    R2  RMSE
#>  Trial 1       MCM  4.495 12.742 0.975 0.725
#>  Trial 1       LRP  9.150  6.960 0.893 1.513
#>  Trial 1       QRP 11.950  7.028 0.918 1.324
#>  Trial 1 Paranaiba  4.794 10.720    NA    NA
#> 
#> Xo ranges from 4.495 to 11.950 (a factor of 2.7).
```

## Which statistics are comparable

`R2` and `RMSE` are **recomputed here** from the residuals on the
original CV scale, so they mean the same thing for every method even
when the fit itself was weighted – a weighted fit stores weighted
statistics, and comparing those across methods would compare different
quantities.

Two deliberate absences:

- **AIC and BIC** are not shown. The MCM has two parameters against the
  plateau models’ three plus a breakpoint, and the breakpoint is not an
  ordinary parameter, so the information criteria are not on a common
  footing.
- The Paranaíba row has `NA` for `R2` and `RMSE` by nature, not by
  omission: it is a closed form over the basic units with no fitted
  residuals.

## Weighting

`weights = TRUE` uses the `n` column – the number of plots each shape
yields – which a grid always provides. It reaches the MCM through that
method’s Federer `df` argument, the same weighted least squares by
another name. The Paranaíba row is unaffected, since it never sees the
CV table.

``` r

compare_methods(grid1, step = 0.05, weights = TRUE)$summary
#> Comparing 4 method(s) on 1 trial(s).
#>     trial    method       Xo      CVxo        R2     RMSE
#> 1 Trial 1       MCM 4.506414 12.726291 0.9753314 0.725023
#> 2 Trial 1       LRP 5.050000  9.156181 0.7937801 2.096258
#> 3 Trial 1       QRP 7.850000  8.682611 0.8502828 1.786139
#> 4 Trial 1 Paranaiba 4.793933 10.719560        NA       NA
```

## Intervals: the point of the exercise

Without intervals the table invites over-reading a difference that may
be noise. `bootstrap = TRUE` adds a confidence interval for each
$`X_o`$, and a breakpoint-existence p-value for the methods that have
one:

``` r

set.seed(1)
compare_methods(grid1, step = 0.05, bootstrap = TRUE, n_boot = 200)
#> Comparing 4 method(s) on 1 trial(s).
#> Plot-size methods compared
#> Trials: 1  | source: raw grid  | weighted: FALSE 
#> 
#>    trial    method     Xo   CVxo    R2  RMSE Xo_lwr Xo_upr p_breakpoint
#>  Trial 1       MCM  4.495 12.742 0.975 0.725  4.171  4.787           NA
#>  Trial 1       LRP  9.150  6.960 0.893 1.513  5.248 16.405        0.005
#>  Trial 1       QRP 11.950  7.028 0.918 1.324  9.786 22.413        0.005
#>  Trial 1 Paranaiba  4.794 10.720    NA    NA     NA     NA           NA
#> 
#> Xo ranges from 4.495 to 11.950 (a factor of 2.7).
#> The intervals do not all overlap: the methods differ beyond their
#>   own uncertainty.
```

The print method reads the intervals for you: if they all share a common
stretch, the methods do not disagree beyond their own uncertainty, and
the choice between them matters less than it looked. If they do not
overlap, the difference is real and worth a decision.

## Several trials at once

A named list of grids compares the methods within each trial:

``` r

grids <- lapply(split(uniformity_trial, uniformity_trial$trial),
                function(d) as.matrix(d[, grep("^col", names(d))]))

compare_methods(grids, step = 0.05)$summary
#> Comparing 4 method(s) on 3 trial(s).
#>    trial    method        Xo      CVxo        R2      RMSE
#> 1     T1       MCM  4.494865 12.741565 0.9753420 0.7248675
#> 2     T1       LRP  9.150000  6.960245 0.8925619 1.5130685
#> 3     T1       QRP 11.950000  7.028148 0.9177681 1.3237304
#> 4     T1 Paranaiba  4.793933 10.719560        NA        NA
#> 5     T2       MCM  4.956314 11.533737 0.8820733 1.9801430
#> 6     T2       LRP 14.400000  3.916190 0.8341009 2.3486180
#> 7     T2       QRP 22.500000  3.375754 0.8665844 2.1061713
#> 8     T2 Paranaiba  4.638860 10.372807        NA        NA
#> 9     T3       MCM  4.904998 12.210404 0.9067767 1.6560095
#> 10    T3       LRP 10.100000  6.101082 0.8019446 2.4137591
#> 11    T3       QRP 14.350000  5.966308 0.8298434 2.2373037
#> 12    T3 Paranaiba  4.936716 11.038833        NA        NA
```

## Choosing the methods

By default every applicable method runs. Restrict them with `methods`;
the Paranaíba method requires a grid and errors if asked for on a CV
table:

``` r

compare_methods(grid1, step = 0.05, methods = c("lrp", "qrp"))$summary
#> Comparing 2 method(s) on 1 trial(s).
#>     trial method    Xo     CVxo        R2     RMSE
#> 1 Trial 1    LRP  9.15 6.960245 0.8925619 1.513069
#> 2 Trial 1    QRP 11.95 7.028148 0.9177681 1.323730
```

## Where this fits

[`compare_methods()`](https://willyanjnr.github.io/trialSizing/reference/compare_methods.md)
brings together the methods documented one at a time in
[`vignette("lrp")`](https://willyanjnr.github.io/trialSizing/articles/lrp.md),
[`vignette("qrp")`](https://willyanjnr.github.io/trialSizing/articles/qrp.md),
[`vignette("mcm")`](https://willyanjnr.github.io/trialSizing/articles/mcm.md)
and
[`vignette("paranaiba")`](https://willyanjnr.github.io/trialSizing/articles/paranaiba.md).
Run
[`vignette("check_trial")`](https://willyanjnr.github.io/trialSizing/articles/check_trial.md)
first to know whether the trial is worth sizing plots against, and
[`vignette("replicates")`](https://willyanjnr.github.io/trialSizing/articles/replicates.md)
afterwards to turn the chosen $`CV_{Xo}`$ into a number of replications.

### References

Cargnelutti Filho, A. et al. (2025). Determinação do tamanho de parcela
para avaliar a massa de parte aérea de grão-de-bico. *Revista
Vivências*, 21(43), 499-513, which compares the same three methods and
reports 4.81, 7.19 and 10.25 m² for the MCM, LRP and QRP respectively.
