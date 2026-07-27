# Compare the plot-size methods on the same trial

Runs the CV-based methods on one set of data and tabulates what each
recommends: the optimal plot size \\X_o\\, the CV at that size, and
comparable fit statistics. Given the raw grid instead of a CV table, it
also builds the table with
[`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md)
and adds
[`calc_paranaiba()`](https://willyanjnr.github.io/trialSizing/reference/calc_paranaiba.md),
which works on the basic units rather than on CV values.

## Usage

``` r
compare_methods(
  .data,
  x = NULL,
  cv = NULL,
  trial = NULL,
  methods = NULL,
  step = 0.001,
  weights = FALSE,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95
)
```

## Arguments

- .data:

  one of: a data frame holding a CV table (with `x` and `cv` columns); a
  matrix, being the raw grid of basic experimental units of one trial;
  or a named list of such matrices. Grids are expanded into CV tables
  with
  [`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md),
  and only then can the Paranaiba method take part.

- x, cv:

  column names in `.data` when it is a data frame (default `"x"` and
  `"cv"`, the names
  [`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md)
  produces).

- trial:

  optional column name identifying the trial.

- methods:

  which methods to run; `NULL` (default) runs every applicable one. Any
  subset of `"mcm"`, `"lrp"`, `"qrp"` and `"paranaiba"`, the last
  requiring a grid.

- step:

  grid step for the breakpoint search of the LRP and QRP.

- weights:

  weighting for the fits, as in
  [`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md).
  `TRUE` uses the `n` column, which a grid always provides. It reaches
  the MCM through that method's Federer `df` argument, the same weighted
  least squares by another name. The Paranaiba row is unaffected: it is
  a closed form over the basic units and never sees the CV table.

- bootstrap:

  logical; add a confidence interval for each \\X_o\\, and the
  breakpoint-existence p-value where the method has one.

- n_boot, conf_level:

  bootstrap size and confidence level.

## Value

An object of class `"method_comparison"`: a list with `summary` (one row
per method and trial: `trial`, `method`, `Xo`, `CVxo`, `R2`, `RMSE`, and
`Xo_lwr`, `Xo_upr`, `p_breakpoint` when `bootstrap = TRUE`) and `meta`.

## What the table is for

The methods disagree systematically: the MCM optimum is the smallest,
the LRP intermediate and the QRP the largest, an ordering reported
across many crops. Seeing the three side by side, with intervals when
`bootstrap = TRUE`, shows whether that ordering is a real difference or
three readings of the same imprecise quantity.

## Which statistics are comparable

`R2` and `RMSE` are recomputed here from the residuals on the original
CV scale, so they mean the same thing for every method even when the fit
itself was weighted. AIC and BIC are deliberately absent: the MCM has
two parameters against the plateau models' three plus a breakpoint, and
the breakpoint is not an ordinary parameter, so the information criteria
are not on a common footing. The Paranaiba estimate is a closed form
with no fitted residuals, so its `R2` and `RMSE` are `NA` by nature, not
by omission.

## References

Cargnelutti Filho, A. et al. (2025). *Revista Vivencias*, 21(43),
499-513, which compares the same three methods and reports 4.81, 7.19
and 10.25 m2 for the MCM, LRP and QRP respectively.

## See also

[`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md),
[`fit_qrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_qrp.md),
[`fit_mcm()`](https://willyanjnr.github.io/trialSizing/reference/fit_mcm.md),
[`calc_paranaiba()`](https://willyanjnr.github.io/trialSizing/reference/calc_paranaiba.md),
[`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md)

## Examples

``` r
## The bundled simulated uniformity trial (see ?uniformity_trial): one
## 8 x 12 grid of 1 m2 basic units per trial.
grid1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
                                    grep("^col", names(uniformity_trial))])

## From a CV table (built here from the grid)
cv_tab <- calc_cv_shapes(list(T1 = grid1))
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

# \donttest{
## From the raw grid: the CV table is built on the way, and the Paranaiba
## method joins because it needs the basic units.
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

## Weighted by the number of plots per shape
compare_methods(grid1, step = 0.05, weights = TRUE)$summary
#> Comparing 4 method(s) on 1 trial(s).
#>     trial    method       Xo      CVxo        R2     RMSE
#> 1 Trial 1       MCM 4.506414 12.726291 0.9753314 0.725023
#> 2 Trial 1       LRP 5.050000  9.156181 0.7937801 2.096258
#> 3 Trial 1       QRP 7.850000  8.682611 0.8502828 1.786139
#> 4 Trial 1 Paranaiba 4.793933 10.719560        NA       NA

## With intervals, which is the point: the spread between methods is often
## smaller than the uncertainty within each one.
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

## Several trials at once
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
# }
```
