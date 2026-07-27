# Fit the Linear Response Plateau (LRP) model by grid search

Fits the linear-plateau (broken-line) model \$\$f(x) = a + b\\x \\
\text{ if } x \le X_0, \qquad f(x) = a + b\\X_0 \\ \text{ if } x \>
X_0\$\$ by profiling the breakpoint \\X_0\\ over a fine grid. For each
candidate breakpoint the linear coefficients are obtained by least
squares and the plateau is set to \\a + b\\X_0\\; the breakpoint
minimizing the residual sum of squares over all observations is
returned. Profiling the breakpoint avoids the starting-value sensitivity
and local-minima of a direct [`nls`](https://rdrr.io/r/stats/nls.html)
fit, so no initial values are required.

## Usage

``` r
fit_lrp(
  .data = NULL,
  x = NULL,
  cv = NULL,
  trial = NULL,
  step = 0.001,
  method = c("segment", "ramp"),
  search_range = NULL,
  start = NULL,
  local_min_tol = 0.1,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95,
  weights = FALSE
)
```

## Arguments

- .data:

  optional data frame. When supplied, `x`, `cv` and `trial` are
  interpreted as column names (character strings). When `NULL` (default)
  the vector interface is used.

- x, cv:

  either numeric vectors (vector interface) or, when `.data` is a data
  frame, the names of the predictor and response columns.

- trial:

  optional name of the column identifying the trial. When given, one
  model is fit per trial.

- step:

  grid step for the breakpoint search. Default `0.001`; larger values
  run faster with a slightly coarser breakpoint.

- method:

  how the linear coefficients are estimated at each candidate
  breakpoint. `"segment"` (default) uses only observations with
  `x <= X0`, reproducing the Paranaiba et al. (2009) procedure. `"ramp"`
  uses all observations on the basis `pmin(x, X0)`, the standard
  least-squares LRP.

- search_range:

  optional numeric `c(lower, upper)` restricting the interval (in units
  of `x`) where the breakpoint is searched. Must fall within the data
  range. Useful when the optimum is known to lie in a region and an
  outlier could otherwise pull the breakpoint outside it. `NULL`
  (default) searches the full feasible interval.

- start:

  optional single breakpoint value. The fit is unchanged, but the result
  also carries a `$compat` element holding the local minimum of the
  basin containing `start`: the solution a gradient-based fitter (`nls`,
  `nlsLM`) seeded there would return. Use it to reproduce published
  results obtained with such implementations, and to see how much worse
  they fit.

- local_min_tol:

  relative SSE tolerance (default 0.10) deciding which local minima
  count as competing. A second basin fitting within 10% of the optimum
  means the breakpoint is not sharply identified. Only labelling is
  affected (the `competing` column of `$local_minima` and the stars in
  [`print()`](https://rdrr.io/r/base/print.html)); the fit itself never
  changes, and no warning is issued.

- bootstrap:

  logical; if `TRUE`, also estimate the uncertainty of the breakpoint by
  resampling (default `FALSE`). Off by default because the published
  procedure reports the point estimate alone. See the section
  "Uncertainty of the breakpoint".

- n_boot:

  number of bootstrap resamples (default 1000). Used only when
  `bootstrap = TRUE`.

- conf_level:

  confidence level of the percentile interval (default 0.95). Used only
  when `bootstrap = TRUE`.

- weights:

  weights for a weighted least-squares fit. `FALSE` (default) fits
  unweighted, as the published procedure does. `TRUE` uses the `n`
  column of `.data`, the number of plots behind each CV, as returned by
  [`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md).
  A single column name or a numeric vector are also accepted. See the
  section "Weighting by the number of plots".

## Value

For a single series, an object of class `"lrp_fit"`: a list with
`coefficients` (a, b); `parameters` (Breakpoint, Breakpoint_Response,
R2, RMSE, AIC, BIC); `fitted`; `residuals`; `data`; `method`; `step`;
`local_min_tol`; `local_minima` (the competing basins, with their SSE
excess over the optimum and a `competing` flag, or `NULL`);
`sse_profile` (the SSE at every candidate breakpoint); `compat` when
`start` was given; and `bootstrap` when `bootstrap = TRUE`, a list with
`ci`, `se`, `p_value` (existence of the breakpoint), `statistic`,
`replicates`, `n_valid` and `conf_level`.  
With `trial`, an object of class `"lrp_multi"`: a list with `summary`
(one row per trial), `fits` (the individual `"lrp_fit"` objects) and
`method`.

## Details

The response is typically the coefficient of variation (CV, percent) and
the predictor the plot size. Supply the individual CV values (one per
basic-unit form), not the means per plot size; repeated `x` values are
expected.

## Competing breakpoints

With few distinct plot sizes the residual sum of squares is a stepped
function of the breakpoint, and it often has several local minima. The
grid search always returns the global optimum, but a second basin may
fit almost as well, in which case the breakpoint is not sharply
identified and different implementations legitimately disagree. Every
local minimum of the profile is reported in `$local_minima`; those
fitting within `local_min_tol` of the optimum are flagged in the
`competing` column and starred by
[`print()`](https://rdrr.io/r/base/print.html). No warning is issued: on
stepped profiles competing basins are common, so a warning would fire on
almost every fit. Inspect `$local_minima` and `$sse_profile` instead,
and lower `local_min_tol` to flag only near-ties. Gradient-based fitters
return whichever basin their starting value lands in; `start` reproduces
that.

## Weighting by the number of plots

The CV values are not equally reliable. A shape of area \\X\\ leaves \\n
= LC/X\\ plots in the grid, so the CV of the largest plot size may rest
on two plots while the smallest rests on dozens. `weights = TRUE` fits
by weighted least squares with \\n\\ as the weight, which is the natural
measure of how much information stands behind each point.

This is off by default because it is not what the published procedure
does, and it moves the answer: the small plot sizes, where the CV is
highest, gain most of the weight, so the fitted line is pulled towards
them and the breakpoint typically falls. Report both if you use it.

Two cautions. Weighting corrects for unequal information, not for
dependence: every CV in the table comes from the same grid of basic
units, so the points are not independent with or without weights. And a
weighted fit's \\R^2\\, RMSE, AIC and BIC follow the
[`lm`](https://rdrr.io/r/stats/lm.html) convention of being computed on
weighted residuals, so they cannot be compared with those of the
unweighted fit.

## Uncertainty of the breakpoint

`bootstrap = TRUE` adds two things the point estimate cannot give.

The first is a percentile confidence interval: the shapes (the rows of
the CV table) are resampled with replacement, the model is refit on each
resample, and the empirical quantiles of the resulting breakpoints form
the interval. Resamples with fewer than three distinct plot sizes cannot
place a breakpoint and are discarded; `$bootstrap$n_valid` reports how
many were kept. The interval is usually wide, which is the honest
reading of a breakpoint estimated from a handful of plot sizes.

The second is a test of whether the breakpoint exists at all. Under the
null hypothesis that the CV falls linearly and never plateaus, the
breakpoint is not identified, so the usual likelihood-ratio statistic
does not have a chi-square distribution (Davies, 1987). The null
distribution is simulated instead: residuals of the straight-line fit
are resampled, the plateau model is refit to each simulated series, and
the p-value is the proportion of simulated SSE reductions that match or
exceed the observed one. A large p-value means a straight line explains
the data as well as the plateau, and the optimal plot size should not be
read off this fit.

Both are computed on the same breakpoint grid as the main fit, so `step`
controls their resolution too. Set the random seed before calling to
make the result reproducible.

## Two ways to call

- Vectors:

  `fit_lrp(x, cv)` fits a single series and returns an `"lrp_fit"`.

- Data frame:

  `fit_lrp(.data, x = "x", cv = "cv", trial = "trial")` takes a data
  frame plus column names. With `trial`, one model is fit per trial and
  an `"lrp_multi"` object (summary table plus the individual fits) is
  returned. Column names default to `"x"`, `"cv"` and `"trial"`; missing
  columns raise a clear error.

## References

Paranaiba, P. F., Ferreira, D. F. & Morais, A. R. (2009). Tamanho otimo
de parcelas experimentais: proposicao de metodos de estimacao. *Revista
Brasileira de Biometria*, 27(2), 255-268.  
Cargnelutti Filho, A. et al. (2025). *Revista Vivencias*, 21(43),
499-513.  
Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is
present only under the alternative. *Biometrika*, 74(1), 33-43.  
Efron, B. & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall, New York.

## See also

[`plot.lrp_fit()`](https://willyanjnr.github.io/trialSizing/reference/plot.lrp_fit.md),
[`predict.lrp_fit()`](https://willyanjnr.github.io/trialSizing/reference/predict.lrp_fit.md)

## Examples

``` r
## A simulated uniformity trial bundled with the package (see
## ?uniformity_trial). Group its 1 m2 basic units into plots of every shape,
## then read one CV per shape -- the input the plateau models expect. Several
## shapes share an area, so the plot sizes repeat.
grid1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
                                    grep("^col", names(uniformity_trial))])
cv_tab <- calc_cv_shapes(list(T1 = grid1))
#> CV by shape for 1 trial(s).
X   <- cv_tab$x    # plot size (m2)
CV1 <- cv_tab$cv   # CV (%) among plots of that shape

## A coarse grid runs fast and already reproduces Xo to two decimals; the
## default step = 0.001 refines the third.
fit <- fit_lrp(X, CV1, step = 0.01)
fit
#> Linear Response Plateau (LRP) fit
#> Method:                  segment  
#> Breakpoint (Xo):         9.160 
#> CV at breakpoint:        6.945 
#> R2: 0.893  RMSE: 1.513  AIC: 92.3  BIC: 96.9 
#> 
#> Local minima of the SSE profile (8):
#>   Xo =   7.400   SSE  +15.3% vs optimum
#>   Xo =  12.820   SSE  +34.5% vs optimum
#>   Xo =   5.750   SSE  +50.1% vs optimum
#>   ... see $local_minima for all
coef(fit)
#>         a         b 
#> 21.049460 -1.539805 
fit$parameters[c("Breakpoint", "Breakpoint_Response")]
#>          Breakpoint Breakpoint_Response 
#>            9.160000            6.944847 

## CV expected at plot sizes that were not evaluated
predict(fit, newx = c(2, 5, 7.5, 15))
#> [1] 17.969850 13.350435  9.500923  6.944847

# \donttest{
## Full precision is the default; it costs about ten times more time and
## only refines the third decimal, so it is shown here rather than used
## throughout:
# fit_lrp(X, CV1)$parameters["Breakpoint"]

## Title and styling belong to plot(), not to the fit
plot(fit, title = "Uniformity trial, T1")


## Weighting by the number of plots -----------------------------------------
## The CV of a large shape rests on very few plots, that of a 1 m2 shape on
## many. Weighting by n (carried in the CV table) pulls the fit towards the
## small sizes, so Xo falls.
c(unweighted = unname(fit_lrp(X, CV1, step = 0.01)$parameters["Breakpoint"]),
  weighted   = unname(fit_lrp(X, CV1, step = 0.01,
                              weights = cv_tab$n)$parameters["Breakpoint"]))
#> unweighted   weighted 
#>       9.16       5.05 

## Straight from the grid table, `weights = TRUE` finds the n column itself
fit_lrp(cv_tab, x = "x", cv = "cv", step = 0.05,
        weights = TRUE)$parameters["Breakpoint"]
#> Using x = 'x', cv = 'cv' (single series).
#> Breakpoint 
#>       5.05 

## Uncertainty of Xo -------------------------------------------------------
## Off by default, because the published procedure reports the point alone.
set.seed(1)
unc <- fit_lrp(X, CV1, step = 0.01, bootstrap = TRUE, n_boot = 200)
unc
#> Linear Response Plateau (LRP) fit
#> Method:                  segment  
#> Breakpoint (Xo):         9.160 
#>   95% CI (percentile):  [5.436, 17.055]   SE 3.025
#>   breakpoint exists:     p = 0.0050  (200 resamples)
#> CV at breakpoint:        6.945 
#> R2: 0.893  RMSE: 1.513  AIC: 92.3  BIC: 96.9 
#> 
#> Local minima of the SSE profile (8):
#>   Xo =   7.400   SSE  +15.3% vs optimum
#>   Xo =  12.820   SSE  +34.5% vs optimum
#>   Xo =   5.750   SSE  +50.1% vs optimum
#>   ... see $local_minima for all
unc$bootstrap$ci
#> [1]  5.43575 17.05475

## p_value tests whether a breakpoint exists at all: a large value means a
## straight line explains the CV just as well, and no plateau should be read.
unc$bootstrap$p_value
#> [1] 0.004975124

## Competing breakpoints ---------------------------------------------------
## Every local minimum of the SSE profile is reported; those fitting within
## local_min_tol of the optimum are flagged as competing.
fit$local_minima
#>   breakpoint       SSE SSE_excess competing
#> 1       7.40  60.68947  0.1526208     FALSE
#> 2      12.82  70.83320  0.3452716     FALSE
#> 3       5.75  79.02721  0.5008931     FALSE
#> 4      16.00  89.28975  0.6958005     FALSE
#> 5      24.00 142.08688  1.6985290     FALSE
#> 6       3.99 158.60770  2.0122943     FALSE
#> 7      32.00 176.35459  2.3493450     FALSE
#> 8       2.99 366.62097  5.9629043     FALSE

## Lower the tolerance to flag only near-ties
fit_lrp(X, CV1, step = 0.01, local_min_tol = 0.02)$local_minima
#>   breakpoint       SSE SSE_excess competing
#> 1       7.40  60.68947  0.1526208     FALSE
#> 2      12.82  70.83320  0.3452716     FALSE
#> 3       5.75  79.02721  0.5008931     FALSE
#> 4      16.00  89.28975  0.6958005     FALSE
#> 5      24.00 142.08688  1.6985290     FALSE
#> 6       3.99 158.60770  2.0122943     FALSE
#> 7      32.00 176.35459  2.3493450     FALSE
#> 8       2.99 366.62097  5.9629043     FALSE

## The whole profile, for inspection
plot(fit$sse_profile, type = "l", xlab = "Breakpoint", ylab = "SSE")
abline(v = fit$parameters["Breakpoint"], col = "forestgreen")


## What a gradient fitter (nls, nlsLM) seeded at 12 would have returned.
## The reported fit does not change; $compat shows the cost in SSE.
fit_lrp(X, CV1, step = 0.01, start = 12)$compat
#> $start
#> [1] 12
#> 
#> $breakpoint
#> [1] 12.82
#> 
#> $coefficients
#>          a          b 
#> 18.8969810 -0.9859507 
#> 
#> $plateau
#> [1] 6.257094
#> 
#> $SSE
#> [1] 70.8332
#> 
#> $SSE_excess
#> [1] 0.3452716
#> 

## Other arguments ---------------------------------------------------------
## Restrict the search when an outlier pulls the breakpoint away
fit_lrp(X, CV1, step = 0.01, search_range = c(5, 12))$parameters["Breakpoint"]
#> Breakpoint 
#>       9.16 

## "ramp" estimates the descending line from every observation instead of
## only those below the breakpoint, which can shift the optimum
fit_lrp(X, CV1, step = 0.01, method = "ramp")$parameters["Breakpoint"]
#> Breakpoint 
#>       9.16 

## Several trials at once --------------------------------------------------
trials <- rbind(
  data.frame(x = X, cv = CV1,        trial = "T1"),
  data.frame(x = X, cv = CV1 * 0.85, trial = "T2")
)
res <- fit_lrp(trials, x = "x", cv = "cv", trial = "trial", step = 0.01)
#> Using x = 'x', cv = 'cv', trial = 'trial' -> 2 trials.
res$summary
#>   trial       a       b breakpoint plateau     R2   RMSE    AIC    BIC n_local
#> 1    T1 21.0495 -1.5398       9.16  6.9448 0.8926 1.5130 92.321 96.863       8
#> 2    T2 17.8920 -1.3088       9.16  5.9031 0.8926 1.2861 84.845 89.387       8

## CVxo feeds the number of replications
calc_replicates(treatments = c(5, 10, 20),
                cv_percent = unname(fit$parameters["Breakpoint_Response"]),
                lsd_percent = c(10, 20))
#> Optimal number of replications
#> Design: CRD  CV: 6.94%  alpha: 0.05 
#> Rows: 6  | Non-converged: 0 
#> 
#>  Treatments CV_percent LSD_percent Alpha Design r_continuous r_optimal df_error
#>           5   6.944847          10  0.05    CRD         7.98         8       35
#>          10   6.944847          10  0.05    CRD        10.15        11      100
#>          20   6.944847          10  0.05    CRD        12.41        13      240
#>           5   6.944847          20  0.05    CRD         2.76         3       10
#>          10   6.944847          20  0.05    CRD         3.02         4       30
#>          20   6.944847          20  0.05    CRD         3.39         4       60
#>  q_tukey converged at_floor
#>    4.066      TRUE    FALSE
#>    4.577      TRUE    FALSE
#>    5.069      TRUE    FALSE
#>    4.654      TRUE    FALSE
#>    4.824      TRUE    FALSE
#>    5.241      TRUE    FALSE
# }
```
