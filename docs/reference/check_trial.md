# Check a uniformity trial before sizing plots

Inspects the raw grid of basic experimental units and reports what
should be settled before any plot-size method is run: whether the grid
is complete and usable, whether the values are sane, and whether the
field has spatial structure worth sizing plots against. The result
carries a [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method drawing the field map.

## Usage

``` r
check_trial(
  .data,
  value = NULL,
  row_id = NULL,
  col_id = NULL,
  trial = NULL,
  cell_size = c(1, 1),
  n_bins = 15,
  max_dist = NULL,
  variogram = TRUE
)
```

## Arguments

- .data:

  a matrix (one trial), a named list of matrices, or a long data frame
  with one row per basic unit.

- value:

  name of the measurement column (long format).

- row_id, col_id:

  names of the row and column index columns (long format).

- trial:

  optional column name identifying the trial.

- cell_size:

  numeric `c(row, col)`: the distance between the centres of adjacent
  basic units, in the units you want distances reported in. Default
  `c(1, 1)`.

- n_bins:

  number of distance classes for the empirical variogram.

- max_dist:

  largest separation used; `NULL` (default) takes half the maximum
  distance in the field, the usual rule.

- variogram:

  logical; fit the variogram (default `TRUE`). Turning it off skips the
  only appreciable computation.

## Value

An object of class `"trial_check"`: a list with `checks` (one element
per trial, each holding the dimensions, statistics, outliers, trend,
spatial measures, fitted `variogram` and any `issues`) and `summary`,
one row per trial.

## What is checked

- Structure:

  Grid dimensions, missing cells, and how many rectangular plot shapes
  the grid admits. A grid whose sides are prime yields almost no shapes,
  which stops the CV-based methods before they start; the count is
  flagged when it falls below six.

- Values:

  Missing, negative, zero and constant values, and outliers by the
  boxplot rule. In a uniformity trial an outlying basic unit is usually
  a harvest failure or a typing error, and it contaminates every shape
  that contains it.

- Trend:

  Means along rows and columns, with the p-value of a monotone trend. A
  gradient in one direction is field fertility, and it is why shapes of
  equal area but different orientation give different CVs.

- Spatial structure:

  Moran's I with its p-value, the first-order autocorrelations used by
  [`calc_paranaiba()`](https://willyanjnr.github.io/trialSizing/reference/calc_paranaiba.md),
  and a fitted variogram.

## Why the variogram matters here

The fitted model gives three numbers that bear directly on plot size.
The *range* is the distance beyond which basic units stop being
correlated: a plot larger than the range is averaging units that are
already independent, which is the spatial reading of the plateau that
[`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md)
and
[`fit_qrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_qrp.md)
estimate empirically. The *nugget-to-sill ratio* is the share of
variance with no spatial structure; following Cambardella et al. (1994)
it is read as strong (below 0.25), moderate (0.25 to 0.75) or weak
(above 0.75) spatial dependence. When it is weak, the field varies
almost at random from unit to unit and no choice of plot size will buy
much precision – worth knowing before fitting five models to the CV
table.

The model is fitted by profiling the range over a grid and solving the
nugget and partial sill exactly at each candidate, under non-negativity.
Spherical, exponential and Gaussian shapes are all tried and the best
weighted fit is kept; exponential and Gaussian use the effective-range
convention, where the range is the distance reaching 95% of the sill.

## References

Matheron, G. (1963). Principles of geostatistics. *Economic Geology*,
58(8), 1246-1266.  
Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
*Biometrika*, 37(1/2), 17-23.  
Cambardella, C. A. et al. (1994). Field-scale variability of soil
properties in central Iowa soils. *Soil Science Society of America
Journal*, 58(5), 1501-1511.  
Cressie, N. (1993). *Statistics for Spatial Data*. Wiley, New York.  
Webster, R. & Oliver, M. A. (2007). *Geostatistics for Environmental
Scientists*, 2nd ed. Wiley, Chichester.

## See also

[`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md),
[`calc_paranaiba()`](https://willyanjnr.github.io/trialSizing/reference/calc_paranaiba.md),
[`plot.trial_check()`](https://willyanjnr.github.io/trialSizing/reference/plot.trial_check.md)

## Examples

``` r
## The bundled simulated uniformity trial (see ?uniformity_trial)
grid1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
                                    grep("^col", names(uniformity_trial))])
chk <- check_trial(grid1)
#> Checking 1 trial(s).
chk
#> Uniformity trial check -- Trial 1 
#>   Grid:      8 x 12 = 96 basic units, 0 missing
#>   Shapes:    23 rectangular plot shapes available
#>   Values:    mean 251.013, sd 58.922, CV 23.47%, range [114.590, 368.440]
#>   Outliers:  0 by the boxplot rule
#>   Trend:     rows p = 0.387, columns p = 0.0336
#>   Moran's I: 0.057 (p = 0.367)   rho row 0.059, col 0.104
#>   Variogram: gaussian | nugget 3223.678, sill 3715.746, range 9.54
#>              nugget/sill 0.87 -> weak spatial dependence
#>              little spatial structure: plot size will buy little
#>              precision here, whatever method is used
#>   Issues:    none
#> 

## the fitted variogram, as numbers rather than as a picture
chk$checks[[1]]$variogram[c("model", "nugget", "sill", "range",
                            "nugget_ratio", "dependence")]
#> $model
#> [1] "gaussian"
#> 
#> $nugget
#> [1] 3223.678
#> 
#> $sill
#> [1] 3715.746
#> 
#> $range
#> [1] 9.543577
#> 
#> $nugget_ratio
#> [1] 0.8675723
#> 
#> $dependence
#> [1] "weak"
#> 

# \donttest{
## the field map: kriged surface with the basic units drawn on top
plot(chk)


## several trials at once
grids <- lapply(split(uniformity_trial, uniformity_trial$trial),
                function(d) as.matrix(d[, grep("^col", names(d))]))
check_trial(grids)$summary
#> Checking 3 trial(s).
#>   trial rows cols n_missing n_shapes     mean       cv n_outliers   morans_i
#> 1    T1    8   12         0       23 251.0129 23.47375          0 0.05689653
#> 2    T2    8   12         0       23 256.1666 22.58366          4 0.10942413
#> 3    T3    8   12         0       23 244.5397 24.61861          1 0.08041496
#>      range nugget_ratio dependence n_issues
#> 1 9.543577   0.86757225       weak        0
#> 2 1.429325   0.02860592     strong        0
#> 3 1.386393   0.57994864   moderate        0

## a grid whose sides are prime admits almost no plot shapes
check_trial(matrix(rnorm(77, 100, 10), nrow = 7))$checks[[1]]$issues
#> Checking 1 trial(s).
#> [1] "only 3 plot shape(s) available from a 7 x 11 grid"
# }
```
