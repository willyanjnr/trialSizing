# CV by plot shape

``` r

library(trialSizing)
```

## The step before the models

[`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md),
[`fit_qrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_qrp.md)
and
[`fit_mcm()`](https://willyanjnr.github.io/trialSizing/reference/fit_mcm.md)
all start from the same thing: a table of coefficients of variation, one
row per planned plot size.
[`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md)
builds that table from the raw grid of basic experimental units (BEU) as
the trial was harvested, so the whole pipeline is grid in, CV table out,
model fit.

## Theory

A uniformity trial is a grid of $`L`$ rows by $`C`$ columns of BEU. A
rectangular **plot shape** groups adjacent units into blocks of $`X_L`$
rows by $`X_C`$ columns, and this is only possible when $`X_L`$ divides
$`L`$ and $`X_C`$ divides $`C`$. Each shape tiles the grid into

``` math
n = \frac{L \, C}{X_L \, X_C}
```

non-overlapping plots. For every shape the plot **totals** are formed,
and their mean, standard deviation and CV are computed:

``` math
CV = 100 \times \frac{s}{\bar{y}}
```

As plots get larger they average over more units, the totals vary
proportionally less, and the CV falls – the decreasing curve every
plot-size method is built to read.

### Shapes are not interchangeable at equal area

A $`1 \times 2`$ shape and a $`2 \times 1`$ shape cover the same two
basic units, but along different directions of the field. Whenever
fertility is not the same in both directions their CVs differ, so both
are reported and the CV table has **repeated plot sizes**
($`x = X_L X_C`$) that are genuinely different rows. This is the
table-level counterpart of the directional autocorrelation that
\[calc_paranaiba()\] uses.

### Fewer plots, less information

A shape of area $`X`$ yields $`n = LC/X`$ plots, so a large plot size
rests on few plots and its CV is estimated with little information: the
last rows of the table are the least reliable. The `n` column is
returned for exactly this reason and is the natural weight for a
weighted fit. Shapes leaving fewer than `min_plots` plots are dropped,
which by default removes only the whole grid (a single plot, which has
no variance).

## Data and basic use

``` r

grid1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
                                    grep("^col", names(uniformity_trial))])
dim(grid1)
#> [1]  8 12
```

``` r

tab <- calc_cv_shapes(grid1)
#> CV by shape for 1 trial(s).
tab
#>      trial X_L X_C  x  n       mean        sd        cv
#> 1  Trial 1   1   1  1 96   251.0129  58.92214 23.473749
#> 2  Trial 1   1   2  2 48   502.0258  91.52677 18.231487
#> 3  Trial 1   2   1  2 48   502.0258  87.47753 17.424907
#> 4  Trial 1   1   3  3 32   753.0388 109.57170 14.550606
#> 5  Trial 1   1   4  4 24  1004.0517 145.29792 14.471160
#> 6  Trial 1   2   2  4 24  1004.0517 135.13136 13.458606
#> 7  Trial 1   4   1  4 24  1004.0517 128.77483 12.825519
#> 8  Trial 1   1   6  6 16  1506.0775 170.84259 11.343546
#> 9  Trial 1   2   3  6 16  1506.0775 177.45778 11.782779
#> 10 Trial 1   2   4  8 12  2008.1033 210.36242 10.475677
#> 11 Trial 1   4   2  8 12  2008.1033 194.66600  9.694023
#> 12 Trial 1   8   1  8 12  2008.1033 173.34725  8.632387
#> 13 Trial 1   1  12 12  8  3012.1550 252.44911  8.381013
#> 14 Trial 1   2   6 12  8  3012.1550 294.88366  9.789790
#> 15 Trial 1   4   3 12  8  3012.1550 247.35836  8.212006
#> 16 Trial 1   4   4 16  6  4016.2067 290.71099  7.238447
#> 17 Trial 1   8   2 16  6  4016.2067 250.68705  6.241886
#> 18 Trial 1   2  12 24  4  6024.3100 453.63761  7.530117
#> 19 Trial 1   4   6 24  4  6024.3100 385.55772  6.400031
#> 20 Trial 1   8   3 24  4  6024.3100 404.50742  6.714585
#> 21 Trial 1   8   4 32  3  8032.4133 393.61571  4.900342
#> 22 Trial 1   4  12 48  2 12048.6200 745.99765  6.191561
#> 23 Trial 1   8   6 48  2 12048.6200 579.14874  4.806764
```

The columns are `trial`, `X_L` and `X_C` (the shape in basic units), `x`
(the plot size $`X_L X_C`$), `n` (number of plots), `mean`, `sd` and
`cv` (percent), ordered by plot size.

### Same area, different orientation

The repeated plot sizes make the anisotropy visible directly:

``` r

tab[tab$x == 2, ]
#>     trial X_L X_C x  n     mean       sd       cv
#> 2 Trial 1   1   2 2 48 502.0258 91.52677 18.23149
#> 3 Trial 1   2   1 2 48 502.0258 87.47753 17.42491
```

Two shapes of 2 m², one running along the rows and one along the
columns, with different CVs. If they disagree strongly, a plot size
chosen from one orientation will not be the plot size chosen from the
other.

### n falls as the plot grows

``` r

unique(tab[, c("x", "n")])
#>     x  n
#> 1   1 96
#> 2   2 48
#> 4   3 32
#> 5   4 24
#> 8   6 16
#> 10  8 12
#> 13 12  8
#> 16 16  6
#> 18 24  4
#> 21 32  3
#> 22 48  2
```

The largest plots are supported by a handful of plots only, which is the
case for weighting the fit.

## Feeding the fitters

The returned columns are already named `x`, `cv` and `trial`, so the
table goes straight into the models with no reshaping:

``` r

fit_lrp(tab, x = "x", cv = "cv", step = 0.05)
#> Using x = 'x', cv = 'cv' (single series).
#> Linear Response Plateau (LRP) fit
#> Method:                  segment  
#> Breakpoint (Xo):         9.150 
#> CV at breakpoint:        6.960 
#> R2: 0.893  RMSE: 1.513  AIC: 92.3  BIC: 96.9 
#> 
#> Local minima of the SSE profile (8):
#>   Xo =   7.400   SSE  +15.3% vs optimum
#>   Xo =  12.800   SSE  +34.5% vs optimum
#>   Xo =   5.750   SSE  +50.1% vs optimum
#>   ... see $local_minima for all
```

Passing `weights = TRUE` to the fitters uses the `n` column, giving the
well-estimated small plots more say than the sparse large ones.

## Several trials at once

A named list of grids produces one block of rows per trial, ready to fit
one model per trial with the `trial` argument:

``` r

grids <- lapply(split(uniformity_trial, uniformity_trial$trial),
                function(d) as.matrix(d[, grep("^col", names(d))]))

tab3 <- calc_cv_shapes(grids)
#> CV by shape for 3 trial(s).
table(tab3$trial)
#> 
#> T1 T2 T3 
#> 23 23 23
```

``` r

fit_lrp(tab3, x = "x", cv = "cv", trial = "trial", step = 0.05)$summary
#> Using x = 'x', cv = 'cv', trial = 'trial' -> 3 trials.
#>   trial       a       b breakpoint plateau     R2   RMSE     AIC     BIC
#> 1    T1 21.0495 -1.5398       9.15  6.9602 0.8926 1.5131  92.322  96.864
#> 2    T2 19.3597 -1.0725      14.40  3.9162 0.8341 2.3486 112.547 117.089
#> 3    T3 21.3271 -1.5075      10.10  6.1011 0.8019 2.4138 113.806 118.348
#>   n_local
#> 1       8
#> 2       7
#> 3       8
```

## Long-format data

When the measurements come out of a spreadsheet, the grid can be rebuilt
from row and column index columns rather than trusting the row order of
the file:

``` r

long <- expand.grid(col = 1:12, row = 1:8)
long$mf <- as.vector(t(grid1))

head(calc_cv_shapes(long, value = "mf", row_id = "row", col_id = "col"), 4)
#> CV by shape for 1 trial(s).
#>     trial X_L X_C x  n     mean        sd       cv
#> 1 Trial 1   1   1 1 96 251.0129  58.92214 23.47375
#> 2 Trial 1   1   2 2 48 502.0258  91.52677 18.23149
#> 3 Trial 1   2   1 2 48 502.0258  87.47753 17.42491
#> 4 Trial 1   1   3 3 32 753.0388 109.57170 14.55061
```

Missing cells are an error rather than a silent gap: every row/column
combination must be present.

## Where this fits

[`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md)
sits between checking the trial
([`vignette("check_trial")`](https://willyanjnr.github.io/trialSizing/articles/check_trial.md))
and fitting a model
([`vignette("lrp")`](https://willyanjnr.github.io/trialSizing/articles/lrp.md),
[`vignette("qrp")`](https://willyanjnr.github.io/trialSizing/articles/qrp.md),
[`vignette("mcm")`](https://willyanjnr.github.io/trialSizing/articles/mcm.md)).
The Paranaíba method
([`vignette("paranaiba")`](https://willyanjnr.github.io/trialSizing/articles/paranaiba.md))
skips it, working on the raw grid directly, and
[`vignette("compare")`](https://willyanjnr.github.io/trialSizing/articles/compare.md)
runs several methods on one trial at once, building this table on the
way when given a grid.

### References

Cargnelutti Filho, A. et al. (2025). Determinação do tamanho de parcela
para avaliar a massa de parte aérea de grão-de-bico. *Revista
Vivências*, 21(43), 499-513.

Paranaíba, P. F., Ferreira, D. F. & Morais, A. R. (2009). Tamanho ótimo
de parcelas experimentais. *Revista Brasileira de Biometria*, 27(2),
255-268.
