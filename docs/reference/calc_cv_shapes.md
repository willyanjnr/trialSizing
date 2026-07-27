# CV by plot shape from a uniformity-trial grid

Builds the coefficient-of-variation table that
[`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md),
[`fit_qrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_qrp.md)
and
[`fit_mcm()`](https://willyanjnr.github.io/trialSizing/reference/fit_mcm.md)
consume, starting from the raw grid of basic experimental units (BEU) as
the trial was harvested. Adjacent BEU are grouped into every rectangular
plot shape the grid admits: for a grid of \\L\\ rows by \\C\\ columns,
each shape is \\X_L \times X_C\\ basic units with \\X_L\\ a divisor of
\\L\\ and \\X_C\\ a divisor of \\C\\. For each shape the plot totals are
formed and their mean, standard deviation and CV reported.

## Usage

``` r
calc_cv_shapes(
  .data,
  value = NULL,
  row_id = NULL,
  col_id = NULL,
  trial = NULL,
  min_plots = 2
)
```

## Arguments

- .data:

  a matrix (one trial), a named list of matrices (several trials), or a
  long data frame with one row per basic unit.

- value:

  name of the column holding the BEU measurement (long format).

- row_id, col_id:

  names of the row and column index columns (long format).

- trial:

  optional column name identifying the trial (long format).

- min_plots:

  minimum number of plots a shape must yield to be kept (default 2).

## Value

A data frame with one row per shape and trial: `trial`, `X_L`, `X_C`
(shape in basic units), `x` (plot size), `n` (number of plots), `mean`,
`sd` and `cv` (percent), ordered by plot size.

## Details

This closes the pipeline: the grid goes in, the CV table comes out, and
it feeds the fitters directly, since the returned columns are already
named `x`, `cv` and `trial`.

## Shapes and plot counts

A shape of area \\X = X_L X_C\\ yields \\n = LC/X\\ plots, so the CV of
a large plot size rests on few plots and is estimated with little
information: the last rows of the table are the least reliable. The `n`
column is returned for exactly this reason and is the natural weight for
a weighted fit. Shapes leaving fewer than `min_plots` plots are dropped,
which by default removes only the whole grid (a single plot, no
variance).

Shapes are not interchangeable at equal area: \\1 \times 2\\ and \\2
\times 1\\ cover the same two basic units but run along different
directions of the field, and their CVs differ whenever fertility is not
isotropic. Both are reported, which is why the CV table has repeated `x`
values.

## References

Cargnelutti Filho, A. et al. (2025). Determinacao do tamanho de parcela
para avaliar a massa de parte aerea de grao-de-bico. *Revista
Vivencias*, 21(43), 499-513.  
Paranaiba, P. F., Ferreira, D. F. & Morais, A. R. (2009). *Revista
Brasileira de Biometria*, 27(2), 255-268.

## See also

[`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md),
[`fit_qrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_qrp.md),
[`fit_mcm()`](https://willyanjnr.github.io/trialSizing/reference/fit_mcm.md),
[`calc_paranaiba()`](https://willyanjnr.github.io/trialSizing/reference/calc_paranaiba.md)

## Examples

``` r
## Three trials of the bundled simulated uniformity trial, each an 8 x 12
## grid of basic units (see ?uniformity_trial).
grids <- lapply(split(uniformity_trial, uniformity_trial$trial),
                function(d) as.matrix(d[, grep("^col", names(d))]))

tab <- calc_cv_shapes(grids)
#> CV by shape for 3 trial(s).
head(tab, 8)
#>   trial X_L X_C x  n      mean        sd       cv
#> 1    T1   1   1 1 96  251.0129  58.92214 23.47375
#> 2    T1   1   2 2 48  502.0258  91.52677 18.23149
#> 3    T1   2   1 2 48  502.0258  87.47753 17.42491
#> 4    T1   1   3 3 32  753.0388 109.57170 14.55061
#> 5    T1   1   4 4 24 1004.0517 145.29792 14.47116
#> 6    T1   2   2 4 24 1004.0517 135.13136 13.45861
#> 7    T1   4   1 4 24 1004.0517 128.77483 12.82552
#> 8    T1   1   6 6 16 1506.0775 170.84259 11.34355

## Same area, different orientation, different CV
tab[tab$trial == "T1" & tab$x == 2, ]
#>   trial X_L X_C x  n     mean       sd       cv
#> 2    T1   1   2 2 48 502.0258 91.52677 18.23149
#> 3    T1   2   1 2 48 502.0258 87.47753 17.42491

## The table feeds the fitters unchanged
fit_lrp(tab[tab$trial == "T1", ], x = "x", cv = "cv", step = 0.05)
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

# \donttest{
## One model per trial, straight from the grids
fit_lrp(tab, x = "x", cv = "cv", trial = "trial", step = 0.01)$summary
#> Using x = 'x', cv = 'cv', trial = 'trial' -> 3 trials.
#>   trial       a       b breakpoint plateau     R2   RMSE     AIC     BIC
#> 1    T1 21.0495 -1.5398       9.16  6.9448 0.8926 1.5130  92.321  96.863
#> 2    T2 19.3597 -1.0725      14.40  3.9162 0.8341 2.3486 112.547 117.089
#> 3    T3 21.3271 -1.5075      10.11  6.0860 0.8019 2.4137 113.805 118.347
#>   n_local
#> 1       8
#> 2       7
#> 3       8

## n falls as the plot grows: the last shapes rest on very few plots
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
# }
```
