# Simulated uniformity trial

A simulated uniformity trial used throughout the package examples and
vignettes. A uniformity trial is a field sown uniformly – one genotype,
one management – and harvested in a fine grid of small *basic
experimental units* (BEU). The spatial variation that remains is
environmental noise, and it is what plot-size methods use to decide how
large a plot must be.

## Usage

``` r
uniformity_trial
```

## Format

A data frame with 24 rows (8 grid rows \\\times\\ 3 trials) and 14
variables:

- trial:

  trial identifier, `"T1"`, `"T2"` or `"T3"`.

- row:

  grid row index, 1-8.

- col01:

  response (g m\\^{-2}\\) at grid column 1.

- col02:

  response (g m\\^{-2}\\) at grid column 2.

- col03:

  response (g m\\^{-2}\\) at grid column 3.

- col04:

  response (g m\\^{-2}\\) at grid column 4.

- col05:

  response (g m\\^{-2}\\) at grid column 5.

- col06:

  response (g m\\^{-2}\\) at grid column 6.

- col07:

  response (g m\\^{-2}\\) at grid column 7.

- col08:

  response (g m\\^{-2}\\) at grid column 8.

- col09:

  response (g m\\^{-2}\\) at grid column 9.

- col10:

  response (g m\\^{-2}\\) at grid column 10.

- col11:

  response (g m\\^{-2}\\) at grid column 11.

- col12:

  response (g m\\^{-2}\\) at grid column 12.

To recover a trial as a numeric matrix (the form the grid functions
expect):


    g1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
                                      grep("^col", names(uniformity_trial))])

## Source

Simulated data for teaching; see `data-raw/uniformity_trial.R`.

## Details

The data are entirely synthetic (a separable first-order autoregressive
field plus independent noise), so they carry no usage restriction. Three
independent trials are provided, each an 8 (rows) by 12 (columns) grid
of 1 m\\^2\\ BEU whose response is a biomass-like measurement in g
m\\^{-2}\\. The parameters were chosen so the coefficient of variation
falls with plot size and levels off at a plateau, as in a real trial.
The generating script is in `data-raw/uniformity_trial.R`.

## See also

[`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md),
[`calc_paranaiba()`](https://willyanjnr.github.io/trialSizing/reference/calc_paranaiba.md),
[`check_trial()`](https://willyanjnr.github.io/trialSizing/reference/check_trial.md),
[`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md)
