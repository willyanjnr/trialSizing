# Changelog

## trialSizing 0.1.0

- First release.

- Optimal plot size from a uniformity trial by four methods: the
  modified maximum curvature method
  ([`fit_mcm()`](https://willyanjnr.github.io/trialSizing/reference/fit_mcm.md)),
  the linear and quadratic response plateau models
  ([`fit_lrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_lrp.md),
  [`fit_qrp()`](https://willyanjnr.github.io/trialSizing/reference/fit_qrp.md)),
  and the closed form of Paranaiba, Ferreira and Morais (2009)
  ([`calc_paranaiba()`](https://willyanjnr.github.io/trialSizing/reference/calc_paranaiba.md)).
  [`compare_methods()`](https://willyanjnr.github.io/trialSizing/reference/compare_methods.md)
  puts them side by side on the same trial.

- [`check_trial()`](https://willyanjnr.github.io/trialSizing/reference/check_trial.md)
  tests a uniformity trial for the spatial structure the sizing methods
  assume, reporting a semivariogram, Moran’s I and a kriged field map of
  the basic units.

- [`calc_cv_shapes()`](https://willyanjnr.github.io/trialSizing/reference/calc_cv_shapes.md)
  groups a grid of basic units into every plot shape the grid admits and
  returns one coefficient of variation per shape, the input the plateau
  models expect.

- [`calc_replicates()`](https://willyanjnr.github.io/trialSizing/reference/calc_replicates.md)
  turns the coefficient of variation at the optimal plot size into the
  number of replications needed to detect a given difference between
  treatment means, for completely randomised and randomised complete
  block designs.

- Fits carry standardised diagnostic statistics, optional bootstrap
  uncertainty for the breakpoint, and
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
