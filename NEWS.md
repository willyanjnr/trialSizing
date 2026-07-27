# trialSizing 0.1.0

* First release.

* Optimal plot size from a uniformity trial by four methods: the modified
  maximum curvature method (`fit_mcm()`), the linear and quadratic response
  plateau models (`fit_lrp()`, `fit_qrp()`), and the closed form of Paranaiba,
  Ferreira and Morais (2009) (`calc_paranaiba()`). `compare_methods()` puts
  them side by side on the same trial.

* `check_trial()` tests a uniformity trial for the spatial structure the
  sizing methods assume, reporting a semivariogram, Moran's I and a kriged
  field map of the basic units.

* `calc_cv_shapes()` groups a grid of basic units into every plot shape the
  grid admits and returns one coefficient of variation per shape, the input
  the plateau models expect.

* `calc_replicates()` turns the coefficient of variation at the optimal plot
  size into the number of replications needed to detect a given difference
  between treatment means, for completely randomised and randomised complete
  block designs.

* Fits carry standardised diagnostic statistics, optional bootstrap
  uncertainty for the breakpoint, and `print()`, `summary()`, `predict()` and
  `plot()` methods.
