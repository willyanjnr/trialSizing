## Regression tests for check_trial() and its field-map plot method.
## The geostatistics is written from the standard definitions, so the values
## below are anchored to the current implementation on chickpea trial 1.

test_that("the object has the documented structure", {
  chk <- check_trial(chickpea_grid1)

  expect_s3_class(chk, "trial_check")
  expect_equal(names(chk$summary),
               c("trial", "rows", "cols", "n_missing", "n_shapes", "mean",
                 "cv", "n_outliers", "morans_i", "range", "nugget_ratio",
                 "dependence", "n_issues"))
  expect_equal(chk$meta$n_trials, 1)
  expect_equal(chk$meta$cell_size, c(1, 1))
  expect_equal(chk$checks[[1]]$matrix, chickpea_grid1)
})

test_that("structure and value summaries match the grid", {
  k <- check_trial(chickpea_grid1)$checks[[1]]

  expect_equal(k$dim[["rows"]], 6)
  expect_equal(k$dim[["cols"]], 6)
  expect_equal(k$n_cells, 36)
  expect_equal(k$n_missing, 0)
  expect_equal(k$n_shapes, 15)          # same 15 shapes calc_cv_shapes uses

  expect_equal(unname(k$stats[["mean"]]), 227.6872, tolerance = 1e-3)
  expect_equal(unname(k$stats[["sd"]]),   69.2260,  tolerance = 1e-3)
  expect_equal(unname(k$stats[["cv"]]),   30.4040,  tolerance = 1e-3)
  expect_equal(unname(k$stats[["min"]]),  92.88)
  expect_equal(unname(k$stats[["max"]]),  346.70)
  expect_null(k$outliers)
  expect_length(k$issues, 0)
})

test_that("the spatial measures are stable", {
  sp <- check_trial(chickpea_grid1)$checks[[1]]$spatial

  expect_equal(sp$morans_i, -0.0202237, tolerance = 1e-4)
  expect_gt(sp$morans_p, 0.9)          # no significant global autocorrelation
  expect_equal(sp$rho_row, -0.0486832, tolerance = 1e-4)
  expect_equal(sp$rho_col, -0.0586492, tolerance = 1e-4)
})

test_that("the fitted variogram is stable and read for dependence", {
  v <- check_trial(chickpea_grid1)$checks[[1]]$variogram

  expect_true(v$model %in% c("spherical", "exponential", "gaussian"))
  expect_equal(v$model, "gaussian")
  expect_equal(v$sill, v$nugget + v$psill)
  expect_equal(v$range, 4.7434, tolerance = 1e-3)
  expect_equal(v$nugget_ratio, 0.6905, tolerance = 1e-3)
  expect_equal(v$dependence, "moderate")   # 0.25 < ratio <= 0.75
})

test_that("transposing the grid swaps the directional autocorrelations", {
  straight   <- check_trial(chickpea_grid1)$checks[[1]]$spatial
  transposed <- check_trial(t(chickpea_grid1))$checks[[1]]$spatial

  expect_equal(transposed$rho_row, straight$rho_col)
  expect_equal(transposed$rho_col, straight$rho_row)
})

test_that("a constant field is flagged and no variogram is fitted", {
  flat <- matrix(100, nrow = 6, ncol = 6)
  chk  <- check_trial(flat)

  expect_match(chk$checks[[1]]$issues, "identical", all = FALSE)
  expect_null(chk$checks[[1]]$variogram)
  expect_equal(chk$summary$cv, 0)
})

test_that("a grid with prime sides admits almost no plot shapes", {
  set.seed(1)
  chk <- check_trial(matrix(rnorm(77, 100, 10), nrow = 7))

  expect_equal(chk$checks[[1]]$n_shapes, 3)
  expect_match(chk$checks[[1]]$issues, "plot shape", all = FALSE)
})

test_that("outliers are reported by the boxplot rule", {
  spiked <- chickpea_grid1
  spiked[1, 1] <- 5000                 # a clear harvest-failure spike
  k <- check_trial(spiked)$checks[[1]]

  expect_false(is.null(k$outliers))
  expect_true(any(k$outliers$value == 5000))
})

test_that("variogram = FALSE skips the fit", {
  chk <- check_trial(chickpea_grid1, variogram = FALSE)

  expect_null(chk$checks[[1]]$variogram)
  expect_true(is.na(chk$summary$range))
  expect_true(is.na(chk$summary$nugget_ratio))
})

test_that("a list of grids gives one row per trial", {
  chk <- check_trial(list(A = chickpea_grid1, B = chickpea_grid2))

  expect_equal(nrow(chk$summary), 2)
  expect_equal(chk$summary$trial, c("A", "B"))
  expect_equal(chk$summary$cv[1],
               check_trial(chickpea_grid1)$summary$cv)
})

test_that("invalid input is rejected", {
  expect_error(check_trial(chickpea_grid1, cell_size = 1), "two positive")
  expect_error(check_trial(chickpea_grid1, cell_size = c(-1, 1)), "two positive")
  expect_error(check_trial(chickpea_grid1, n_bins = 2), "at least 3")
  expect_error(check_trial("not a grid"), "matrix, a list of matrices")
})

test_that("plot, print and summary work", {
  skip_if_not_installed("ggplot2")
  chk <- check_trial(list(A = chickpea_grid1, B = chickpea_grid2))

  expect_s3_class(plot(chk), "ggplot")
  expect_s3_class(plot(chk, palette = "blues", point_values = FALSE), "ggplot")
  expect_s3_class(plot(chk, surface = FALSE), "ggplot")
  expect_error(plot(chk, resolution = 2), "at least 10")
  expect_output(print(chk), "Uniformity trial check")
  expect_output(summary(chk), "nugget_ratio")
})
