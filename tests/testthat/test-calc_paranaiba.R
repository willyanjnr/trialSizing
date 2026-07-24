## Regression tests for calc_paranaiba()
## References:
##   Paranaiba, Ferreira & Morais (2009), Revista Brasileira de Biometria
##   27(2), 255-268 (original method);
##   Cargnelutti Filho et al. (2014), Ciencia Rural 44(10), 1732-1739,
##   Table 1, which prints the formulas and the values reproduced below.

test_that("the published formulas reproduce Table 1 of the black oat article", {
  ## Xo = 10 * (2 (1 - rho^2) s^2 m)^(1/3) / m
  xo <- 10 * (2 * (1 - oat_rho^2) * oat_s2 * oat_mean)^(1 / 3) / oat_mean
  expect_equal(round(xo, 2), oat_Xo_published, tolerance = 0.011)

  ## CVxo = 100 * sqrt((1 - rho^2) s^2 / m^2) / sqrt(Xo)
  cvxo <- 100 * sqrt((1 - oat_rho^2) * oat_s2 / oat_mean^2) / sqrt(xo)
  expect_equal(round(cvxo, 2), oat_CVxo_published, tolerance = 0.031)

  ## the article's reported mean for the first evaluation date
  expect_equal(round(mean(xo), 2), 4.12)
})

test_that("calc_paranaiba estimates are stable to three decimals", {
  fit <- calc_paranaiba(chickpea_grid1)

  expect_equal(fit$summary$mean, 227.687, tolerance = 1e-3)
  expect_equal(fit$summary$variance, 4792.240, tolerance = 1e-3)
  expect_equal(fit$summary$rho_row, -0.018972, tolerance = 1e-4)
  expect_equal(fit$summary$rho_col,  0.129026, tolerance = 1e-4)
  expect_equal(fit$summary$Xo, 5.696, tolerance = 1e-3)
  expect_equal(fit$summary$CVxo, 12.737, tolerance = 1e-3)
})

test_that("the closed forms hold for the fitted object", {
  fit <- calc_paranaiba(chickpea_grid1)
  s <- fit$summary

  expect_equal(s$Xo,
               10 * (2 * (1 - s$rho^2) * s$variance * s$mean)^(1 / 3) / s$mean)
  expect_equal(s$CVxo,
               100 * sqrt((1 - s$rho^2) * s$variance / s$mean^2) / sqrt(s$Xo))

  ## the CV of the trial is the usual one
  expect_equal(s$CV, 100 * sqrt(s$variance) / s$mean)
})

test_that("the object has the documented structure", {
  fit <- calc_paranaiba(chickpea_grid1)

  expect_s3_class(fit, "paranaiba_fit")
  expect_true(all(c("trial", "mean", "variance", "CV", "rho_row", "rho_col",
                    "rho", "Xo", "CVxo", "valid") %in% names(fit$summary)))
  expect_true(fit$summary$valid)
  expect_equal(fit$meta$rho_direction, "row")
  expect_equal(fit$matrices[[1]], chickpea_grid1)
})

test_that("rho_direction selects which walk is used", {
  by_row  <- calc_paranaiba(chickpea_grid1, rho_direction = "row")
  by_col  <- calc_paranaiba(chickpea_grid1, rho_direction = "col")
  by_mean <- calc_paranaiba(chickpea_grid1, rho_direction = "mean")

  ## both directions are always reported, whichever is selected
  expect_equal(by_row$summary$rho_row, by_col$summary$rho_row)
  expect_equal(by_row$summary$rho_col, by_col$summary$rho_col)

  ## the selected one drives the estimate
  expect_equal(by_row$summary$rho,  by_row$summary$rho_row)
  expect_equal(by_col$summary$rho,  by_col$summary$rho_col)
  expect_equal(by_mean$summary$rho,
               mean(c(by_row$summary$rho_row, by_row$summary$rho_col)))

  expect_equal(by_col$summary$Xo, 5.665, tolerance = 1e-3)
  expect_equal(by_mean$summary$Xo, 5.691, tolerance = 1e-3)
})

test_that("Xo is largest when there is no spatial autocorrelation", {
  ## rho enters only through (1 - rho^2), so |rho| > 0 shrinks Xo
  m <- 250; s2 <- 5000
  xo <- function(rho) 10 * (2 * (1 - rho^2) * s2 * m)^(1 / 3) / m

  expect_gt(xo(0), xo(0.5))
  expect_gt(xo(0), xo(-0.5))
  expect_equal(xo(0.5), xo(-0.5))
})

test_that("a list of matrices gives one row per trial", {
  res <- calc_paranaiba(list(`Trial 1` = chickpea_grid1,
                             `Trial 2` = chickpea_grid2))

  expect_equal(nrow(res$summary), 2)
  expect_equal(res$summary$trial, c("Trial 1", "Trial 2"))
  expect_equal(res$meta$n_trials, 2)

  ## the first trial matches the same grid fitted alone
  alone <- calc_paranaiba(chickpea_grid1)
  expect_equal(res$summary$Xo[1], alone$summary$Xo)
})

test_that("the long data-frame interface rebuilds the grid from the indices", {
  long <- expand.grid(col = 1:6, row = 1:6)
  long$mf <- as.vector(t(chickpea_grid1))

  from_long   <- calc_paranaiba(long, value = "mf", row_id = "row",
                                col_id = "col")
  from_matrix <- calc_paranaiba(chickpea_grid1)

  expect_equal(from_long$summary$Xo,   from_matrix$summary$Xo)
  expect_equal(from_long$summary$rho,  from_matrix$summary$rho)
})

test_that("row order in the data frame does not affect the result", {
  ## the grid comes from row_id/col_id, not from the order of the rows
  long <- expand.grid(col = 1:6, row = 1:6)
  long$mf <- as.vector(t(chickpea_grid1))

  set.seed(1)
  shuffled <- long[sample(nrow(long)), ]

  expect_equal(
    calc_paranaiba(shuffled, value = "mf", row_id = "row",
                   col_id = "col")$summary$Xo,
    calc_paranaiba(long, value = "mf", row_id = "row",
                   col_id = "col")$summary$Xo
  )
})

test_that("a trial column splits the data frame", {
  long1 <- expand.grid(col = 1:6, row = 1:6)
  long1$mf <- as.vector(t(chickpea_grid1)); long1$trial <- "Trial 1"
  long2 <- expand.grid(col = 1:6, row = 1:6)
  long2$mf <- as.vector(t(chickpea_grid2)); long2$trial <- "Trial 2"

  res <- calc_paranaiba(rbind(long1, long2), value = "mf", row_id = "row",
                        col_id = "col", trial = "trial")
  expect_equal(nrow(res$summary), 2)
  expect_equal(res$summary$Xo[1], calc_paranaiba(chickpea_grid1)$summary$Xo)
})

test_that("transposing the grid changes the row-wise estimate", {
  ## rho is directional, so the layout must be preserved
  straight  <- calc_paranaiba(chickpea_grid1)
  transposed <- calc_paranaiba(t(chickpea_grid1))

  expect_equal(straight$summary$rho_row, transposed$summary$rho_col)
  expect_false(isTRUE(all.equal(straight$summary$Xo, transposed$summary$Xo)))
})

test_that("degenerate trials are flagged rather than returning NaN", {
  flat <- matrix(100, nrow = 6, ncol = 6)
  expect_warning(res <- calc_paranaiba(flat), "zero variance")

  expect_false(res$summary$valid)
  expect_true(is.na(res$summary$Xo))
  expect_true(is.na(res$summary$CVxo))
})

test_that("invalid input is rejected", {
  expect_error(calc_paranaiba("not a grid"), "matrix, a list of matrices")
  expect_error(calc_paranaiba(data.frame(a = 1:4)), "give `value`")

  long <- expand.grid(col = 1:6, row = 1:6)
  long$mf <- as.vector(t(chickpea_grid1))
  expect_error(calc_paranaiba(long, value = "MF", row_id = "row",
                              col_id = "col"), "not found")
  expect_error(calc_paranaiba(chickpea_grid1, n_row = 8),
               "does not match n_row/n_col")
})

test_that("plot, print and summary work", {
  skip_if_not_installed("ggplot2")

  res <- calc_paranaiba(list(`Trial 1` = chickpea_grid1,
                             `Trial 2` = chickpea_grid2))
  expect_s3_class(plot(res), "ggplot")
  expect_s3_class(plot(res, y_var = "CVxo"), "ggplot")
  expect_output(print(res), "Paranaiba optimal plot size")
  expect_output(summary(res), "Optimal plot size")
})
