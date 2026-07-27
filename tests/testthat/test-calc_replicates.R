## Regression tests for calc_replicates()
## Reference: Cargnelutti Filho et al. (2014), Ciencia Rural 44(10), 1732-1739,
## Tables 2 (CRD) and 3 (RCBD), computed with CVxo = 9.25%.

test_that("Table 2 (CRD) is reproduced", {
  ## selected cells, d = 10%
  published <- c(`3` = 10.46, `4` = 12.18, `5` = 13.50, `6` = 14.58,
                 `7` = 15.49, `8` = 16.28, `9` = 16.98, `10` = 17.60,
                 `25` = 23.13, `50` = 27.41)

  for (i in names(published)) {
    r <- calc_replicates(as.numeric(i), oat_CVxo, 10, design = "CRD")
    expect_equal(round(r$data$r_continuous, 2), unname(published[i]),
                 info = paste("i =", i))
  }
})

test_that("Table 3 (RCBD) is reproduced", {
  published <- c(`3` = 10.95, `10` = 17.65, `50` = 27.42)

  for (i in names(published)) {
    r <- calc_replicates(as.numeric(i), oat_CVxo, 10, design = "RCBD")
    expect_equal(round(r$data$r_continuous, 2), unname(published[i]),
                 tolerance = 0.011, info = paste("i =", i))
  }
})

test_that("other columns of the published tables are reproduced", {
  ## d = 12%, 20% and 30% for i = 10 and i = 50 (CRD)
  expect_equal(round(calc_replicates(10, oat_CVxo, 12)$data$r_continuous, 2),
               12.38)
  expect_equal(round(calc_replicates(50, oat_CVxo, 20)$data$r_continuous, 2),
               6.97)
  expect_equal(round(calc_replicates(50, oat_CVxo, 30)$data$r_continuous, 2),
               3.21)
  expect_equal(round(calc_replicates(50, oat_CVxo, 30,
                                     design = "RCBD")$data$r_continuous, 2),
               3.22)
})

test_that("the fixed point satisfies its own equation", {
  ## r = (q(df(r)) * CV / d)^2 must hold at the returned r_continuous
  for (i in c(3, 10, 50)) {
    for (d in c(10, 20, 30)) {
      r  <- calc_replicates(i, oat_CVxo, d, design = "CRD")$data
      df <- i * (r$r_continuous - 1)
      q  <- stats::qtukey(0.95, nmeans = i, df = df)
      expect_equal(r$r_continuous, (q * oat_CVxo / d)^2, tolerance = 1e-5,
                   info = paste("i =", i, "d =", d))
    }
  }
})

test_that("r_optimal is the ceiling of r_continuous, floored at two", {
  res <- calc_replicates(c(3, 10, 50), oat_CVxo, c(10, 20, 30))
  expect_equal(res$data$r_optimal,
               pmax(2, ceiling(res$data$r_continuous)))

  ## a large LSD drives the requirement below 2, where the floor applies:
  ## the result is reported at the floor and flagged, never as NaN
  easy <- calc_replicates(3, 5, 100)$data
  expect_true(easy$at_floor)
  expect_equal(easy$r_continuous, 2)
  expect_equal(easy$r_optimal, 2)
  expect_gt(easy$df_error, 0)
  expect_false(is.nan(easy$q_tukey))
  expect_true(easy$converged)
})

test_that("the object has the documented structure", {
  res <- calc_replicates(c(3, 10), oat_CVxo, c(10, 20))

  expect_s3_class(res, "replicates_fit")
  expect_equal(nrow(res$data), 4)
  expect_true(all(c("Treatments", "CV_percent", "LSD_percent", "Alpha",
                    "Design", "r_continuous", "r_optimal", "df_error",
                    "q_tukey", "converged", "at_floor") %in% names(res$data)))
  expect_true(all(res$data$converged))
  expect_false(any(res$data$at_floor))
  expect_equal(res$meta$design, "CRD")
})

test_that("error degrees of freedom follow the design", {
  crd  <- calc_replicates(10, oat_CVxo, 20, design = "CRD")$data
  rcbd <- calc_replicates(10, oat_CVxo, 20, design = "RCBD")$data

  expect_equal(crd$df_error, 10 * (crd$r_optimal - 1))
  expect_equal(rcbd$df_error, (10 - 1) * (rcbd$r_optimal - 1))
})

test_that("RCBD needs at least as many replications as CRD", {
  for (i in c(3, 10, 50)) {
    crd  <- calc_replicates(i, oat_CVxo, 10, design = "CRD")$data
    rcbd <- calc_replicates(i, oat_CVxo, 10, design = "RCBD")$data
    expect_gte(rcbd$r_continuous, crd$r_continuous)
  }
})

test_that("the requirement moves in the expected direction", {
  ## more treatments -> more replications
  r <- calc_replicates(c(3, 10, 50), oat_CVxo, 10)$data
  expect_true(all(diff(r$r_continuous) > 0))

  ## a larger detectable difference -> fewer replications
  d <- calc_replicates(10, oat_CVxo, c(10, 20, 30))$data
  expect_true(all(diff(d$r_continuous) < 0))

  ## a higher CV -> more replications
  expect_gt(calc_replicates(10, 20, 10)$data$r_continuous,
            calc_replicates(10, 10, 10)$data$r_continuous)

  ## a stricter alpha -> more replications
  expect_gt(calc_replicates(10, oat_CVxo, 20, alpha = 0.01)$data$r_continuous,
            calc_replicates(10, oat_CVxo, 20, alpha = 0.05)$data$r_continuous)
})

test_that("the grid covers every combination", {
  res <- calc_replicates(3:6, oat_CVxo, c(10, 20, 30))
  expect_equal(nrow(res$data), 4 * 3)
  expect_equal(sort(unique(res$data$Treatments)), 3:6)
  expect_equal(sort(unique(res$data$LSD_percent)), c(10, 20, 30))
})

test_that("invalid input is rejected", {
  expect_error(calc_replicates(1, oat_CVxo, 10), "2 or greater")
  expect_error(calc_replicates(10, -1, 10), "single positive number")
  expect_error(calc_replicates(10, c(9, 10), 10), "single positive number")
  expect_error(calc_replicates(10, oat_CVxo, -5), "strictly positive")
  expect_error(calc_replicates(10, oat_CVxo, 10, alpha = 0), "between 0 and 1")
  expect_error(calc_replicates(10, oat_CVxo, 10, alpha = 1), "between 0 and 1")
})

test_that("the CVxo of a plot-size fit can be fed straight in", {
  cvxo <- unname(fit_lrp(chickpea_X, chickpea_CV1,
                         step = 0.01)$parameters["Breakpoint_Response"])
  res <- calc_replicates(10, cvxo, 20, design = "RCBD")

  expect_equal(res$data$CV_percent, cvxo)
  expect_true(res$data$converged)
  expect_gte(res$data$r_optimal, 2)
})

test_that("plot, print and summary work", {
  skip_if_not_installed("ggplot2")

  res <- calc_replicates(3:10, oat_CVxo, c(10, 20, 30))
  expect_s3_class(plot(res), "ggplot")
  expect_s3_class(plot(res, y_var = "r_continuous"), "ggplot")
  expect_output(print(res), "Optimal number of replications")
  expect_output(summary(res), "by LSD")
})
