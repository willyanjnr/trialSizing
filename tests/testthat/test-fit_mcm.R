## Regression tests for fit_mcm()
## References:
##   Meier & Lessman (1971), Crop Science 11(5), 648-650 (original method);
##   Cargnelutti Filho et al. (2025), Revista Vivencias 21(43), 499-513,
##   Table 2 (chickpea, MMC method) for the default estimator.

test_that("fit_mcm reproduces the published chickpea optima", {
  pub <- chickpea_published$MCM

  ## the article rounds intermediate quantities, so trial 1 lands one
  ## hundredth away (5.7549 -> 5.75 here, 5.76 as printed); a tolerance of
  ## 0.011 accepts that while still catching any real drift
  for (i in seq_along(chickpea_cv_list)) {
    fit <- fit_mcm(chickpea_X, chickpea_cv_list[[i]])
    expect_equal(unname(fit$parameters["Breakpoint"]), pub$Xo[i],
                 tolerance = 0.011, info = paste("trial", i))
    expect_equal(unname(fit$parameters["Breakpoint_Response"]), pub$CVxo[i],
                 tolerance = 0.011, info = paste("trial", i))
  }
})

test_that("the published a and b coefficients are reproduced", {
  ## Table 2 of the article: a = 30.692, b = 0.511 for trial 1
  fit <- fit_mcm(chickpea_X, chickpea_CV1)
  expect_equal(round(unname(fit$coefficients["a"]), 2), 30.69)
  expect_equal(round(unname(fit$coefficients["b"]), 3), 0.511)
})

test_that("the mean optimum matches the article's overall value", {
  xo <- vapply(chickpea_cv_list,
               function(y) unname(fit_mcm(chickpea_X, y)$parameters["Breakpoint"]),
               numeric(1))
  expect_equal(round(mean(xo), 2), chickpea_published$MCM$mean_Xo)
})

test_that("estimates are stable to three decimals", {
  fit <- fit_mcm(chickpea_X, chickpea_CV1)
  expect_equal(unname(fit$parameters["Breakpoint"]), 5.755, tolerance = 1e-3)
  expect_equal(unname(fit$parameters["Breakpoint_Response"]), 12.542,
               tolerance = 1e-3)
  expect_equal(unname(fit$parameters["R2"]), 0.775, tolerance = 1e-3)
})

test_that("the maximum-curvature formula holds", {
  ## Xo = [a^2 b^2 (2b + 1) / (b + 2)]^(1 / (2b + 2))
  fit <- fit_mcm(chickpea_X, chickpea_CV1)
  a <- unname(fit$coefficients["a"])
  b <- unname(fit$coefficients["b"])

  expect_equal(((a^2 * b^2 * (2 * b + 1)) / (b + 2))^(1 / (2 * b + 2)),
               unname(fit$parameters["Breakpoint"]))

  ## CVxo = a * Xo^(-b)
  expect_equal(a * unname(fit$parameters["Breakpoint"])^(-b),
               unname(fit$parameters["Breakpoint_Response"]))
})

test_that("the crambe paper's worked example is reproduced", {
  ## Meier & Lessman (1971): a' = 21.08, b' = 0.2678 give Xc = 3.36 basic units
  a <- 21.08; b <- 0.2678
  xc <- ((a^2 * b^2 * (2 * b + 1)) / (b + 2))^(1 / (2 * b + 2))
  expect_equal(round(xc, 2), 3.36)
})

test_that("the object has the documented structure", {
  fit <- fit_mcm(chickpea_X, chickpea_CV1)

  expect_s3_class(fit, "mcm_fit")
  expect_named(fit$coefficients, c("a", "b"))
  expect_equal(fit$method, "nls")
  expect_false(fit$weighted)

  ## AIC/BIC are deliberately absent: the log-scale option would not be
  ## comparable with the plateau models
  expect_false(any(c("AIC", "BIC") %in% names(fit$parameters)))
})

test_that("the power curve decreases and is never flat", {
  fit <- fit_mcm(chickpea_X, chickpea_CV1)

  expect_gt(unname(fit$coefficients["b"]), 0)
  expect_gt(predict(fit, 1), predict(fit, 6))
  expect_gt(predict(fit, 6), predict(fit, 18))

  ## unlike the plateau models, it keeps falling past the optimum
  xo <- unname(fit$parameters["Breakpoint"])
  expect_gt(predict(fit, xo), predict(fit, xo + 5))
})

test_that("method = 'loglinear' differs from the nls default", {
  nls_fit <- fit_mcm(chickpea_X, chickpea_CV1, method = "nls")
  log_fit <- fit_mcm(chickpea_X, chickpea_CV1, method = "loglinear")

  expect_equal(log_fit$method, "loglinear")
  expect_equal(unname(log_fit$parameters["Breakpoint"]), 6.026,
               tolerance = 1e-3)
  expect_false(isTRUE(all.equal(unname(nls_fit$parameters["Breakpoint"]),
                                unname(log_fit$parameters["Breakpoint"]))))
})

test_that("the Federer degrees-of-freedom weighting changes the fit", {
  unweighted <- fit_mcm(chickpea_X, chickpea_CV1, method = "loglinear")
  weighted   <- fit_mcm(chickpea_X, chickpea_CV1, method = "loglinear",
                        df = chickpea_n - 1)

  expect_true(weighted$weighted)
  expect_equal(unname(weighted$parameters["Breakpoint"]), 5.859,
               tolerance = 1e-3)
  expect_false(isTRUE(all.equal(unname(unweighted$parameters["Breakpoint"]),
                                unname(weighted$parameters["Breakpoint"]))))
})

test_that("MCM gives a smaller optimum than LRP and QRP", {
  for (y in chickpea_cv_list) {
    mcm <- unname(fit_mcm(chickpea_X, y)$parameters["Breakpoint"])
    lrp <- unname(fit_lrp(chickpea_X, y, step = 0.01)$parameters["Breakpoint"])
    qrp <- unname(fit_qrp(chickpea_X, y, step = 0.01)$parameters["Breakpoint"])
    expect_lt(mcm, lrp)
    expect_lt(lrp, qrp)
  }
})

test_that("the data-frame interface fits one model per trial", {
  res <- fit_mcm(chickpea_long, x = "x", cv = "cv", trial = "trial")

  expect_s3_class(res, "mcm_multi")
  expect_equal(nrow(res$summary), 3)
  expect_equal(res$summary$breakpoint, chickpea_published$MCM$Xo,
               tolerance = 0.011)
})

test_that("df can be supplied as a column in the data-frame interface", {
  long <- chickpea_long
  long$df <- rep(chickpea_n - 1, 3)

  res <- fit_mcm(long, x = "x", cv = "cv", df = "df", trial = "trial",
                 method = "loglinear")
  expect_true(all(res$summary$weighted))
  expect_equal(round(res$summary$breakpoint[1], 3), 5.859)
})

test_that("invalid input is rejected", {
  expect_error(fit_mcm(chickpea_X, chickpea_CV1[-1]), "same length")
  expect_error(fit_mcm(c(0, chickpea_X[-1]), chickpea_CV1),
               "strictly positive")
  expect_error(fit_mcm(chickpea_X, c(0, chickpea_CV1[-1])),
               "strictly positive")
  expect_error(fit_mcm(chickpea_X, chickpea_CV1, df = c(1, 2)),
               "matching `x` in length")
  expect_error(fit_mcm(chickpea_long, x = "x", cv = "CV"), "not found")
})

test_that("plot, print and summary work", {
  skip_if_not_installed("ggplot2")

  fit <- fit_mcm(chickpea_X, chickpea_CV1)
  expect_s3_class(plot(fit, title = "Trial 1"), "ggplot")
  expect_s3_class(plot(fit, decimal_mark = ","), "ggplot")
  expect_output(print(fit), "Modified Maximum Curvature")
  expect_output(summary(fit), "Model coefficients")
})
