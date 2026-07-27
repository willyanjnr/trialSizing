## Regression tests for fit_lrp()
## Reference: Cargnelutti Filho et al. (2025), Revista Vivencias 21(43),
## 499-513, Table 2 (chickpea, LRP method).

test_that("fit_lrp reproduces the published chickpea optima", {
  pub <- chickpea_published$LRP

  ## step = 0.001 (the default) is required here: at step = 0.01 the
  ## breakpoints still round to the published Xo, but CVxo moves by a
  ## hundredth and no longer matches Table 2 exactly
  for (i in seq_along(chickpea_cv_list)) {
    fit <- fit_lrp(chickpea_X, chickpea_cv_list[[i]], step = 0.001)
    expect_equal(round(unname(fit$parameters["Breakpoint"]), 2),
                 pub$Xo[i],
                 info = paste("trial", i))
    expect_equal(round(unname(fit$parameters["Breakpoint_Response"]), 2),
                 pub$CVxo[i],
                 info = paste("trial", i))
  }
})

test_that("the mean optimum matches the article's overall recommendation", {
  xo <- vapply(chickpea_cv_list,
               function(y) unname(fit_lrp(chickpea_X, y, step = 0.01)$parameters["Breakpoint"]),
               numeric(1))
  expect_equal(round(mean(xo), 2), chickpea_published$LRP$mean_Xo)
})

test_that("estimates are stable to three decimals", {
  ## the default step = 0.001 is what pins the third decimal; a coarser
  ## grid cannot resolve it
  fit <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.001)
  expect_equal(unname(fit$parameters["Breakpoint"]), 7.436, tolerance = 1e-6)
  expect_equal(unname(fit$parameters["Breakpoint_Response"]), 7.881,
               tolerance = 1e-3)
  expect_equal(unname(fit$parameters["R2"]), 0.715, tolerance = 1e-3)
  expect_equal(unname(fit$coefficients["a"]), 28.0066, tolerance = 1e-3)
  expect_equal(unname(fit$coefficients["b"]), -2.7066, tolerance = 1e-3)
})

test_that("the object has the documented structure", {
  fit <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.01)

  expect_s3_class(fit, "lrp_fit")
  expect_named(fit$coefficients, c("a", "b"))
  expect_named(fit$parameters,
               c("Breakpoint", "Breakpoint_Response", "R2", "RMSE",
                 "AIC", "BIC"))
  expect_equal(nrow(fit$data), length(chickpea_X))
  expect_equal(fit$method, "segment")
})

test_that("fitted values and residuals are internally consistent", {
  fit <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.01)

  expect_equal(fit$fitted + fit$residuals, chickpea_CV1)
  expect_equal(predict(fit), fit$fitted)

  ## R2 recomputed from the residuals
  r2 <- 1 - sum(fit$residuals^2) /
        sum((chickpea_CV1 - mean(chickpea_CV1))^2)
  expect_equal(unname(fit$parameters["R2"]), r2)
})

test_that("the model is flat past the breakpoint", {
  fit <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.01)
  xo  <- unname(fit$parameters["Breakpoint"])
  plateau <- unname(fit$parameters["Breakpoint_Response"])

  expect_equal(predict(fit, xo + 1), plateau)
  expect_equal(predict(fit, xo + 50), plateau)
  expect_lt(predict(fit, xo - 1), predict(fit, xo - 2))
})

test_that("the vector interface accepts positional and named arguments", {
  positional <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.01)
  named      <- fit_lrp(x = chickpea_X, cv = chickpea_CV1, step = 0.01)
  expect_equal(positional$parameters, named$parameters)
})

test_that("the data-frame interface fits one model per trial", {
  res <- fit_lrp(chickpea_long, x = "x", cv = "cv", trial = "trial", step = 0.01)

  expect_s3_class(res, "lrp_multi")
  expect_equal(nrow(res$summary), 3)
  expect_length(res$fits, 3)
  expect_equal(round(res$summary$breakpoint, 2), chickpea_published$LRP$Xo)

  ## a single trial from the group matches the same data fitted alone
  alone <- fit_lrp(chickpea_X, chickpea_CV2, step = 0.01)
  expect_equal(res$fits[["Trial 2"]]$parameters, alone$parameters)
})

test_that("a data frame without a trial column fits a single series", {
  one <- chickpea_long[chickpea_long$trial == "Trial 1", ]
  res <- fit_lrp(one, x = "x", cv = "cv", step = 0.01)
  expect_s3_class(res, "lrp_fit")
  expect_equal(round(unname(res$parameters["Breakpoint"]), 2), 7.44)
})

test_that("method = 'ramp' runs and stays in a plausible range", {
  ramp <- fit_lrp(chickpea_X, chickpea_CV1, method = "ramp", step = 0.01)
  expect_equal(ramp$method, "ramp")
  expect_gt(unname(ramp$parameters["Breakpoint"]), min(chickpea_X))
  expect_lt(unname(ramp$parameters["Breakpoint"]), max(chickpea_X))
})

test_that("a coarser step still reproduces the published value", {
  coarse <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.01)
  expect_equal(round(unname(coarse$parameters["Breakpoint"]), 2), 7.44)
})

test_that("search_range restricts the breakpoint and is validated", {
  ## a range containing the optimum leaves it unchanged
  inside <- fit_lrp(chickpea_X, chickpea_CV1, search_range = c(5, 12), step = 0.01)
  expect_equal(round(unname(inside$parameters["Breakpoint"]), 2), 7.44)

  ## a restricting range pushes the breakpoint to its edge, with a warning
  expect_warning(
    high <- fit_lrp(chickpea_X, chickpea_CV1, search_range = c(10, 18), step = 0.01),
    "edge of the search range"
  )
  expect_gte(unname(high$parameters["Breakpoint"]), 10)

  ## outside the data range, or inverted, is an error
  expect_error(fit_lrp(chickpea_X, chickpea_CV1, search_range = c(20, 40), step = 0.01),
               "within the data range")
  expect_error(fit_lrp(chickpea_X, chickpea_CV1, search_range = c(12, 5), step = 0.01),
               "lower < upper")
})

test_that("invalid input is rejected with a clear message", {
  expect_error(fit_lrp(chickpea_X, chickpea_CV1[-1]), "same length")
  expect_error(fit_lrp(c(chickpea_X[-1], NA), chickpea_CV1), "missing values")
  expect_error(fit_lrp(1:3, c(1, 2, 3)), "At least 4 observations")
  expect_error(fit_lrp(rep(1, 10), rnorm(10)), "distinct")
  expect_error(fit_lrp(chickpea_long, x = "x", cv = "CV", trial = "trial"),
               "not found")
})

test_that("a non-decreasing series warns about the slope", {
  increasing <- chickpea_X * 1.5 + 3
  expect_warning(fit_lrp(chickpea_X, increasing, step = 0.01), "slope is non-negative")
})

test_that("plot methods return ggplot objects without writing files", {
  skip_if_not_installed("ggplot2")

  fit <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.01)
  g <- plot(fit, title = "Trial 1")
  expect_s3_class(g, "ggplot")

  ## the annotation block can be turned off
  expect_s3_class(plot(fit, annotate_model = FALSE), "ggplot")

  ## localized annotations still parse as plotmath
  expect_s3_class(plot(fit, decimal_mark = ",", cond_word = "se"), "ggplot")
})

test_that("print and summary run without error", {
  fit <- fit_lrp(chickpea_X, chickpea_CV1, step = 0.01)
  expect_output(print(fit), "Linear Response Plateau")
  expect_output(summary(fit), "Model coefficients")

  res <- fit_lrp(chickpea_long, x = "x", cv = "cv", trial = "trial", step = 0.01)
  expect_output(print(res), "3 trials")
})
