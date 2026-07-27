## Regression tests for fit_qrp()
## Reference: Cargnelutti Filho et al. (2025), Revista Vivencias 21(43),
## 499-513, Table 2 (chickpea, QRP method).

test_that("fit_qrp reproduces the published chickpea optima", {
  pub <- chickpea_published$QRP

  ## the article rounds intermediate quantities, so some cells land one
  ## hundredth away; a tolerance of 0.011 accepts that while still catching
  ## any real drift
  for (i in seq_along(chickpea_cv_list)) {
    fit <- fit_qrp(chickpea_X, chickpea_cv_list[[i]], step = 0.01)
    expect_equal(unname(fit$parameters["Breakpoint"]), pub$Xo[i],
                 tolerance = 0.011, info = paste("trial", i))
    expect_equal(unname(fit$parameters["Breakpoint_Response"]), pub$CVxo[i],
                 tolerance = 0.011, info = paste("trial", i))
  }
})

test_that("the mean optimum matches the article's overall value", {
  xo <- vapply(chickpea_cv_list,
               function(y) unname(fit_qrp(chickpea_X, y, step = 0.01)$parameters["Breakpoint"]),
               numeric(1))
  expect_equal(round(mean(xo), 2), chickpea_published$QRP$mean_Xo)
})

test_that("estimates are stable to three decimals", {
  ## the default step = 0.001 is what pins the third decimal; a coarser
  ## grid cannot resolve it
  fit <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.001)
  expect_equal(unname(fit$parameters["Breakpoint"]), 10.965, tolerance = 1e-6)
  expect_equal(unname(fit$parameters["Breakpoint_Response"]), 7.500,
               tolerance = 1e-3)
  expect_equal(unname(fit$parameters["R2"]), 0.741, tolerance = 1e-3)
})

test_that("the object has the documented structure", {
  fit <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.01)

  expect_s3_class(fit, "qrp_fit")
  expect_named(fit$coefficients, c("a", "b", "c"))
  expect_true(all(c("Breakpoint", "Breakpoint_Response", "R2", "R2_adj",
                    "RMSE", "MAE", "AIC", "BIC", "SSE", "MSE") %in%
                  names(fit$parameters)))
})

test_that("the breakpoint equals the vertex of the parabola", {
  fit <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.01)
  cf  <- fit$coefficients

  ## Xo = -b / (2c)
  expect_equal(unname(-cf["b"] / (2 * cf["c"])),
               unname(fit$parameters["Breakpoint"]))

  ## plateau = a - b^2 / (4c)
  expect_equal(unname(cf["a"] - cf["b"]^2 / (4 * cf["c"])),
               unname(fit$parameters["Breakpoint_Response"]),
               tolerance = 1e-8)
})

test_that("the descending arm is convex and the plateau is flat", {
  fit <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.01)
  xo  <- unname(fit$parameters["Breakpoint"])
  plateau <- unname(fit$parameters["Breakpoint_Response"])

  expect_lt(unname(fit$coefficients["b"]), 0)   # descends
  expect_gt(unname(fit$coefficients["c"]), 0)   # opens upward
  expect_equal(predict(fit, xo + 5), plateau)
  expect_gt(predict(fit, 1), predict(fit, 5))
})

test_that("fitted values and residuals are internally consistent", {
  fit <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.01)

  expect_equal(fit$fitted + fit$residuals, chickpea_CV1)
  expect_equal(unname(fit$parameters["SSE"]), sum(fit$residuals^2))
  expect_equal(unname(fit$parameters["RMSE"]),
               sqrt(unname(fit$parameters["MSE"])))
})

test_that("the grid search matches a converged nls fit", {
  ## same model, fitted by nls starting from the grid-search solution:
  ## if the grid found the least-squares optimum, the SSE must agree
  fit <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.01)
  cf  <- fit$coefficients

  qrp_model <- function(x, a, b, c) {
    bp <- -b / (2 * c)
    ifelse(x <= bp, a + b * x + c * x^2, a - b^2 / (4 * c))
  }

  nls_fit <- try(
    stats::nls(cv ~ qrp_model(x, a, b, c),
               data = data.frame(x = chickpea_X, cv = chickpea_CV1),
               start = list(a = unname(cf["a"]), b = unname(cf["b"]),
                            c = unname(cf["c"]))),
    silent = TRUE)

  skip_if(inherits(nls_fit, "try-error"), "nls did not converge")

  expect_equal(unname(fit$parameters["SSE"]),
               sum(stats::residuals(nls_fit)^2),
               tolerance = 1e-4)
})

test_that("QRP gives a larger optimum than LRP on the same data", {
  for (y in chickpea_cv_list) {
    qrp <- unname(fit_qrp(chickpea_X, y, step = 0.01)$parameters["Breakpoint"])
    lrp <- unname(fit_lrp(chickpea_X, y, step = 0.01)$parameters["Breakpoint"])
    expect_gt(qrp, lrp)
  }
})

test_that("the data-frame interface fits one model per trial", {
  res <- fit_qrp(chickpea_long, x = "x", cv = "cv", trial = "trial", step = 0.01)

  expect_s3_class(res, "qrp_multi")
  expect_equal(nrow(res$summary), 3)
  expect_true("c" %in% names(res$summary))
  expect_equal(res$summary$breakpoint, chickpea_published$QRP$Xo,
               tolerance = 0.011)
})

test_that("search_range and step behave as documented", {
  ## a range containing the optimum must leave it untouched, so this one is
  ## compared against the full-precision value and needs step = 0.001
  inside <- fit_qrp(chickpea_X, chickpea_CV1, search_range = c(6, 15),
                    step = 0.001)
  expect_equal(unname(inside$parameters["Breakpoint"]), 10.965,
               tolerance = 1e-3)

  coarse <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.01)
  expect_equal(round(unname(coarse$parameters["Breakpoint"]), 1), 11.0)

  expect_error(fit_qrp(chickpea_X, chickpea_CV1, search_range = c(20, 40), step = 0.01),
               "within the data range")
})

test_that("invalid input is rejected", {
  expect_error(fit_qrp(chickpea_X, chickpea_CV1[-1]), "same length")
  expect_error(fit_qrp(c(chickpea_X[-1], NA), chickpea_CV1), "missing values")
  expect_error(fit_qrp(1:3, c(1, 2, 3)), "At least 4 observations")
  expect_error(fit_qrp(chickpea_long, x = "x", cv = "CV"), "not found")
})

test_that("plot, print and summary work", {
  skip_if_not_installed("ggplot2")

  fit <- fit_qrp(chickpea_X, chickpea_CV1, step = 0.01)
  expect_s3_class(plot(fit, title = "Trial 1"), "ggplot")
  expect_s3_class(plot(fit, decimal_mark = ",", cond_word = "se"), "ggplot")
  expect_output(print(fit), "Quadratic Response Plateau")
  expect_output(summary(fit), "Model coefficients")
})
