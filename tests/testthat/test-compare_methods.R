## Regression tests for compare_methods().
## The recommendations should reproduce the per-trial optima published in
## Cargnelutti Filho et al. (2025), and preserve the ordering MCM < LRP < QRP.

test_that("the object has the documented structure", {
  cmp <- compare_methods(data.frame(x = chickpea_X, cv = chickpea_CV1),
                         step = 0.05)

  expect_s3_class(cmp, "method_comparison")
  expect_equal(names(cmp$summary),
               c("trial", "method", "Xo", "CVxo", "R2", "RMSE"))
  expect_equal(cmp$summary$method, c("MCM", "LRP", "QRP"))
  expect_false(cmp$meta$from_grid)
  expect_false(cmp$meta$weighted)
  expect_equal(cmp$meta$n_trials, 1)
})

test_that("a CV table reproduces the published per-trial optima", {
  cmp <- compare_methods(data.frame(x = chickpea_X, cv = chickpea_CV1),
                         step = 0.05)
  xo <- setNames(cmp$summary$Xo, cmp$summary$method)

  ## Table 2 of the article, trial 1
  expect_equal(unname(xo["MCM"]), chickpea_published$MCM$Xo[1], tolerance = 0.01)
  expect_equal(unname(xo["LRP"]), chickpea_published$LRP$Xo[1], tolerance = 0.01)
  expect_equal(unname(xo["QRP"]), chickpea_published$QRP$Xo[1], tolerance = 0.01)
})

test_that("the methods keep the ordering MCM < LRP < QRP", {
  cmp <- compare_methods(data.frame(x = chickpea_X, cv = chickpea_CV1),
                         step = 0.05)
  xo <- setNames(cmp$summary$Xo, cmp$summary$method)

  expect_lt(xo["MCM"], xo["LRP"])
  expect_lt(xo["LRP"], xo["QRP"])
})

test_that("R2 and RMSE are recomputed for the fitted methods", {
  cmp <- compare_methods(data.frame(x = chickpea_X, cv = chickpea_CV1),
                         step = 0.05)
  fitted <- cmp$summary[cmp$summary$method %in% c("MCM", "LRP", "QRP"), ]

  expect_true(all(is.finite(fitted$R2)))
  expect_true(all(is.finite(fitted$RMSE)))
  expect_true(all(fitted$R2 > 0 & fitted$R2 < 1))
})

test_that("a raw grid adds the Paranaiba method", {
  cmp <- compare_methods(chickpea_grid1, step = 0.05)

  expect_true(cmp$meta$from_grid)
  expect_equal(cmp$summary$method, c("MCM", "LRP", "QRP", "Paranaiba"))

  ## the closed form has no residuals, so no R2 or RMSE
  para <- cmp$summary[cmp$summary$method == "Paranaiba", ]
  expect_true(is.na(para$R2))
  expect_true(is.na(para$RMSE))
  expect_equal(para$Xo, calc_paranaiba(chickpea_grid1)$summary$Xo,
               tolerance = 1e-4)
})

test_that("bootstrap adds interval and p-value columns", {
  no_boot <- compare_methods(chickpea_grid1, step = 0.05)
  expect_false(any(c("Xo_lwr", "Xo_upr", "p_breakpoint") %in%
                     names(no_boot$summary)))

  set.seed(1)
  boot <- compare_methods(chickpea_grid1, step = 0.05,
                          bootstrap = TRUE, n_boot = 50)
  expect_true(all(c("Xo_lwr", "Xo_upr", "p_breakpoint") %in%
                    names(boot$summary)))
  expect_true(boot$meta$bootstrap)
})

test_that("weighting is recorded and passed through", {
  cmp <- compare_methods(chickpea_grid1, step = 0.05, weights = TRUE)
  expect_true(cmp$meta$weighted)
})

test_that("the methods argument selects a subset", {
  cmp <- compare_methods(chickpea_grid1, step = 0.05,
                         methods = c("lrp", "qrp"))
  expect_equal(cmp$summary$method, c("LRP", "QRP"))
})

test_that("a list of grids compares within each trial", {
  cmp <- compare_methods(list(A = chickpea_grid1, B = chickpea_grid2),
                         step = 0.05)

  expect_equal(sort(unique(cmp$summary$trial)), c("A", "B"))
  expect_equal(nrow(cmp$summary), 8)   # 4 methods x 2 trials
})

test_that("invalid input is rejected", {
  ## Paranaiba needs the raw grid, not a CV table
  expect_error(
    compare_methods(data.frame(x = chickpea_X, cv = chickpea_CV1),
                    methods = "paranaiba"),
    "raw grid")
  expect_error(
    compare_methods(chickpea_grid1, methods = "banana"),
    "Unknown method")
  expect_error(compare_methods(1:10), "data frame")
  expect_error(
    compare_methods(data.frame(a = 1, b = 2), x = "x", cv = "cv"),
    "not found")
})

test_that("print and summary work", {
  cmp <- compare_methods(chickpea_grid1, step = 0.05)
  expect_output(print(cmp), "methods compared")
  expect_output(print(cmp), "Xo ranges from")
  expect_output(summary(cmp), "MCM")
})
