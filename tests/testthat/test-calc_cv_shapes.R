## Regression tests for calc_cv_shapes()
## The CV table it builds from the raw grid must reproduce Table 1 of
## Cargnelutti Filho et al. (2025), Revista Vivencias 21(43), 499-513, from
## which chickpea_CV1 is taken.

test_that("the CV table reproduces the published Table 1", {
  tab <- calc_cv_shapes(chickpea_grid1)

  ## the shapes a 6 x 6 grid admits, and the plots each yields
  expect_equal(nrow(tab), 15)
  expect_equal(tab$x, chickpea_X)
  expect_equal(tab$n, chickpea_n)

  ## the headline: the CVs match the article to two decimals
  expect_equal(round(tab$cv, 2), chickpea_CV1, tolerance = 0.011)
})

test_that("the object has the documented structure", {
  tab <- calc_cv_shapes(chickpea_grid1)

  expect_s3_class(tab, "data.frame")
  expect_equal(names(tab),
               c("trial", "X_L", "X_C", "x", "n", "mean", "sd", "cv"))
  expect_equal(tab$x, tab$X_L * tab$X_C)
  expect_equal(tab$n, 36 / tab$x)
  expect_true(all(tab$cv == 100 * tab$sd / tab$mean))
  expect_true(all(tab$trial == "Trial 1"))
})

test_that("shapes of equal area but different orientation are both kept", {
  tab <- calc_cv_shapes(chickpea_grid1)
  two <- tab[tab$x == 2, ]

  ## 1 x 2 and 2 x 1 cover the same units along different directions
  expect_equal(nrow(two), 2)
  expect_setequal(paste(two$X_L, two$X_C), c("1 2", "2 1"))
  expect_false(isTRUE(all.equal(two$cv[1], two$cv[2])))
})

test_that("the table feeds the fitters unchanged", {
  tab <- calc_cv_shapes(chickpea_grid1)
  fit <- fit_lrp(tab, x = "x", cv = "cv", step = 0.05)

  expect_s3_class(fit, "lrp_fit")
  expect_equal(nrow(fit$data), 15)
})

test_that("a list of grids gives one block of rows per trial", {
  tab <- calc_cv_shapes(list(A = chickpea_grid1, B = chickpea_grid2))

  expect_equal(sort(unique(tab$trial)), c("A", "B"))
  expect_equal(as.integer(table(tab$trial)), c(15L, 15L))

  ## the first trial matches the same grid built alone
  alone <- calc_cv_shapes(chickpea_grid1)
  expect_equal(tab$cv[tab$trial == "A"], alone$cv)
})

test_that("the long interface rebuilds the grid from the indices", {
  long <- expand.grid(col = 1:6, row = 1:6)
  long$mf <- as.vector(t(chickpea_grid1))

  from_long   <- calc_cv_shapes(long, value = "mf", row_id = "row",
                                col_id = "col")
  from_matrix <- calc_cv_shapes(chickpea_grid1)
  expect_equal(from_long$cv, from_matrix$cv)

  ## the grid comes from the indices, not the row order of the file
  set.seed(1)
  shuffled <- long[sample(nrow(long)), ]
  expect_equal(
    calc_cv_shapes(shuffled, value = "mf", row_id = "row", col_id = "col")$cv,
    from_matrix$cv)
})

test_that("min_plots controls which shapes are kept", {
  ## the whole grid (one plot, no variance) is dropped by default
  tab <- calc_cv_shapes(chickpea_grid1)
  expect_false(36 %in% tab$x)
  expect_true(all(tab$n >= 2))

  ## a higher floor drops more of the large, sparsely-estimated shapes
  fewer <- calc_cv_shapes(chickpea_grid1, min_plots = 6)
  expect_true(all(fewer$n >= 6))
  expect_lt(nrow(fewer), nrow(tab))
})

test_that("invalid input is rejected", {
  expect_error(calc_cv_shapes(chickpea_grid1, min_plots = 1),
               "at least 2")
  expect_error(calc_cv_shapes("not a grid"),
               "matrix, a list of matrices")
  expect_error(calc_cv_shapes(data.frame(a = 1:4)), "give `value`")

  na_grid <- chickpea_grid1; na_grid[1, 1] <- NA
  expect_error(calc_cv_shapes(na_grid), "missing values")
})
