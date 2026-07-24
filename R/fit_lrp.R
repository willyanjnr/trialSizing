## ============================================================================
## trialSize :: Linear Response Plateau (LRP) model by grid search
## ============================================================================

## ----------------------------------------------------------------------------
## Weights (internal, shared with fit_qrp)
## ----------------------------------------------------------------------------
## Resolves the `weights` argument to a numeric vector or NULL. The CV of a
## large plot size rests on few plots, so weighting by the number of plots is
## the natural correction; `weights = TRUE` reaches for the `n` column that
## calc_cv_shapes() produces, which is where that count lives.
.resolve_weights <- function(weights, .data = NULL, n_obs = NULL) {
  if (is.null(weights) || identical(weights, FALSE)) return(NULL)

  if (isTRUE(weights)) {
    if (is.null(.data))
      stop("`weights = TRUE` takes the `n` column of the data frame. With the ",
           "vector interface, pass the weights themselves, for example ",
           "`weights = tab$n`.", call. = FALSE)
    if (!"n" %in% names(.data))
      stop("`weights = TRUE` needs a column named `n` (the number of plots per ",
           "shape, as returned by calc_cv_shapes()).\n  Available columns: ",
           paste(names(.data), collapse = ", "), ".", call. = FALSE)
    w <- .data[["n"]]

  } else if (is.character(weights)) {
    if (length(weights) != 1)
      stop("`weights` must be a single column name.", call. = FALSE)
    if (is.null(.data) || !weights %in% names(.data))
      stop(sprintf("Weight column '%s' not found in `.data`.", weights),
           call. = FALSE)
    w <- .data[[weights]]

  } else if (is.numeric(weights)) {
    w <- weights

  } else {
    stop("`weights` must be TRUE/FALSE, a column name, or a numeric vector.",
         call. = FALSE)
  }

  if (!is.numeric(w) || anyNA(w) || any(w <= 0))
    stop("Weights must be numeric, complete and strictly positive.",
         call. = FALSE)
  if (!is.null(n_obs) && length(w) != n_obs)
    stop(sprintf("`weights` has length %d but there are %d observations.",
                 length(w), n_obs), call. = FALSE)
  as.numeric(w)
}

## ----------------------------------------------------------------------------
## Fast SSE profile (internal)
## ----------------------------------------------------------------------------
## Same quantity the grid search needs -- the residual sum of squares at every
## candidate breakpoint -- but computed from cumulative cross-products instead
## of one lm() call per candidate. Both methods are ordinary least squares on a
## basis that changes only through "which points lie below x0", so the sums can
## be accumulated once and reused. This is what makes the bootstrap affordable:
## a replicate costs about a millisecond instead of a second.
##
## Returns a numeric vector aligned with `grid`, with Inf where the candidate
## is infeasible (fewer than two points below it, or no spread in x).
## `w = NULL` is the unweighted fit; otherwise every sum below is a weighted
## sum, which is exactly what weighted least squares needs, and the counts of
## points become sums of weights.
.lrp_ss_profile <- function(x, cv, grid, method = "segment", w = NULL) {
  o  <- order(x)
  xs <- x[o]; ys <- cv[o]
  n  <- length(xs)
  ws <- if (is.null(w)) rep(1, n) else w[o]

  ## cumulative (weighted) sums over the points with x <= x0, indexed 0..n
  z0   <- 0
  cSw  <- c(z0, cumsum(ws))
  cSx  <- c(z0, cumsum(ws * xs));      cSy  <- c(z0, cumsum(ws * ys))
  cSxx <- c(z0, cumsum(ws * xs * xs)); cSxy <- c(z0, cumsum(ws * xs * ys))
  cSyy <- c(z0, cumsum(ws * ys * ys))
  Tw <- cSw[n + 1]; Ty <- cSy[n + 1]; Tyy <- cSyy[n + 1]

  k  <- findInterval(grid, xs)   # number of observations with x <= x0
  ki <- k + 1L                   # index into the padded cumulative sums

  Sw <- cSw[ki]; Sx <- cSx[ki]; Sy <- cSy[ki]
  Sxx <- cSxx[ki]; Sxy <- cSxy[ki]; Syy <- cSyy[ki]
  Sw_r <- Tw - Sw                 # weight above the breakpoint
  Sy_r <- Ty - Sy; Syy_r <- Tyy - Syy

  if (method == "segment") {
    ## line fitted on the points below the breakpoint only
    den <- Sxx - Sx^2 / Sw
    b   <- (Sxy - Sx * Sy / Sw) / den
    a   <- Sy / Sw - b * Sx / Sw
    ss_below <- Syy - a * Sy - b * Sxy
    p  <- a + b * grid
    ss <- ss_below + (Syy_r - 2 * p * Sy_r + Sw_r * p^2)
    ss[k < 2 | !is.finite(den) | den <= 0] <- Inf
  } else {
    ## line fitted on every point, on the basis z = pmin(x, x0)
    Sz  <- Sx + Sw_r * grid
    Szz <- Sxx + Sw_r * grid^2
    Szy <- Sxy + grid * Sy_r
    den <- Szz - Sz^2 / Tw
    b   <- (Szy - Sz * Ty / Tw) / den
    a   <- Ty / Tw - b * Sz / Tw
    ss  <- Tyy - a * Ty - b * Szy
    ss[!is.finite(den) | den <= 0] <- Inf
  }

  ss[!is.finite(ss)] <- Inf
  ## floating-point cancellation can push a zero-residual fit slightly negative
  pmax(ss, 0)
}

## ----------------------------------------------------------------------------
## Bootstrap of the breakpoint (internal)
## ----------------------------------------------------------------------------
## Two questions, one machinery:
##   * how precise is Xo -- resample the shapes (the rows of the CV table) with
##     replacement and refit, giving a percentile interval;
##   * is there a breakpoint at all -- the null "CV falls linearly and never
##     plateaus" leaves the breakpoint unidentified, so the usual likelihood
##     ratio has no chi-square distribution (Davies' problem). Simulate the
##     null instead: resample the residuals of the straight line, refit the
##     plateau model, and see how often the SSE drop matches the observed one.
.lrp_bootstrap <- function(x, cv, grid, method, ss_best, n_boot, conf_level,
                           w = NULL) {
  n <- length(x)

  ## ---- percentile interval for Xo ----
  reps <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    xb  <- x[idx]
    ## a resample needs enough distinct plot sizes to place a breakpoint
    if (length(unique(xb)) < 3) next
    ssb <- .lrp_ss_profile(xb, cv[idx], grid, method, w[idx])
    if (!any(is.finite(ssb))) next
    reps[i] <- grid[which.min(ssb)]
  }
  ok <- !is.na(reps)

  ## ---- bootstrap test for the existence of the breakpoint ----
  w0       <- if (is.null(w)) rep(1, n) else w
  lin      <- if (is.null(w)) stats::lm(cv ~ x) else stats::lm(cv ~ x, weights = w)
  ss_lin   <- sum(w0 * stats::residuals(lin)^2)
  ## relative SSE drop of the plateau model over the straight line
  stat_obs <- (ss_lin - ss_best) / ss_best
  fit_lin  <- stats::fitted(lin)
  ## weights say the residuals have unequal variance, so resample them on the
  ## common scale and put each back on its own; with no weights this is the
  ## plain residual bootstrap
  res_std  <- stats::residuals(lin) * sqrt(w0)

  stat_null <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    y0  <- fit_lin + sample(res_std, n, replace = TRUE) / sqrt(w0)
    ss0 <- .lrp_ss_profile(x, y0, grid, method, w)
    if (!any(is.finite(ss0))) next
    ss0_best <- min(ss0)
    l0 <- if (is.null(w)) stats::lm(y0 ~ x) else stats::lm(y0 ~ x, weights = w)
    ss0_lin <- sum(w0 * stats::residuals(l0)^2)
    stat_null[i] <- (ss0_lin - ss0_best) / ss0_best
  }
  ok_null <- !is.na(stat_null)

  alpha <- (1 - conf_level) / 2
  list(
    ci         = if (sum(ok) >= 2)
      unname(stats::quantile(reps[ok], c(alpha, 1 - alpha))) else c(NA, NA),
    se         = if (sum(ok) >= 2) stats::sd(reps[ok]) else NA_real_,
    replicates = reps[ok],
    n_valid    = sum(ok),
    n_boot     = n_boot,
    conf_level = conf_level,
    ## +1 in numerator and denominator: a p-value of exactly zero is not
    ## attainable with a finite number of replicates
    p_value    = if (any(ok_null))
      (1 + sum(stat_null[ok_null] >= stat_obs)) / (1 + sum(ok_null)) else NA_real_,
    statistic  = stat_obs
  )
}

## ----------------------------------------------------------------------------
## Internal single-series engine (not exported)
## ----------------------------------------------------------------------------
#' Fit one LRP series (internal)
#'
#' Core grid-search fitter for a single \code{x}/\code{cv} series. Users should
#' call [fit_lrp()] instead.
#'
#' @param x,cv numeric vectors of equal length.
#' @param step grid step for the breakpoint search.
#' @param method \code{"segment"} or \code{"ramp"}.
#' @return An object of class \code{"lrp_fit"}.
#' @keywords internal
#' @noRd
.lrp_fit_one <- function(x, cv, step = 0.001, method = "segment",
                         search_range = NULL, start = NULL,
                         local_min_tol = 0.10, bootstrap = FALSE,
                         n_boot = 1000, conf_level = 0.95, w = NULL) {

  if (!is.numeric(x) || !is.numeric(cv))
    stop("`x` and `cv` must be numeric.", call. = FALSE)
  if (length(x) != length(cv))
    stop("`x` and `cv` must have the same length.", call. = FALSE)
  if (anyNA(x) || anyNA(cv))
    stop("`x` and `cv` cannot contain missing values (NA).", call. = FALSE)
  if (length(x) < 4)
    stop("At least 4 observations are required.", call. = FALSE)
  if (length(unique(x)) < 3)
    stop("At least 3 distinct `x` values are required for the breakpoint grid.",
         call. = FALSE)
  if (!is.numeric(local_min_tol) || length(local_min_tol) != 1 ||
      is.na(local_min_tol) || local_min_tol < 0)
    stop("`local_min_tol` must be a single non-negative number.", call. = FALSE)
  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap))
    stop("`bootstrap` must be TRUE or FALSE.", call. = FALSE)
  if (isTRUE(bootstrap)) {
    if (!is.numeric(n_boot) || length(n_boot) != 1 || is.na(n_boot) ||
        n_boot < 2)
      stop("`n_boot` must be a single number of at least 2.", call. = FALSE)
    if (!is.numeric(conf_level) || length(conf_level) != 1 ||
        is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
      stop("`conf_level` must be a single number strictly between 0 and 1.",
           call. = FALSE)
  }

  if (!is.null(w)) {
    if (!is.numeric(w) || length(w) != length(x) || anyNA(w) || any(w <= 0))
      stop("Weights must be numeric, complete, strictly positive and as long ",
           "as `x`.", call. = FALSE)
  }
  wv <- if (is.null(w)) rep(1, length(x)) else w

  x_unique <- sort(unique(x))
  feas_lo  <- x_unique[2]
  feas_hi  <- x_unique[length(x_unique)]

  ## optional restriction of the breakpoint search range
  if (is.null(search_range)) {
    grid_lo <- feas_lo; grid_hi <- feas_hi
  } else {
    if (!is.numeric(search_range) || length(search_range) != 2 ||
        anyNA(search_range) || search_range[1] >= search_range[2])
      stop("`search_range` must be numeric c(lower, upper) with lower < upper.",
           call. = FALSE)
    if (search_range[1] < min(x) || search_range[2] > max(x))
      stop(sprintf("`search_range` must fall within the data range [%.4g, %.4g].",
                   min(x), max(x)), call. = FALSE)
    grid_lo <- max(search_range[1], feas_lo)
    grid_hi <- search_range[2]
    if (grid_lo > grid_hi)
      stop(sprintf(paste0("`search_range` upper bound is below the feasible ",
                          "breakpoint minimum (%.4g); at least 2 points must ",
                          "lie below the breakpoint."), feas_lo), call. = FALSE)
  }

  ## residual sum of squares for a fixed breakpoint
  ss_at <- function(x0) {
    if (method == "segment") {
      m <- x <= x0
      if (sum(m) < 2) return(list(ss = Inf))
      cf <- stats::coef(if (is.null(w)) stats::lm(cv[m] ~ x[m])
                        else stats::lm(cv[m] ~ x[m], weights = w[m]))
    } else {
      z  <- pmin(x, x0)
      cf <- stats::coef(if (is.null(w)) stats::lm(cv ~ z)
                        else stats::lm(cv ~ z, weights = w))
    }
    a <- unname(cf[1]); b <- unname(cf[2]); p <- a + b * x0
    fit <- ifelse(x <= x0, a + b * x, p)
    list(ss = sum(wv * (cv - fit)^2), a = a, b = b, p = p)
  }

  ## grid search over the breakpoint
  grid <- seq(grid_lo, grid_hi, by = step)
  ss_profile <- vapply(grid, function(x0) ss_at(x0)$ss, numeric(1))

  if (!any(is.finite(ss_profile)))
    stop("Grid search failed to find a valid breakpoint.", call. = FALSE)

  ## global optimum
  i_best <- which.min(ss_profile)
  best <- ss_at(grid[i_best])
  a <- best$a; b <- best$b; x0 <- grid[i_best]; p <- best$p
  ss_best <- ss_profile[i_best]

  ## ---- local minima of the SSE profile -------------------------------------
  ## Competing basins matter: with few distinct plot sizes the profile is
  ## stepped, and gradient-based fitters (nls, nlsLM) settle in whichever basin
  ## their starting value falls into. Reporting them shows the user that more
  ## than one breakpoint is defensible.
  local_minima <- NULL
  fin <- is.finite(ss_profile)
  if (sum(fin) > 3) {
    gg <- grid[fin]; pp <- ss_profile[fin]
    dsign <- sign(diff(pp))
    dsign[dsign == 0] <- 1                    # flat steps count as ascending
    idx <- which(diff(dsign) == 2) + 1        # falling then rising
    if (length(idx)) {
      lm_x0 <- gg[idx]; lm_ss <- pp[idx]
      keep <- abs(lm_x0 - x0) > 10 * step     # drop the global one itself
      lm_x0 <- lm_x0[keep]; lm_ss <- lm_ss[keep]
      if (length(lm_x0)) {
        ord <- order(lm_ss)
        excess <- lm_ss[ord] / ss_best - 1
        local_minima <- data.frame(
          breakpoint = lm_x0[ord],
          SSE        = lm_ss[ord],
          SSE_excess = excess,
          competing  = excess <= local_min_tol,
          row.names  = NULL)
      }
    }
  }

  ## ---- optional compatibility fit ------------------------------------------
  ## `start` reproduces what a gradient fitter seeded at that breakpoint would
  ## return: the local minimum of the basin containing `start`.
  compat <- NULL
  if (!is.null(start)) {
    if (!is.numeric(start) || length(start) != 1 || is.na(start))
      stop("`start` must be a single numeric breakpoint value.", call. = FALSE)
    if (start < min(x) || start > max(x))
      stop(sprintf("`start` must fall within the data range [%.4g, %.4g].",
                   min(x), max(x)), call. = FALSE)

    ## walk downhill from the grid point nearest `start`
    i <- which.min(abs(grid - start))
    repeat {
      lo <- if (i > 1) ss_profile[i - 1] else Inf
      hi <- if (i < length(grid)) ss_profile[i + 1] else Inf
      cur <- ss_profile[i]
      if (!is.finite(cur)) { i <- i + 1; next }
      if (is.finite(lo) && lo < cur) { i <- i - 1
      } else if (is.finite(hi) && hi < cur) { i <- i + 1
      } else break
    }
    cs <- ss_at(grid[i])
    compat <- list(start = start, breakpoint = grid[i],
                   coefficients = c(a = cs$a, b = cs$b),
                   plateau = cs$p, SSE = cs$ss,
                   SSE_excess = cs$ss / ss_best - 1)
  }

  ## Fitted values and fit statistics. With weights these follow the lm()
  ## convention: residual sums are weighted and the log-likelihood carries the
  ## 0.5 sum(log w) term. A weighted fit's R2, RMSE, AIC and BIC are therefore
  ## NOT comparable with the unweighted fit of the same data.
  fitted    <- ifelse(x <= x0, a + b * x, p)
  residuals <- cv - fitted
  n   <- length(cv)
  rss <- sum(wv * residuals^2)
  cv_bar    <- sum(wv * cv) / sum(wv)
  r_squared <- 1 - rss / sum(wv * (cv - cv_bar)^2)
  rmse      <- sqrt(rss / sum(wv))

  ## AIC / BIC via the Gaussian log-likelihood (MLE variance = rss / n).
  ## Parameter count k = a, b, Xo, sigma = 4, matching the nls / AIC.default
  ## convention. Only comparable with a QRP fit that uses the same count.
  loglik <- -0.5 * n * (log(2 * pi) + log(rss / n) + 1) + 0.5 * sum(log(wv))
  k <- 4L
  aic <- -2 * loglik + 2 * k
  bic <- -2 * loglik + log(n) * k

  if (b >= 0)
    warning("Estimated slope is non-negative; the decreasing-plateau ",
            "interpretation may not hold for these data.", call. = FALSE)
  if (isTRUE(all.equal(x0, grid_lo)) || isTRUE(all.equal(x0, grid_hi)))
    warning("Breakpoint is at the edge of the search range; the data (or the ",
            "supplied `search_range`) may not bracket the true breakpoint.",
            call. = FALSE)

  ## Competing local minima are reported, never warned about: on stepped SSE
  ## profiles they are the rule rather than the exception, so a warning fires on
  ## nearly every fit and stops carrying information. `local_min_tol` flags them
  ## in $local_minima and in print() instead; the user sets the threshold.

  ## optional uncertainty of the breakpoint
  boot <- NULL
  if (isTRUE(bootstrap))
    boot <- .lrp_bootstrap(x, cv, grid, method, ss_best, n_boot, conf_level, w)

  structure(
    list(
      coefficients = c(a = a, b = b),
      parameters   = c(Breakpoint = x0, Breakpoint_Response = p,
                       R2 = r_squared, RMSE = rmse, AIC = aic, BIC = bic),
      fitted = fitted, residuals = residuals,
      data = data.frame(x = x, cv = cv), method = method, step = step,
      weights = w, search_range = search_range, local_min_tol = local_min_tol,
      local_minima = local_minima, compat = compat, bootstrap = boot,
      sse_profile = data.frame(breakpoint = grid, SSE = ss_profile)
    ),
    class = "lrp_fit"
  )
}

## ----------------------------------------------------------------------------
## Public fitter
## ----------------------------------------------------------------------------
#' Fit the Linear Response Plateau (LRP) model by grid search
#'
#' Fits the linear-plateau (broken-line) model
#' \deqn{f(x) = a + b\,x \ \text{ if } x \le X_0, \qquad
#'       f(x) = a + b\,X_0 \ \text{ if } x > X_0}
#' by profiling the breakpoint \eqn{X_0} over a fine grid. For each candidate
#' breakpoint the linear coefficients are obtained by least squares and the
#' plateau is set to \eqn{a + b\,X_0}; the breakpoint minimizing the residual
#' sum of squares over all observations is returned. Profiling the breakpoint
#' avoids the starting-value sensitivity and local-minima of a direct
#' \code{\link[stats]{nls}} fit, so no initial values are required.
#'
#' The response is typically the coefficient of variation (CV, percent) and the
#' predictor the plot size. Supply the individual CV values (one per basic-unit
#' form), not the means per plot size; repeated \code{x} values are expected.
#'
#' @section Competing breakpoints:
#' With few distinct plot sizes the residual sum of squares is a stepped
#' function of the breakpoint, and it often has several local minima. The grid
#' search always returns the global optimum, but a second basin may fit almost
#' as well, in which case the breakpoint is not sharply identified and different
#' implementations legitimately disagree. Every local minimum of the profile is
#' reported in \code{$local_minima}; those fitting within \code{local_min_tol}
#' of the optimum are flagged in the \code{competing} column and starred by
#' \code{print()}. No warning is issued: on stepped profiles competing basins
#' are common, so a warning would fire on almost every fit. Inspect
#' \code{$local_minima} and \code{$sse_profile} instead, and lower
#' \code{local_min_tol} to flag only near-ties. Gradient-based fitters return
#' whichever basin their starting value lands in; \code{start} reproduces that.
#'
#' @section Weighting by the number of plots:
#' The CV values are not equally reliable. A shape of area \eqn{X} leaves
#' \eqn{n = LC/X} plots in the grid, so the CV of the largest plot size may rest
#' on two plots while the smallest rests on dozens. \code{weights = TRUE} fits
#' by weighted least squares with \eqn{n} as the weight, which is the natural
#' measure of how much information stands behind each point.
#'
#' This is off by default because it is not what the published procedure does,
#' and it moves the answer: the small plot sizes, where the CV is highest, gain
#' most of the weight, so the fitted line is pulled towards them and the
#' breakpoint typically falls. Report both if you use it.
#'
#' Two cautions. Weighting corrects for unequal information, not for
#' dependence: every CV in the table comes from the same grid of basic units,
#' so the points are not independent with or without weights. And a weighted
#' fit's \eqn{R^2}, RMSE, AIC and BIC follow the \code{\link[stats]{lm}}
#' convention of being computed on weighted residuals, so they cannot be
#' compared with those of the unweighted fit.
#'
#' @section Uncertainty of the breakpoint:
#' \code{bootstrap = TRUE} adds two things the point estimate cannot give.
#'
#' The first is a percentile confidence interval: the shapes (the rows of the
#' CV table) are resampled with replacement, the model is refit on each
#' resample, and the empirical quantiles of the resulting breakpoints form the
#' interval. Resamples with fewer than three distinct plot sizes cannot place a
#' breakpoint and are discarded; \code{$bootstrap$n_valid} reports how many were
#' kept. The interval is usually wide, which is the honest reading of a
#' breakpoint estimated from a handful of plot sizes.
#'
#' The second is a test of whether the breakpoint exists at all. Under the null
#' hypothesis that the CV falls linearly and never plateaus, the breakpoint is
#' not identified, so the usual likelihood-ratio statistic does not have a
#' chi-square distribution (Davies, 1987). The null distribution is simulated
#' instead: residuals of the straight-line fit are resampled, the plateau model
#' is refit to each simulated series, and the p-value is the proportion of
#' simulated SSE reductions that match or exceed the observed one. A large
#' p-value means a straight line explains the data as well as the plateau, and
#' the optimal plot size should not be read off this fit.
#'
#' Both are computed on the same breakpoint grid as the main fit, so
#' \code{step} controls their resolution too. Set the random seed before
#' calling to make the result reproducible.
#'
#' @section Two ways to call:
#' \describe{
#'   \item{Vectors}{\code{fit_lrp(x, cv)} fits a single series and returns an
#'     \code{"lrp_fit"}.}
#'   \item{Data frame}{\code{fit_lrp(.data, x = "x", cv = "cv", trial = "trial")}
#'     takes a data frame plus column names. With \code{trial}, one model is fit
#'     per trial and an \code{"lrp_multi"} object (summary table plus the
#'     individual fits) is returned. Column names default to \code{"x"},
#'     \code{"cv"} and \code{"trial"}; missing columns raise a clear error.}
#' }
#'
#' @param .data optional data frame. When supplied, \code{x}, \code{cv} and
#'   \code{trial} are interpreted as column names (character strings). When
#'   \code{NULL} (default) the vector interface is used.
#' @param x,cv either numeric vectors (vector interface) or, when \code{.data}
#'   is a data frame, the names of the predictor and response columns.
#' @param trial optional name of the column identifying the trial. When given,
#'   one model is fit per trial.
#' @param step grid step for the breakpoint search. Default \code{0.001}; larger
#'   values run faster with a slightly coarser breakpoint.
#' @param method how the linear coefficients are estimated at each candidate
#'   breakpoint. \code{"segment"} (default) uses only observations with
#'   \code{x <= X0}, reproducing the Paranaiba et al. (2009) procedure.
#'   \code{"ramp"} uses all observations on the basis \code{pmin(x, X0)}, the
#'   standard least-squares LRP.
#' @param search_range optional numeric \code{c(lower, upper)} restricting the
#'   interval (in units of \code{x}) where the breakpoint is searched. Must fall
#'   within the data range. Useful when the optimum is known to lie in a region
#'   and an outlier could otherwise pull the breakpoint outside it. \code{NULL}
#'   (default) searches the full feasible interval.
#' @param start optional single breakpoint value. The fit is unchanged, but the
#'   result also carries a \code{$compat} element holding the local minimum of
#'   the basin containing \code{start}: the solution a gradient-based fitter
#'   (\code{nls}, \code{nlsLM}) seeded there would return. Use it to reproduce
#'   published results obtained with such implementations, and to see how much
#'   worse they fit.
#' @param local_min_tol relative SSE tolerance (default 0.10) deciding which
#'   local minima count as competing. A second basin fitting within 10\% of the
#'   optimum means the breakpoint is not sharply identified. Only labelling is
#'   affected (the \code{competing} column of \code{$local_minima} and the stars
#'   in \code{print()}); the fit itself never changes, and no warning is issued.
#' @param bootstrap logical; if \code{TRUE}, also estimate the uncertainty of
#'   the breakpoint by resampling (default \code{FALSE}). Off by default because
#'   the published procedure reports the point estimate alone. See the section
#'   "Uncertainty of the breakpoint".
#' @param n_boot number of bootstrap resamples (default 1000). Used only when
#'   \code{bootstrap = TRUE}.
#' @param conf_level confidence level of the percentile interval (default 0.95).
#'   Used only when \code{bootstrap = TRUE}.
#' @param weights weights for a weighted least-squares fit. \code{FALSE}
#'   (default) fits unweighted, as the published procedure does. \code{TRUE}
#'   uses the \code{n} column of \code{.data}, the number of plots behind each
#'   CV, as returned by [calc_cv_shapes()]. A single column name or a numeric
#'   vector are also accepted. See the section "Weighting by the number of
#'   plots".
#'
#' @return
#' For a single series, an object of class \code{"lrp_fit"}: a list with
#' \code{coefficients} (a, b); \code{parameters} (Breakpoint, Breakpoint_Response,
#' R2, RMSE, AIC, BIC); \code{fitted}; \code{residuals}; \code{data};
#' \code{method}; \code{step}; \code{local_min_tol}; \code{local_minima} (the
#' competing basins, with their SSE excess over the optimum and a
#' \code{competing} flag, or \code{NULL}); \code{sse_profile} (the
#' SSE at every candidate breakpoint); \code{compat} when \code{start} was
#' given; and \code{bootstrap} when \code{bootstrap = TRUE}, a list with
#' \code{ci}, \code{se}, \code{p_value} (existence of the breakpoint),
#' \code{statistic}, \code{replicates}, \code{n_valid} and \code{conf_level}.
#' \cr
#' With \code{trial}, an object of class \code{"lrp_multi"}: a list with
#' \code{summary} (one row per trial), \code{fits} (the individual
#' \code{"lrp_fit"} objects) and \code{method}.
#'
#' @references
#' Paranaiba, P. F., Ferreira, D. F. & Morais, A. R. (2009). Tamanho otimo de
#' parcelas experimentais: proposicao de metodos de estimacao.
#' \emph{Revista Brasileira de Biometria}, 27(2), 255-268. \cr
#' Cargnelutti Filho, A. et al. (2025). \emph{Revista Vivencias}, 21(43), 499-513. \cr
#' Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is present
#' only under the alternative. \emph{Biometrika}, 74(1), 33-43. \cr
#' Efron, B. & Tibshirani, R. J. (1993). \emph{An Introduction to the
#' Bootstrap}. Chapman & Hall, New York.
#'
#' @examples
#' ## Chickpea uniformity trial, trial 1 (Cargnelutti Filho et al., 2025).
#' ## One CV per basic-unit form: plot sizes repeat because several shapes
#' ## give the same area.
#' X   <- c(1, 2, 2, 3, 3, 4, 6, 6, 6, 6, 9, 12, 12, 18, 18)
#' CV1 <- c(30.40, 19.51, 23.72, 12.89, 21.32, 16.69, 6.71, 10.75,
#'          17.58, 14.94, 11.93, 3.18, 8.63, 4.25, 11.41)
#'
#' ## A coarse grid runs fast and already reproduces the published Xo to two
#' ## decimals; the default step = 0.001 refines the third.
#' fit <- fit_lrp(X, CV1, step = 0.01)
#' fit
#' coef(fit)
#' fit$parameters[c("Breakpoint", "Breakpoint_Response")]
#'
#' ## CV expected at plot sizes that were not evaluated
#' predict(fit, newx = c(2, 5, 7.5, 15))
#'
#' \donttest{
#' ## Full precision (about ten times slower)
#' fit_lrp(X, CV1)$parameters["Breakpoint"]
#'
#' ## Title and styling belong to plot(), not to the fit
#' plot(fit, title = "Chickpea, trial 1")
#'
#' ## Weighting by the number of plots -----------------------------------------
#' ## The CV of an 18 m2 shape rests on 2 plots, that of a 1 m2 shape on 36.
#' ## Weighting by n pulls the fit towards the small sizes, so Xo falls.
#' n_plots <- c(36, 18, 18, 12, 12, 9, 6, 6, 6, 6, 4, 3, 3, 2, 2)
#' c(unweighted = unname(fit_lrp(X, CV1, step = 0.01)$parameters["Breakpoint"]),
#'   weighted   = unname(fit_lrp(X, CV1, step = 0.01,
#'                               weights = n_plots)$parameters["Breakpoint"]))
#'
#' ## Straight from a grid, `weights = TRUE` finds the n column itself
#' grid1 <- as.matrix(dados_ensaio_C1[dados_ensaio_C1$Rep == 1,
#'                                    paste0("C", 1:8)])
#' tab <- calc_cv_shapes(grid1)
#' fit_lrp(tab, x = "x", cv = "cv", step = 0.05, weights = TRUE)$parameters["Breakpoint"]
#'
#' ## Uncertainty of Xo -------------------------------------------------------
#' ## Off by default, because the published procedure reports the point alone.
#' set.seed(1)
#' unc <- fit_lrp(X, CV1, step = 0.01, bootstrap = TRUE, n_boot = 500)
#' unc
#' unc$bootstrap$ci
#'
#' ## p_value tests whether a breakpoint exists at all: a large value means a
#' ## straight line explains the CV just as well, and no plateau should be read.
#' unc$bootstrap$p_value
#'
#' ## Competing breakpoints ---------------------------------------------------
#' ## Every local minimum of the SSE profile is reported; those fitting within
#' ## local_min_tol of the optimum are flagged as competing.
#' fit$local_minima
#'
#' ## Lower the tolerance to flag only near-ties
#' fit_lrp(X, CV1, step = 0.01, local_min_tol = 0.02)$local_minima
#'
#' ## The whole profile, for inspection
#' plot(fit$sse_profile, type = "l", xlab = "Breakpoint", ylab = "SSE")
#' abline(v = fit$parameters["Breakpoint"], col = "forestgreen")
#'
#' ## What a gradient fitter (nls, nlsLM) seeded at 12 would have returned.
#' ## The reported fit does not change; $compat shows the cost in SSE.
#' fit_lrp(X, CV1, step = 0.01, start = 12)$compat
#'
#' ## Other arguments ---------------------------------------------------------
#' ## Restrict the search when an outlier pulls the breakpoint away
#' fit_lrp(X, CV1, step = 0.01, search_range = c(5, 12))$parameters["Breakpoint"]
#'
#' ## "ramp" estimates the descending line from every observation instead of
#' ## only those below the breakpoint, which can shift the optimum
#' fit_lrp(X, CV1, step = 0.01, method = "ramp")$parameters["Breakpoint"]
#'
#' ## Several trials at once --------------------------------------------------
#' trials <- rbind(
#'   data.frame(x = X, cv = CV1,        trial = "T1"),
#'   data.frame(x = X, cv = CV1 * 0.85, trial = "T2")
#' )
#' res <- fit_lrp(trials, x = "x", cv = "cv", trial = "trial", step = 0.01)
#' res$summary
#'
#' ## CVxo feeds the number of replications
#' calc_replicates(treatments = c(5, 10, 20),
#'                 cv_percent = unname(fit$parameters["Breakpoint_Response"]),
#'                 lsd_percent = c(10, 20))
#' }
#'
#' @seealso [plot.lrp_fit()], [predict.lrp_fit()]
#' @export
fit_lrp <- function(.data = NULL, x = NULL, cv = NULL, trial = NULL,
                    step = 0.001, method = c("segment", "ramp"),
                    search_range = NULL, start = NULL,
                    local_min_tol = 0.10, bootstrap = FALSE,
                    n_boot = 1000, conf_level = 0.95, weights = FALSE) {

  method <- match.arg(method)

  ## Positional vector call: fit_lrp(x_vec, cv_vec)
  if (!is.null(.data) && !is.data.frame(.data)) {
    cv <- x; x <- .data; .data <- NULL
  }

  ## ---- vector interface ----
  if (is.null(.data)) {
    if (is.null(x) || is.null(cv))
      stop("Provide numeric `x` and `cv`, or a data frame as `.data`.",
           call. = FALSE)
    return(.lrp_fit_one(x, cv, step, method, search_range, start,
                        local_min_tol, bootstrap, n_boot, conf_level,
                        .resolve_weights(weights, NULL, length(x))))
  }

  ## ---- data-frame interface ----
  if (!is.data.frame(.data))
    stop("`.data` must be a data frame.", call. = FALSE)

  xcol  <- if (is.null(x))  "x"  else x
  cvcol <- if (is.null(cv)) "cv" else cv
  need  <- c(xcol, cvcol, if (!is.null(trial)) trial)
  missing_cols <- setdiff(need, names(.data))
  if (length(missing_cols))
    stop(sprintf("Column(s) not found in `.data`: %s.\n  Available columns: %s.",
                 paste(missing_cols, collapse = ", "),
                 paste(names(.data), collapse = ", ")), call. = FALSE)

  wvec <- .resolve_weights(weights, .data, nrow(.data))

  if (is.null(trial)) {
    message(sprintf("Using x = '%s', cv = '%s' (single series).", xcol, cvcol))
    return(.lrp_fit_one(.data[[xcol]], .data[[cvcol]], step, method,
                        search_range, start, local_min_tol, bootstrap,
                        n_boot, conf_level, wvec))
  }

  groups <- split(.data, .data[[trial]])
  ## the weights follow their rows into each trial
  wsplit <- if (is.null(wvec)) NULL else split(wvec, .data[[trial]])
  message(sprintf("Using x = '%s', cv = '%s', trial = '%s' -> %d trials.",
                  xcol, cvcol, trial, length(groups)))

  fits <- Map(function(g, nm) .lrp_fit_one(g[[xcol]], g[[cvcol]],
                                           step, method, search_range,
                                           start, local_min_tol,
                                           bootstrap, n_boot, conf_level,
                                           if (is.null(wsplit)) NULL else wsplit[[nm]]),
              groups, names(groups))
  summ <- do.call(rbind, Map(function(f, nm) data.frame(
    trial      = nm,
    a          = round(unname(f$coefficients["a"]), 4),
    b          = round(unname(f$coefficients["b"]), 4),
    breakpoint = round(unname(f$parameters["Breakpoint"]), 4),
    plateau    = round(unname(f$parameters["Breakpoint_Response"]), 4),
    R2         = round(unname(f$parameters["R2"]), 4),
    RMSE       = round(unname(f$parameters["RMSE"]), 4),
    AIC        = round(unname(f$parameters["AIC"]), 3),
    BIC        = round(unname(f$parameters["BIC"]), 3),
    n_local    = if (is.null(f$local_minima)) 0L else nrow(f$local_minima),
    stringsAsFactors = FALSE), fits, names(fits)))
  rownames(summ) <- NULL

  ## the interval only earns its columns when it was actually computed
  if (isTRUE(bootstrap)) {
    summ$Xo_lwr <- round(vapply(fits, function(f) f$bootstrap$ci[1], 0), 4)
    summ$Xo_upr <- round(vapply(fits, function(f) f$bootstrap$ci[2], 0), 4)
    summ$p_breakpoint <- round(vapply(fits, function(f) f$bootstrap$p_value, 0), 4)
  }

  structure(list(fits = fits, summary = summ, method = method),
            class = "lrp_multi")
}

## ----------------------------------------------------------------------------
## Methods
## ----------------------------------------------------------------------------
#' Predictions from an LRP fit
#'
#' @param object an object of class \code{"lrp_fit"}.
#' @param newx numeric vector of predictor values. Defaults to the fitted data.
#' @param ... ignored.
#' @return A numeric vector of predicted responses.
#' @export
predict.lrp_fit <- function(object, newx = NULL, ...) {
  a  <- object$coefficients["a"]; b <- object$coefficients["b"]
  xo <- object$parameters["Breakpoint"]
  if (is.null(newx)) newx <- object$data$x
  unname(ifelse(newx < xo, a + b * newx, a + b * xo))
}

#' Print an LRP fit
#'
#' @param x an object of class \code{"lrp_fit"}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.lrp_fit <- function(x, ...) {
  cat("Linear Response Plateau (LRP) fit\n")
  cat("Method:                 ", x$method,
      if (is.null(x$weights)) "" else " (weighted)", "\n")
  cat("Breakpoint (Xo):        ", sprintf("%.3f", x$parameters["Breakpoint"]), "\n")
  if (!is.null(x$bootstrap)) {
    b <- x$bootstrap
    cat(sprintf("  %.0f%% CI (percentile):  [%.3f, %.3f]   SE %.3f\n",
                100 * b$conf_level, b$ci[1], b$ci[2], b$se))
    cat(sprintf("  breakpoint exists:     p = %.4f  (%d resamples)\n",
                b$p_value, b$n_valid))
  }
  cat("CV at breakpoint:       ", sprintf("%.3f", x$parameters["Breakpoint_Response"]), "\n")
  cat("R2:", sprintf("%.3f", x$parameters["R2"]),
      " RMSE:", sprintf("%.3f", x$parameters["RMSE"]),
      " AIC:", sprintf("%.1f", x$parameters["AIC"]),
      " BIC:", sprintf("%.1f", x$parameters["BIC"]), "\n")

  if (!is.null(x$local_minima)) {
    lmin <- x$local_minima
    tol  <- if (is.null(x$local_min_tol)) 0.10 else x$local_min_tol
    close <- sum(lmin$competing)
    cat("\nLocal minima of the SSE profile (", nrow(lmin), "):\n", sep = "")
    show <- utils::head(lmin, 3)
    for (i in seq_len(nrow(show)))
      cat(sprintf("  Xo = %7.3f   SSE %+6.1f%% vs optimum%s\n",
                  show$breakpoint[i], 100 * show$SSE_excess[i],
                  if (isTRUE(show$competing[i])) "  *" else ""))
    if (nrow(lmin) > 3) cat("  ... see $local_minima for all\n")
    if (close)
      cat(sprintf(paste0("  * fits within %.0f%% of the optimum (local_min_tol);",
                         " breakpoint not sharply identified\n"), 100 * tol))
  }

  if (!is.null(x$compat))
    cat(sprintf("\nCompatibility fit (start = %.3f): Xo = %.3f, SSE %+.1f%%\n",
                x$compat$start, x$compat$breakpoint,
                100 * x$compat$SSE_excess))
  invisible(x)
}

#' Summarize an LRP fit
#'
#' @param object an object of class \code{"lrp_fit"}.
#' @param ... ignored.
#' @return \code{object}, invisibly.
#' @export
summary.lrp_fit <- function(object, ...) {
  cat("Model coefficients:\n"); print(object$coefficients)
  cat("\nGoodness of fit:\n");   print(object$parameters)
  invisible(object)
}

## Save helper (internal) ----------------------------------------------------
.save_lrp <- function(g, file, base, format, dpi, width, height, units,
                      compression) {
  format <- match.arg(format, c("tiff", "png", "jpeg", "pdf", "eps"))
  if (is.null(file)) {
    slug <- gsub("[^A-Za-z0-9]+", "_", base)
    slug <- gsub("^_|_$", "", slug)
    if (!nzchar(slug)) slug <- "lrp_plot"
    file <- paste0(slug, ".", format)
  }
  args <- list(filename = file, plot = g, dpi = dpi,
               width = width, height = height, units = units)
  if (format %in% c("tiff", "tif")) args$compression <- compression
  do.call(ggplot2::ggsave, args)
  message("Saved: ", file)
  invisible(file)
}

## Article-style plot --------------------------------------------------------
#' Plot an LRP fit (publication style)
#'
#' Draws the observed points, the fitted broken line, dotted guides to the
#' breakpoint and the plateau-model annotations, in the layout used in
#' plot-size articles: the model block (equations and \eqn{R^2}) centred at the
#' top and \code{Xo}/\code{CVxo} next to the breakpoint. Axis limits come from
#' the data.
#'
#' @param x an object of class \code{"lrp_fit"}.
#' @param title plot title.
#' @param annotate_model logical; draw the model equations and statistics.
#' @param xlab,ylab axis titles. \code{xlab = NULL} uses "Plot size (m^2)".
#' @param decimal_mark decimal separator for the annotations, "." or ",".
#' @param cond_word conditional word in the equations (e.g. "if" or "se").
#' @param digits_coef,digits_stat decimals for the coefficients (a, b) and for
#'   the statistics (Xo, CVxo, R2).
#' @param point_size,line_size sizes of the points and the fitted line.
#' @param point_colour,line_colour,bp_colour colours of the points, the fitted
#'   line and the breakpoint marker.
#' @param base_size base font size for the theme; axis titles and text scale
#'   from it.
#' @param label_size size of the annotation text (equations, Xo, CVxo).
#' @param title_size size of the plot title. \code{NULL} uses \code{base_size}.
#' @param family font family for the theme and annotations.
#' @param theme optional \pkg{ggplot2} theme object used instead of the default
#'   \code{theme_classic(base_size)}; the axis-line, tick and centred-title
#'   tweaks are applied on top.
#' @param save logical; if \code{TRUE}, write the figure to disk (default FALSE).
#' @param file output file name; if \code{NULL}, derived from \code{title}.
#' @param format one of "tiff", "png", "jpeg", "pdf", "eps". TIFF/PDF/EPS are
#'   preferred for journals; TIFF is written with LZW compression.
#' @param dpi resolution for raster formats.
#' @param width,height,units figure size.
#' @param compression TIFF compression (default "lzw").
#' @param ... ignored.
#' @return A \code{ggplot} object (invisibly when saved).
#' @export
plot.lrp_fit <- function(x, title = "Linear Plateau", annotate_model = TRUE,
                         xlab = NULL, ylab = "CV (%)",
                         decimal_mark = c(".", ","), cond_word = "if",
                         digits_coef = 3, digits_stat = 2,
                         point_size = 2.3, line_size = 0.8,
                         point_colour = "black", line_colour = "black",
                         bp_colour = "red",
                         base_size = 13, label_size = 4.6, title_size = NULL,
                         family = "serif", theme = NULL,
                         save = FALSE, file = NULL,
                         format = c("tiff", "png", "jpeg", "pdf", "eps"),
                         dpi = 300, width = 18, height = 12, units = "cm",
                         compression = "lzw", ...) {

  decimal_mark <- match.arg(decimal_mark)
  format       <- match.arg(format)
  if (is.null(title_size)) title_size <- base_size

  fmtn <- function(v, d) {
    s <- formatC(v, format = "f", digits = d)
    if (decimal_mark == ",") s <- gsub("\\.", ",", s)
    s
  }

  d  <- x$data
  a  <- unname(x$coefficients["a"]); b <- unname(x$coefficients["b"])
  xo <- unname(x$parameters["Breakpoint"])
  p  <- unname(x$parameters["Breakpoint_Response"])
  r2 <- unname(x$parameters["R2"])

  xmax <- max(d$x) * 1.05
  ymax <- max(d$cv) * 1.12
  curve <- data.frame(x = seq(0, max(d$x), length.out = 400))
  curve$cv <- predict(x, curve$x)

  if (is.null(xlab)) xlab <- expression("Plot size (" * m^2 * ")")

  g <- ggplot2::ggplot(d, ggplot2::aes(x, cv)) +
    ggplot2::annotate("segment", x = 0, xend = xo, y = p, yend = p,
                      linetype = 3, linewidth = 0.5) +
    ggplot2::annotate("segment", x = xo, xend = xo, y = 0, yend = p,
                      linetype = 3, linewidth = 0.5) +
    ggplot2::geom_line(data = curve, ggplot2::aes(x, cv),
                       linewidth = line_size, colour = line_colour) +
    ggplot2::geom_point(size = point_size, colour = point_colour) +
    ggplot2::annotate("point", x = xo, y = p, colour = bp_colour,
                      size = point_size * 1.3)

  if (annotate_model) {
    sgn  <- if (b < 0) "-" else "+"
    a_s  <- fmtn(a, digits_coef);  b_s <- fmtn(abs(b), digits_coef)
    xo_s <- fmtn(xo, digits_stat); p_s <- fmtn(p, digits_stat)
    r2_s <- fmtn(r2, digits_stat)

    eq1 <- paste0("CV[(x)]=='", a_s, "'", sgn, "'", b_s, "'*X~~'",
                  cond_word, "'~~X<='", xo_s, "'")
    eq2 <- paste0("CV[(x)]=='", p_s, "'~~'", cond_word, "'~~X>'", xo_s, "'")
    eq3 <- paste0("R^2=='", r2_s, "'")
    lXo <- paste0("X[o]=='", xo_s, "'")
    lCV <- paste0("CV[Xo]=='", p_s, "'")

    ## model block: centred, top
    xm <- 0.55 * xmax
    ym <- ymax * c(0.97, 0.89, 0.81)
    for (i in seq_len(3))
      g <- g + ggplot2::annotate("text", x = xm, y = ym[i],
                                 label = c(eq1, eq2, eq3)[i], parse = TRUE,
                                 hjust = 0.5, size = label_size, family = family)

    ## breakpoint block: right of the vertical guide, below the plateau
    xb <- xo + 0.03 * max(d$x)
    yb <- c(max(p - 0.11 * ymax, 0.10 * ymax),
            max(p - 0.20 * ymax, 0.02 * ymax))
    for (i in 1:2)
      g <- g + ggplot2::annotate("text", x = xb, y = yb[i],
                                 label = c(lXo, lCV)[i], parse = TRUE,
                                 hjust = 0, size = label_size, family = family)
  }

  g <- g +
    ggplot2::scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    (if (is.null(theme))
      ggplot2::theme_classic(base_size = base_size, base_family = family)
     else theme) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = title_size),
      axis.line  = ggplot2::element_line(colour = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(colour = "black"),
      axis.text  = ggplot2::element_text(colour = "black")
    )

  if (save) {
    .save_lrp(g, file, title, format, dpi, width, height, units, compression)
    return(invisible(g))
  }
  g
}

#' Print LRP fits for several trials
#'
#' @param x an object of class \code{"lrp_multi"}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.lrp_multi <- function(x, ...) {
  cat(sprintf("LRP fits for %d trials (method = %s)\n\n",
              length(x$fits), x$method))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Summarize LRP fits for several trials
#'
#' @param object an object of class \code{"lrp_multi"}.
#' @param ... ignored.
#' @return \code{object}, invisibly.
#' @export
summary.lrp_multi <- function(object, ...) {
  print(object$summary, row.names = FALSE)
  invisible(object)
}

#' Plot LRP fits for several trials (paginated grid)
#'
#' Arranges the per-trial plots in a grid of at most 6 panels (3 rows by 2
#' columns). With more than 6 trials the panels are paginated: each page is a
#' separate 3x2 figure, and with \code{save = TRUE} each page is written to its
#' own file. Requires the \pkg{patchwork} package.
#'
#' @param x an object of class \code{"lrp_multi"}.
#' @param save logical; write the page(s) to disk (default FALSE).
#' @param file base file name; page number and extension are appended when there
#'   is more than one page.
#' @param format,dpi,width,height,units,compression passed to the saver.
#' @param ... styling arguments forwarded to [plot.lrp_fit()]
#'   (for example \code{decimal_mark}, \code{cond_word}, \code{label_size}).
#' @return A list of page objects, invisibly.
#' @export
plot.lrp_multi <- function(x, save = FALSE, file = NULL,
                           format = c("tiff", "png", "jpeg", "pdf", "eps"),
                           dpi = 300, width = 18, height = 24, units = "cm",
                           compression = "lzw", ...) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required to arrange multiple trials. ",
         "Install it with install.packages('patchwork').", call. = FALSE)
  format <- match.arg(format)

  plots <- Map(function(f, nm) plot(f, title = as.character(nm), save = FALSE, ...),
               x$fits, names(x$fits))

  idx    <- seq_along(plots)
  chunks <- split(idx, ceiling(idx / 6))
  pages  <- lapply(chunks, function(ix)
    patchwork::wrap_plots(plots[ix], nrow = 3, ncol = 2))

  if (save) {
    base <- if (is.null(file)) "lrp_trials" else sub("\\.[^.]*$", "", file)
    n <- length(pages)
    for (i in seq_len(n)) {
      fn <- if (n == 1) paste0(base, ".", format)
      else sprintf("%s_p%d.%s", base, i, format)
      .save_lrp(pages[[i]], fn, base, format, dpi, width, height, units,
                compression)
    }
  } else {
    for (pg in pages) print(pg)
  }
  invisible(pages)
}
