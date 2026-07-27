## ============================================================================
## trialSizing :: Quadratic Response Plateau (QRP) model by grid search
## ----------------------------------------------------------------------------
## Mirrors fit_lrp(): same interfaces, S3 classes, publication plot and saving.
## The quadratic-plateau joins smoothly (vertex at the breakpoint), so fixing
## the breakpoint X0 makes the model linear in (A, c):
##   x <= X0:  cv = A + c (x - X0)^2      x > X0:  cv = A   (plateau)
## The breakpoint is profiled over a grid, avoiding nls starting values.
## Relation to the a + b x + c x^2 form: a = A + c X0^2, b = -2 c X0.
## Uses the shared internal saver .save_lrp() defined in fit_lrp.R.
## ============================================================================

## Fast SSE profile (internal) -----------------------------------------------
## Counterpart of .lrp_ss_profile() for the quadratic-plateau basis. Fixing X0
## makes the model an ordinary regression of cv on z = (x - X0)^2 truncated at
## zero, so expanding the square turns every cross-product into a polynomial in
## X0 whose coefficients are cumulative sums over the points below X0:
##   Sz  = Sxx - 2 X0 Sx + k X0^2
##   Szz = Sx4 - 4 X0 Sx3 + 6 X0^2 Sxx - 4 X0^3 Sx + k X0^4
##   Szy = Sxxy - 2 X0 Sxy + X0^2 Sy
## One pass over the data then gives the whole profile, which is what makes the
## bootstrap affordable. Inf marks candidates with no usable spread in z.
.qrp_ss_profile <- function(x, cv, grid, w = NULL) {
  o  <- order(x)
  xs <- x[o]; ys <- cv[o]
  n  <- length(xs)
  ws <- if (is.null(w)) rep(1, n) else w[o]

  z0    <- 0
  cSw   <- c(z0, cumsum(ws))
  cSx   <- c(z0, cumsum(ws * xs));        cSy   <- c(z0, cumsum(ws * ys))
  cSxx  <- c(z0, cumsum(ws * xs^2));      cSxy  <- c(z0, cumsum(ws * xs * ys))
  cSx3  <- c(z0, cumsum(ws * xs^3));      cSx4  <- c(z0, cumsum(ws * xs^4))
  cSxxy <- c(z0, cumsum(ws * xs^2 * ys))
  Tw  <- cSw[n + 1]; Ty <- cSy[n + 1]; Tyy <- sum(ws * ys^2)

  k  <- findInterval(grid, xs)
  ki <- k + 1L
  Sw <- cSw[ki]
  Sx <- cSx[ki]; Sy <- cSy[ki]; Sxx <- cSxx[ki]; Sxy <- cSxy[ki]
  Sx3 <- cSx3[ki]; Sx4 <- cSx4[ki]; Sxxy <- cSxxy[ki]

  g <- grid
  Sz  <- Sxx - 2 * g * Sx + Sw * g^2
  Szz <- Sx4 - 4 * g * Sx3 + 6 * g^2 * Sxx - 4 * g^3 * Sx + Sw * g^4
  Szy <- Sxxy - 2 * g * Sxy + g^2 * Sy

  den <- Szz - Sz^2 / Tw
  cc  <- (Szy - Sz * Ty / Tw) / den
  A   <- Ty / Tw - cc * Sz / Tw
  ss  <- Tyy - A * Ty - cc * Szy

  ss[!is.finite(den) | den <= 0 | !is.finite(ss)] <- Inf
  pmax(ss, 0)
}

## Bootstrap of the breakpoint (internal) -------------------------------------
## Same two questions as .lrp_bootstrap(), but the null of the existence test
## differs. A quadratic-plateau that never plateaus is simply a quadratic, so
## the null here is the unconstrained quadratic fitted on every point, not the
## straight line used for the LRP.
.qrp_bootstrap <- function(x, cv, grid, ss_best, n_boot, conf_level, w = NULL) {
  n  <- length(x)
  w0 <- if (is.null(w)) rep(1, n) else w

  reps <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    xb  <- x[idx]
    if (length(unique(xb)) < 3) next
    ssb <- .qrp_ss_profile(xb, cv[idx], grid, w[idx])
    if (!any(is.finite(ssb))) next
    reps[i] <- grid[which.min(ssb)]
  }
  ok <- !is.na(reps)

  quad <- if (is.null(w)) stats::lm(cv ~ x + I(x^2))
  else stats::lm(cv ~ x + I(x^2), weights = w)
  ss_quad  <- sum(w0 * stats::residuals(quad)^2)
  stat_obs <- (ss_quad - ss_best) / ss_best
  fit_q    <- stats::fitted(quad)
  res_std  <- stats::residuals(quad) * sqrt(w0)

  stat_null <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    y0  <- fit_q + sample(res_std, n, replace = TRUE) / sqrt(w0)
    ss0 <- .qrp_ss_profile(x, y0, grid, w)
    if (!any(is.finite(ss0))) next
    q0 <- if (is.null(w)) stats::lm(y0 ~ x + I(x^2))
    else stats::lm(y0 ~ x + I(x^2), weights = w)
    ss0_quad <- sum(w0 * stats::residuals(q0)^2)
    stat_null[i] <- (ss0_quad - min(ss0)) / min(ss0)
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
    p_value    = if (any(ok_null))
      (1 + sum(stat_null[ok_null] >= stat_obs)) / (1 + sum(ok_null)) else NA_real_,
    statistic  = stat_obs,
    null_model = "quadratic"
  )
}

## Internal single-series engine ---------------------------------------------
.qrp_fit_one <- function(x, cv, step = 0.001, search_range = NULL,
                         start = NULL, local_min_tol = 0.10,
                         bootstrap = FALSE, n_boot = 1000, conf_level = 0.95,
                         w = NULL) {

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
                          "breakpoint minimum (%.4g)."), feas_lo), call. = FALSE)
  }

  ## residual sum of squares for a fixed breakpoint (linear in A and c)
  ss_at <- function(x0) {
    z  <- ifelse(x <= x0, (x - x0)^2, 0)
    cf <- stats::coef(if (is.null(w)) stats::lm(cv ~ z)
                      else stats::lm(cv ~ z, weights = w))
    A  <- unname(cf[1]); cc <- unname(cf[2])
    fit <- A + ifelse(x <= x0, cc * (x - x0)^2, 0)
    list(ss = sum(wv * (cv - fit)^2), A = A, cc = cc)
  }

  grid <- seq(grid_lo, grid_hi, by = step)
  ss_profile <- vapply(grid, function(x0) ss_at(x0)$ss, numeric(1))

  if (!any(is.finite(ss_profile)))
    stop("Grid search failed to find a valid breakpoint.", call. = FALSE)

  i_best <- which.min(ss_profile)
  best <- ss_at(grid[i_best])
  x0 <- grid[i_best]; A <- best$A; cc <- best$cc
  ss_best <- ss_profile[i_best]

  ## ---- local minima of the SSE profile -------------------------------------
  ## The quadratic-plateau joins smoothly, so this profile is usually a single
  ## smooth basin (unlike the LRP, whose kink makes it stepped). Competing
  ## minima are still reported when they occur.
  local_minima <- NULL
  fin <- is.finite(ss_profile)
  if (sum(fin) > 3) {
    gg <- grid[fin]; pp <- ss_profile[fin]
    dsign <- sign(diff(pp)); dsign[dsign == 0] <- 1
    idx <- which(diff(dsign) == 2) + 1
    if (length(idx)) {
      lm_x0 <- gg[idx]; lm_ss <- pp[idx]
      keep <- abs(lm_x0 - x0) > 10 * step
      lm_x0 <- lm_x0[keep]; lm_ss <- lm_ss[keep]
      if (length(lm_x0)) {
        ord <- order(lm_ss)
        excess <- lm_ss[ord] / ss_best - 1
        local_minima <- data.frame(
          breakpoint = lm_x0[ord], SSE = lm_ss[ord],
          SSE_excess = excess, competing = excess <= local_min_tol,
          row.names = NULL)
      }
    }
  }

  ## ---- optional compatibility fit ------------------------------------------
  compat <- NULL
  if (!is.null(start)) {
    if (!is.numeric(start) || length(start) != 1 || is.na(start))
      stop("`start` must be a single numeric breakpoint value.", call. = FALSE)
    if (start < min(x) || start > max(x))
      stop(sprintf("`start` must fall within the data range [%.4g, %.4g].",
                   min(x), max(x)), call. = FALSE)
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
                   plateau = cs$A, SSE = cs$ss,
                   SSE_excess = cs$ss / ss_best - 1)
  }
  a <- A + cc * x0^2; b <- -2 * cc * x0; c_par <- cc
  breakpoint <- x0; plateau <- A

  fitted    <- A + ifelse(x <= x0, cc * (x - x0)^2, 0)
  residuals <- cv - fitted
  n   <- length(cv)
  p   <- 3L
  ## with weights these follow the lm() convention; see fit_lrp
  sse <- sum(wv * residuals^2)
  mse <- sse / sum(wv)
  rmse <- sqrt(mse)
  mae  <- sum(wv * abs(residuals)) / sum(wv)
  cv_bar  <- sum(wv * cv) / sum(wv)
  r2      <- 1 - sse / sum(wv * (cv - cv_bar)^2)
  r2_adj  <- 1 - (1 - r2) * (n - 1) / (n - p)
  loglik  <- -0.5 * n * (log(2 * pi) + log(sse / n) + 1) + 0.5 * sum(log(wv))
  k   <- 4L
  aic <- -2 * loglik + 2 * k
  bic <- -2 * loglik + log(n) * k

  if (cc <= 0)
    warning("Estimated quadratic coefficient is non-positive; the ",
            "decreasing-then-plateau shape may not hold for these data.",
            call. = FALSE)
  if (isTRUE(all.equal(x0, grid_lo)) || isTRUE(all.equal(x0, grid_hi)))
    warning("Breakpoint is at the edge of the search range; the data (or the ",
            "supplied `search_range`) may not bracket the true breakpoint.",
            call. = FALSE)

  ## Competing local minima are reported, never warned about; see fit_lrp.
  ## `local_min_tol` only flags them in $local_minima and in print().

  boot <- NULL
  if (isTRUE(bootstrap))
    boot <- .qrp_bootstrap(x, cv, grid, ss_best, n_boot, conf_level, w)

  structure(
    list(
      coefficients = c(a = a, b = b, c = c_par),
      parameters   = c(Breakpoint = breakpoint, Breakpoint_Response = plateau,
                       R2 = r2, R2_adj = r2_adj, RMSE = rmse, MAE = mae,
                       AIC = aic, BIC = bic, SSE = sse, MSE = mse),
      fitted = fitted, residuals = residuals,
      data = data.frame(x = x, cv = cv), step = step, weights = w,
      search_range = search_range, local_min_tol = local_min_tol,
      local_minima = local_minima, compat = compat, bootstrap = boot,
      sse_profile = data.frame(breakpoint = grid, SSE = ss_profile)
    ),
    class = "qrp_fit"
  )
}

## Public fitter --------------------------------------------------------------
#' Fit the Quadratic Response Plateau (QRP) model by grid search
#'
#' Fits the quadratic-plateau (smooth broken-line) model
#' \deqn{f(x) = a + b x + c x^2 \ \text{ if } x \le X_0, \qquad
#'       f(x) = a - b^2/(4c) \ \text{ if } x > X_0,}
#' with \eqn{X_0 = -b/(2c)}. The breakpoint is profiled over a grid: for each
#' candidate \eqn{X_0} the model is linear in the plateau level and the
#' curvature, so it is fit by least squares with no starting values, and the
#' \eqn{X_0} minimizing the residual sum of squares is returned. This mirrors
#' [fit_lrp()] and avoids the convergence problems of a direct nls fit.
#'
#' @section Two ways to call:
#' \describe{
#'   \item{Vectors}{\code{fit_qrp(x, cv)} returns a \code{"qrp_fit"}.}
#'   \item{Data frame}{\code{fit_qrp(.data, x = "x", cv = "cv", trial = "trial")};
#'     with \code{trial}, one model per trial is fit and a \code{"qrp_multi"}
#'     object is returned. Column names default to "x", "cv", "trial".}
#' }
#'
#' @section Uncertainty of the breakpoint:
#' \code{bootstrap = TRUE} resamples the shapes (the rows of the CV table) with
#' replacement, refits on each resample, and returns a percentile interval for
#' \eqn{X_o} in \code{$bootstrap$ci}, together with a bootstrap standard error.
#'
#' It also tests whether the plateau is warranted. The null here is not the one
#' used by [fit_lrp()]: a quadratic-plateau that never plateaus is simply a
#' quadratic, so the null model is the unconstrained quadratic fitted on every
#' observation, and the p-value is the proportion of null resamples whose SSE
#' reduction matches or exceeds the observed one. A large p-value means the
#' plateau segment buys nothing over a plain quadratic, and \eqn{X_o} should not
#' be read as an optimal plot size.
#'
#' Set the random seed before calling to make the result reproducible.
#'
#' @param .data optional data frame; when supplied, \code{x}/\code{cv}/\code{trial}
#'   are column names.
#' @param x,cv numeric vectors, or column names when \code{.data} is a data frame.
#' @param trial optional column name identifying the trial.
#' @param step grid step for the breakpoint search (default 0.001).
#' @param search_range optional \code{c(lower, upper)} restricting the breakpoint
#'   search; must fall within the data range.
#' @param start optional single breakpoint value. The fit is unchanged, but the
#'   result also carries \code{$compat}: the local minimum of the basin
#'   containing \code{start}, i.e. what a gradient-based fitter seeded there
#'   would return. Rarely needed for the QRP, whose SSE profile is usually a
#'   single smooth basin.
#' @param local_min_tol relative SSE tolerance (default 0.10) deciding which
#'   local minima count as competing. Only labelling is affected (the
#'   \code{competing} column of \code{$local_minima} and the stars in
#'   \code{print()}); the fit never changes, and no warning is issued.
#' @param bootstrap logical; if \code{TRUE}, also estimate the uncertainty of
#'   the breakpoint by resampling (default \code{FALSE}). See the section
#'   "Uncertainty of the breakpoint".
#' @param n_boot number of bootstrap resamples (default 1000), used only when
#'   \code{bootstrap = TRUE}.
#' @param conf_level confidence level of the percentile interval (default 0.95),
#'   used only when \code{bootstrap = TRUE}.
#' @param weights weights for a weighted least-squares fit. \code{FALSE}
#'   (default) fits unweighted, as the published procedure does. \code{TRUE}
#'   uses the \code{n} column of \code{.data} (the number of plots behind each
#'   CV, as returned by [calc_cv_shapes()]); a column name or a numeric vector
#'   are also accepted. The caveats are the same as for [fit_lrp()]: the
#'   breakpoint typically falls, the points remain dependent, and the weighted
#'   fit statistics are not comparable with the unweighted ones.
#' @return A \code{"qrp_fit"} (single series) or \code{"qrp_multi"} (per trial).
#'   The fit also carries \code{local_minima} (competing basins with their SSE
#'   excess over the optimum and a \code{competing} flag, or \code{NULL}),
#'   \code{local_min_tol}, \code{sse_profile},
#'   \code{compat} when \code{start} was given, and \code{bootstrap} when
#'   \code{bootstrap = TRUE} (a list with \code{ci}, \code{se}, \code{p_value},
#'   \code{statistic}, \code{replicates}, \code{n_valid}, \code{conf_level} and
#'   \code{null_model}). Because the quadratic-plateau
#'   joins smoothly, this profile is typically a single basin, unlike the
#'   stepped profile of [fit_lrp()].
#' @seealso [fit_lrp()], [plot.qrp_fit()]
#' @examples
#' ## CV per plot shape from the bundled simulated uniformity trial
#' ## (see ?uniformity_trial).
#' grid1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
#'                                     grep("^col", names(uniformity_trial))])
#' cv_tab <- calc_cv_shapes(list(T1 = grid1))
#' X   <- cv_tab$x
#' CV1 <- cv_tab$cv
#'
#' ## A coarse grid runs fast; the default step = 0.001 refines the third decimal
#' fit <- fit_qrp(X, CV1, step = 0.01)
#' fit
#'
#' ## The quadratic joins its plateau smoothly, so the optimum is larger than
#' ## the one the broken-line model gives on the same data
#' coef(fit)
#' predict(fit, newx = c(2, 5, 7.5, 15))
#'
#' \donttest{
#' ## Full precision is the default; it costs about eight times more time and
#' ## only refines the third decimal, so it is shown here rather than used
#' ## throughout:
#' # fit_qrp(X, CV1)$parameters["Breakpoint"]
#'
#' plot(fit, title = "Uniformity trial, T1")
#'
#' ## The smooth join makes the SSE profile a single basin, unlike the stepped
#' ## profile of fit_lrp(): competing minima are rare and much worse.
#' plot(fit$sse_profile, type = "l", xlab = "Breakpoint", ylab = "SSE")
#' fit$local_minima
#'
#' ## Uncertainty of Xo, off by default. The p-value asks whether the plateau
#' ## buys anything over a plain quadratic.
#' set.seed(1)
#' unc <- fit_qrp(X, CV1, step = 0.01, bootstrap = TRUE, n_boot = 200)
#' unc
#' unc$bootstrap$ci
#'
#' ## Restrict the breakpoint search
#' fit_qrp(X, CV1, step = 0.01, search_range = c(5, 15))$parameters["Breakpoint"]
#'
#' ## One model per trial, with a summary table
#' trials <- rbind(
#'   data.frame(x = X, cv = CV1,        trial = "T1"),
#'   data.frame(x = X, cv = CV1 * 0.85, trial = "T2")
#' )
#' fit_qrp(trials, x = "x", cv = "cv", trial = "trial", step = 0.01)$summary
#' }
#' @export
fit_qrp <- function(.data = NULL, x = NULL, cv = NULL, trial = NULL,
                    step = 0.001, search_range = NULL, start = NULL,
                    local_min_tol = 0.10, bootstrap = FALSE,
                    n_boot = 1000, conf_level = 0.95, weights = FALSE) {

  if (!is.null(.data) && !is.data.frame(.data)) {
    cv <- x; x <- .data; .data <- NULL
  }

  if (is.null(.data)) {
    if (is.null(x) || is.null(cv))
      stop("Provide numeric `x` and `cv`, or a data frame as `.data`.",
           call. = FALSE)
    return(.qrp_fit_one(x, cv, step, search_range, start, local_min_tol,
                        bootstrap, n_boot, conf_level,
                        .resolve_weights(weights, NULL, length(x))))
  }

  if (!is.data.frame(.data)) stop("`.data` must be a data frame.", call. = FALSE)

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
    return(.qrp_fit_one(.data[[xcol]], .data[[cvcol]], step, search_range,
                        start, local_min_tol, bootstrap, n_boot, conf_level,
                        wvec))
  }

  groups <- split(.data, .data[[trial]])
  message(sprintf("Using x = '%s', cv = '%s', trial = '%s' -> %d trials.",
                  xcol, cvcol, trial, length(groups)))

  wsplit <- if (is.null(wvec)) NULL else split(wvec, .data[[trial]])
  fits <- Map(function(g, nm) .qrp_fit_one(g[[xcol]], g[[cvcol]],
                                           step, search_range, start,
                                           local_min_tol, bootstrap,
                                           n_boot, conf_level,
                                           if (is.null(wsplit)) NULL else wsplit[[nm]]),
              groups, names(groups))
  summ <- do.call(rbind, Map(function(f, nm) data.frame(
    trial      = nm,
    a          = round(unname(f$coefficients["a"]), 4),
    b          = round(unname(f$coefficients["b"]), 4),
    c          = round(unname(f$coefficients["c"]), 4),
    breakpoint = round(unname(f$parameters["Breakpoint"]), 4),
    plateau    = round(unname(f$parameters["Breakpoint_Response"]), 4),
    R2         = round(unname(f$parameters["R2"]), 4),
    RMSE       = round(unname(f$parameters["RMSE"]), 4),
    AIC        = round(unname(f$parameters["AIC"]), 3),
    BIC        = round(unname(f$parameters["BIC"]), 3),
    n_local    = if (is.null(f$local_minima)) 0L else nrow(f$local_minima),
    stringsAsFactors = FALSE), fits, names(fits)))
  rownames(summ) <- NULL

  if (isTRUE(bootstrap)) {
    summ$Xo_lwr <- round(vapply(fits, function(f) f$bootstrap$ci[1], 0), 4)
    summ$Xo_upr <- round(vapply(fits, function(f) f$bootstrap$ci[2], 0), 4)
    summ$p_breakpoint <- round(vapply(fits, function(f) f$bootstrap$p_value, 0), 4)
  }

  structure(list(fits = fits, summary = summ), class = "qrp_multi")
}

## Methods --------------------------------------------------------------------
#' Predictions from a QRP fit
#' @param object a \code{"qrp_fit"} object.
#' @param newx numeric predictor values; defaults to the fitted data.
#' @param ... ignored.
#' @return numeric vector of predicted responses.
#' @export
predict.qrp_fit <- function(object, newx = NULL, ...) {
  a  <- object$coefficients["a"]; b <- object$coefficients["b"]
  cc <- object$coefficients["c"]
  bp <- object$parameters["Breakpoint"]; pl <- object$parameters["Breakpoint_Response"]
  if (is.null(newx)) newx <- object$data$x
  unname(ifelse(newx <= bp, a + b * newx + cc * newx^2, pl))
}

#' @export
print.qrp_fit <- function(x, ...) {
  cat("Quadratic Response Plateau (QRP) fit\n")
  cat("Breakpoint (Xo):        ", sprintf("%.3f", x$parameters["Breakpoint"]), "\n")
  if (!is.null(x$bootstrap)) {
    b <- x$bootstrap
    cat(sprintf("  %.0f%% CI (percentile):  [%.3f, %.3f]   SE %.3f\n",
                100 * b$conf_level, b$ci[1], b$ci[2], b$se))
    cat(sprintf("  plateau vs quadratic:  p = %.4f  (%d resamples)\n",
                b$p_value, b$n_valid))
  }
  cat("CV at breakpoint:       ", sprintf("%.3f", x$parameters["Breakpoint_Response"]), "\n")
  cat("R2:", sprintf("%.3f", x$parameters["R2"]),
      " R2 adj:", sprintf("%.3f", x$parameters["R2_adj"]),
      " RMSE:", sprintf("%.3f", x$parameters["RMSE"]),
      " MAE:", sprintf("%.3f", x$parameters["MAE"]), "\n")

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

#' @export
summary.qrp_fit <- function(object, ...) {
  cat("Model coefficients:\n"); print(object$coefficients)
  cat("\nGoodness of fit:\n");   print(object$parameters)
  invisible(object)
}

#' Plot a QRP fit (publication style)
#'
#' Same layout and options as [plot.lrp_fit()], with the quadratic descending
#' arm and the quadratic equation in the annotation.
#'
#' @param x a \code{"qrp_fit"} object.
#' @param title plot title.
#' @param annotate_model draw the model equations and statistics.
#' @param xlab,ylab axis titles. \code{xlab = NULL} uses "Plot size (m^2)".
#' @param decimal_mark decimal separator, "." or ",".
#' @param cond_word conditional word in the equations (e.g. "if" or "se").
#' @param digits_coef,digits_c,digits_stat decimals for a/b, for c, and for the
#'   statistics.
#' @param point_size,line_size sizes of points and fitted line.
#' @param point_colour,line_colour,bp_colour point, line and breakpoint colours.
#' @param base_size,label_size,title_size,family theme and annotation sizes and
#'   font family.
#' @param theme optional ggplot2 theme used instead of the default.
#' @param save,file,format,dpi,width,height,units,compression saving controls;
#'   see [plot.lrp_fit()].
#' @param ... ignored.
#' @return A \code{ggplot} object (invisibly when saved).
#' @export
plot.qrp_fit <- function(x, title = "Quadratic Plateau", annotate_model = TRUE,
                         xlab = NULL, ylab = "CV (%)",
                         decimal_mark = c(".", ","), cond_word = "if",
                         digits_coef = 3, digits_c = 4, digits_stat = 2,
                         point_size = 2.3, line_size = 0.8,
                         point_colour = "black", line_colour = "black",
                         bp_colour = "red",
                         base_size = 12, label_size = 4.6, title_size = NULL,
                         family = "sans", theme = NULL,
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
  cc <- unname(x$coefficients["c"])
  xo <- unname(x$parameters["Breakpoint"])
  p  <- unname(x$parameters["Breakpoint_Response"])
  r2 <- unname(x$parameters["R2"])

  xmax <- max(d$x) * 1.05
  ## 1.15 leaves room for the quadratic term's exponent on the top annotation
  ## line, which would otherwise be clipped by the panel edge.
  ymax <- max(d$cv) * 1.15
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
    sgn_b <- if (b  < 0) "-" else "+"
    sgn_c <- if (cc < 0) "-" else "+"
    a_s  <- fmtn(a, digits_coef);      b_s <- fmtn(abs(b), digits_coef)
    c_s  <- fmtn(abs(cc), digits_c)
    xo_s <- fmtn(xo, digits_stat); p_s <- fmtn(p, digits_stat)
    r2_s <- fmtn(r2, digits_stat)

    eq1 <- paste0("CV[(x)]=='", a_s, "'", sgn_b, "'", b_s, "'*X", sgn_c,
                  "'", c_s, "'*X^2~~'", cond_word, "'~~X<='", xo_s, "'")
    eq2 <- paste0("CV[(x)]=='", p_s, "'~~'", cond_word, "'~~X>'", xo_s, "'")
    eq3 <- paste0("R^2=='", r2_s, "'")
    lXo <- paste0("X[o]=='", xo_s, "'")
    lCV <- paste0("CV[Xo]=='", p_s, "'")

    xm <- 0.55 * xmax
    ym <- ymax * c(0.93, 0.85, 0.77)
    for (i in seq_len(3))
      g <- g + ggplot2::annotate("text", x = xm, y = ym[i],
                                 label = c(eq1, eq2, eq3)[i], parse = TRUE,
                                 hjust = 0.5, size = label_size, family = family)

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
      .theme_trialsizing(base_size = base_size, family = family)
     else theme) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = title_size))

  if (save) {
    .save_lrp(g, file, title, format, dpi, width, height, units, compression)
    return(invisible(g))
  }
  g
}

#' @export
print.qrp_multi <- function(x, ...) {
  cat(sprintf("QRP fits for %d trials\n\n", length(x$fits)))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
summary.qrp_multi <- function(object, ...) {
  print(object$summary, row.names = FALSE)
  invisible(object)
}

#' Plot QRP fits for several trials (paginated grid)
#'
#' Arranges the per-trial plots in a grid of at most 6 panels (3 x 2), paginating
#' when there are more. Requires \pkg{patchwork}. See [plot.lrp_multi()].
#'
#' @param x a \code{"qrp_multi"} object.
#' @param save,file,format,dpi,width,height,units,compression saving controls.
#' @param ... styling arguments forwarded to [plot.qrp_fit()].
#' @return A list of page objects, invisibly.
#' @export
plot.qrp_multi <- function(x, save = FALSE, file = NULL,
                           format = c("tiff", "png", "jpeg", "pdf", "eps"),
                           dpi = 300, width = 18, height = 24, units = "cm",
                           compression = "lzw", ...) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required to arrange multiple trials.",
         call. = FALSE)
  format <- match.arg(format)

  plots <- Map(function(f, nm) plot(f, title = as.character(nm), save = FALSE, ...),
               x$fits, names(x$fits))
  idx    <- seq_along(plots)
  chunks <- split(idx, ceiling(idx / 6))
  pages  <- lapply(chunks, function(ix)
    patchwork::wrap_plots(plots[ix], nrow = 3, ncol = 2))

  if (save) {
    base <- if (is.null(file)) "qrp_trials" else sub("\\.[^.]*$", "", file)
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
