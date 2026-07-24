## ============================================================================
## trialSize :: pre-flight check of a uniformity trial
## ----------------------------------------------------------------------------
## Everything that should be looked at before a plot-size method is run: is the
## grid complete and usable, are the values sane, and does the field actually
## have spatial structure worth sizing plots against.
##
## The geostatistics here is written from the standard definitions rather than
## delegating to a spatial package, which keeps the dependency list at ggplot2
## and keeps the variogram parameters -- nugget, sill, range -- available as
## diagnostics rather than as internals of a plotting routine.
## References: Matheron (1963); Cressie (1993); Cambardella et al. (1994);
## Moran (1950); Webster & Oliver (2007).
## ============================================================================

## Variogram model shapes (internal) ------------------------------------------
## Each returns the structured part on a 0-1 scale, so gamma(h) = c0 + c1 g(h).
## Exponential and Gaussian use the effective-range convention: `a` is the
## distance at which the model reaches 95% of its sill.
.vg_shape <- function(h, a, model) {
  switch(model,
    spherical   = ifelse(h >= a, 1, 1.5 * (h / a) - 0.5 * (h / a)^3),
    exponential = 1 - exp(-3 * h / a),
    gaussian    = 1 - exp(-3 * (h / a)^2),
    stop("Unknown variogram model.", call. = FALSE))
}

## Empirical variogram (internal) ---------------------------------------------
## gamma(h) = mean squared difference / 2, over pairs binned by separation.
.vg_empirical <- function(z, xy, n_bins, max_dist = NULL) {
  d  <- stats::dist(xy)
  dz <- stats::dist(z)
  d  <- as.numeric(d); sq <- as.numeric(dz)^2

  if (is.null(max_dist)) max_dist <- max(d) / 2   # the usual half-extent rule
  keep <- d > 0 & d <= max_dist
  d <- d[keep]; sq <- sq[keep]
  if (!length(d))
    stop("No pairs of points within the variogram range.", call. = FALSE)

  brk <- seq(0, max_dist, length.out = n_bins + 1)
  bin <- cut(d, brk, include.lowest = TRUE)
  agg <- data.frame(
    dist  = tapply(d, bin, mean),
    gamma = tapply(sq, bin, mean) / 2,
    n     = tapply(sq, bin, length))
  agg <- agg[!is.na(agg$dist) & agg$n > 0, , drop = FALSE]
  rownames(agg) <- NULL
  agg
}

## Variogram model fit (internal) ---------------------------------------------
## For a fixed range the model is linear in the nugget and the partial sill, so
## the range is profiled over a grid and the linear part solved exactly at each
## candidate -- the same device the plateau fitters use, and for the same
## reason: no starting values and no convergence failures.
.vg_fit <- function(emp, models = c("spherical", "exponential", "gaussian"),
                    n_grid = 200) {
  y <- emp$gamma; w <- emp$n; h <- emp$dist
  if (nrow(emp) < 3)
    return(NULL)

  ## non-negative weighted least squares in (c0, c1) for a fixed shape
  solve_linear <- function(g) {
    S <- c(sum(w), sum(w * g), sum(w * g), sum(w * g * g))
    b <- c(sum(w * y), sum(w * g * y))
    det <- S[1] * S[4] - S[2] * S[3]
    cand <- list()
    if (is.finite(det) && abs(det) > 1e-12) {
      c0 <- (b[1] * S[4] - b[2] * S[2]) / det
      c1 <- (S[1] * b[2] - S[3] * b[1]) / det
      if (c0 >= 0 && c1 >= 0) cand[[length(cand) + 1L]] <- c(c0, c1)
    }
    ## boundary solutions: pure nugget, or no nugget
    cand[[length(cand) + 1L]] <- c(max(sum(w * y) / sum(w), 0), 0)
    den <- sum(w * g * g)
    if (den > 0) cand[[length(cand) + 1L]] <- c(0, max(sum(w * g * y) / den, 0))
    sse <- vapply(cand, function(p) sum(w * (y - p[1] - p[2] * g)^2), numeric(1))
    cand[[which.min(sse)]]
  }

  rng <- seq(min(h), max(h) * 1.5, length.out = n_grid)
  best <- NULL
  for (mod in models) for (a in rng) {
    p  <- solve_linear(.vg_shape(h, a, mod))
    ss <- sum(w * (y - p[1] - p[2] * .vg_shape(h, a, mod))^2)
    if (is.null(best) || ss < best$sse)
      best <- list(model = mod, nugget = p[1], psill = p[2], range = a, sse = ss)
  }

  best$sill <- best$nugget + best$psill
  ## Cambardella et al. (1994): the share of variance that is not spatially
  ## structured, and the usual reading of it
  best$nugget_ratio <- if (best$sill > 0) best$nugget / best$sill else NA_real_
  best$dependence <- if (is.na(best$nugget_ratio)) NA_character_
  else if (best$nugget_ratio < 0.25) "strong"
  else if (best$nugget_ratio <= 0.75) "moderate" else "weak"
  best
}

## Ordinary kriging (internal) ------------------------------------------------
## The left-hand matrix does not depend on the prediction location, so it is
## factorized once and every target solved in a single pass.
.krige <- function(z, xy, newxy, vg) {
  n <- nrow(xy)
  gam <- function(h) {
    g <- vg$nugget + vg$psill * .vg_shape(h, vg$range, vg$model)
    g[h == 0] <- 0
    g
  }
  D <- as.matrix(stats::dist(xy))
  A <- matrix(0, n + 1, n + 1)
  A[seq_len(n), seq_len(n)] <- gam(D)
  A[n + 1, seq_len(n)] <- 1
  A[seq_len(n), n + 1] <- 1

  ## distances from every observation to every target
  d0 <- outer(xy[, 1], newxy[, 1], "-")^2 + outer(xy[, 2], newxy[, 2], "-")^2
  B  <- rbind(gam(sqrt(d0)), 1)

  lam <- tryCatch(solve(A, B), error = function(e) NULL)
  if (is.null(lam))                       # singular: fall back to the mean
    return(rep(mean(z), nrow(newxy)))
  as.numeric(crossprod(lam[seq_len(n), , drop = FALSE], z))
}

## Moran's I (internal) -------------------------------------------------------
## Rook adjacency on the grid, normal-approximation p-value.
.morans_i <- function(m) {
  nr <- nrow(m); nc <- ncol(m); n <- nr * nc
  z  <- as.numeric(m) - mean(m)
  idx <- function(i, j) (j - 1) * nr + i

  num <- 0; S0 <- 0
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    for (nb in list(c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1))) {
      if (nb[1] < 1 || nb[1] > nr || nb[2] < 1 || nb[2] > nc) next
      num <- num + z[idx(i, j)] * z[idx(nb[1], nb[2])]
      S0  <- S0 + 1
    }
  }
  den <- sum(z^2)
  if (den == 0 || S0 == 0) return(list(I = NA_real_, p_value = NA_real_))

  I  <- (n / S0) * num / den
  EI <- -1 / (n - 1)
  ## binary symmetric weights: S1 = 2 S0; S2 needs the row sums
  deg <- numeric(n)
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    d <- 0
    for (nb in list(c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1)))
      if (!(nb[1] < 1 || nb[1] > nr || nb[2] < 1 || nb[2] > nc)) d <- d + 1
    deg[idx(i, j)] <- d
  }
  S1 <- 2 * S0
  S2 <- sum((2 * deg)^2)
  VI <- (n^2 * S1 - n * S2 + 3 * S0^2) / ((n^2 - 1) * S0^2) - EI^2
  zscore <- if (VI > 0) (I - EI) / sqrt(VI) else NA_real_
  list(I = I, expected = EI,
       p_value = if (is.na(zscore)) NA_real_ else 2 * stats::pnorm(-abs(zscore)))
}

## One trial (internal) -------------------------------------------------------
.check_one <- function(m, trial_id, cell_size, n_bins, max_dist, fit_variogram) {
  nr <- nrow(m); nc <- ncol(m)
  v  <- as.numeric(m)
  ok <- !is.na(v)

  divisors <- function(k) seq_len(k)[k %% seq_len(k) == 0]
  n_shapes <- length(divisors(nr)) * length(divisors(nc)) - 1L

  ## ---- values ----
  issues <- character(0)
  if (any(is.na(v)))
    issues <- c(issues, sprintf("%d missing value(s)", sum(is.na(v))))
  if (any(v[ok] < 0))
    issues <- c(issues, sprintf("%d negative value(s)", sum(v[ok] < 0)))
  if (any(v[ok] == 0))
    issues <- c(issues, sprintf("%d zero(s)", sum(v[ok] == 0)))
  if (sum(ok) > 1 && stats::sd(v[ok]) == 0)
    issues <- c(issues, "all values identical (CV is undefined)")
  if (n_shapes < 6)
    issues <- c(issues, sprintf("only %d plot shape(s) available from a %d x %d grid",
                                n_shapes, nr, nc))

  mu <- mean(v[ok]); s <- stats::sd(v[ok])
  ## outliers by the boxplot rule, which does not assume normality
  q  <- stats::quantile(v[ok], c(0.25, 0.75))
  iqr <- q[2] - q[1]
  out_idx <- which(ok & (v < q[1] - 1.5 * iqr | v > q[2] + 1.5 * iqr))
  outliers <- if (length(out_idx))
    data.frame(row = ((out_idx - 1) %% nr) + 1, col = ((out_idx - 1) %/% nr) + 1,
               value = v[out_idx], z = (v[out_idx] - mu) / s)
  else NULL

  ## ---- trend along rows and columns ----
  row_means <- rowMeans(m, na.rm = TRUE)
  col_means <- colMeans(m, na.rm = TRUE)
  trend_p <- function(y) {
    if (length(y) < 3 || stats::sd(y) == 0) return(NA_real_)
    stats::cor.test(seq_along(y), y)$p.value
  }

  ## ---- spatial structure ----
  mi <- if (anyNA(m)) list(I = NA_real_, p_value = NA_real_) else .morans_i(m)

  ## the serpentine autocorrelations the Paranaiba method uses
  ## the serpentine walk needs an unbroken sequence, so a gap makes the
  ## first-order autocorrelation undefined rather than merely awkward
  rho_of <- function(u) {
    u <- as.numeric(u)
    if (length(u) < 3 || anyNA(u) || stats::sd(u) == 0) return(NA_real_)
    unname(stats::acf(u, lag.max = 1, plot = FALSE)$acf[2])
  }
  rho_row <- rho_of(t(m)); rho_col <- rho_of(m)

  xy <- cbind(rep(seq_len(nr), times = nc) * cell_size[1],
              rep(seq_len(nc), each  = nr) * cell_size[2])
  emp <- vg <- NULL
  if (isTRUE(fit_variogram) && sum(ok) >= 9 && s > 0) {
    emp <- .vg_empirical(v[ok], xy[ok, , drop = FALSE], n_bins, max_dist)
    vg  <- .vg_fit(emp)
  }

  list(
    trial = trial_id, matrix = m, coords = xy,
    dim = c(rows = nr, cols = nc), n_cells = nr * nc, n_missing = sum(!ok),
    n_shapes = n_shapes,
    stats = c(mean = mu, sd = s, cv = 100 * s / mu,
              min = min(v[ok]), max = max(v[ok])),
    outliers = outliers,
    trend = list(row_means = row_means, col_means = col_means,
                 p_row = trend_p(row_means), p_col = trend_p(col_means)),
    spatial = list(morans_i = mi$I, morans_p = mi$p_value,
                   rho_row = rho_row, rho_col = rho_col),
    variogram = vg, variogram_empirical = emp,
    issues = issues)
}

## Public function ------------------------------------------------------------
#' Check a uniformity trial before sizing plots
#'
#' Inspects the raw grid of basic experimental units and reports what should be
#' settled before any plot-size method is run: whether the grid is complete and
#' usable, whether the values are sane, and whether the field has spatial
#' structure worth sizing plots against. The result carries a \code{plot()}
#' method drawing the field map.
#'
#' @section What is checked:
#' \describe{
#'   \item{Structure}{Grid dimensions, missing cells, and how many rectangular
#'     plot shapes the grid admits. A grid whose sides are prime yields almost
#'     no shapes, which stops the CV-based methods before they start; the count
#'     is flagged when it falls below six.}
#'   \item{Values}{Missing, negative, zero and constant values, and outliers by
#'     the boxplot rule. In a uniformity trial an outlying basic unit is usually
#'     a harvest failure or a typing error, and it contaminates every shape that
#'     contains it.}
#'   \item{Trend}{Means along rows and columns, with the p-value of a monotone
#'     trend. A gradient in one direction is field fertility, and it is why
#'     shapes of equal area but different orientation give different CVs.}
#'   \item{Spatial structure}{Moran's I with its p-value, the first-order
#'     autocorrelations used by [calc_paranaiba()], and a fitted variogram.}
#' }
#'
#' @section Why the variogram matters here:
#' The fitted model gives three numbers that bear directly on plot size. The
#' \emph{range} is the distance beyond which basic units stop being correlated:
#' a plot larger than the range is averaging units that are already
#' independent, which is the spatial reading of the plateau that [fit_lrp()]
#' and [fit_qrp()] estimate empirically. The \emph{nugget-to-sill ratio} is the
#' share of variance with no spatial structure; following Cambardella et al.
#' (1994) it is read as strong (below 0.25), moderate (0.25 to 0.75) or weak
#' (above 0.75) spatial dependence. When it is weak, the field varies almost at
#' random from unit to unit and no choice of plot size will buy much precision
#' -- worth knowing before fitting five models to the CV table.
#'
#' The model is fitted by profiling the range over a grid and solving the
#' nugget and partial sill exactly at each candidate, under non-negativity.
#' Spherical, exponential and Gaussian shapes are all tried and the best
#' weighted fit is kept; exponential and Gaussian use the effective-range
#' convention, where the range is the distance reaching 95\% of the sill.
#'
#' @param .data a matrix (one trial), a named list of matrices, or a long data
#'   frame with one row per basic unit.
#' @param value name of the measurement column (long format).
#' @param row_id,col_id names of the row and column index columns (long format).
#' @param trial optional column name identifying the trial.
#' @param cell_size numeric \code{c(row, col)}: the distance between the centres
#'   of adjacent basic units, in the units you want distances reported in.
#'   Default \code{c(1, 1)}.
#' @param n_bins number of distance classes for the empirical variogram.
#' @param max_dist largest separation used; \code{NULL} (default) takes half the
#'   maximum distance in the field, the usual rule.
#' @param variogram logical; fit the variogram (default \code{TRUE}). Turning it
#'   off skips the only appreciable computation.
#'
#' @return An object of class \code{"trial_check"}: a list with \code{checks}
#'   (one element per trial, each holding the dimensions, statistics, outliers,
#'   trend, spatial measures, fitted \code{variogram} and any \code{issues}) and
#'   \code{summary}, one row per trial.
#'
#' @references
#' Matheron, G. (1963). Principles of geostatistics. \emph{Economic Geology},
#' 58(8), 1246-1266. \cr
#' Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
#' \emph{Biometrika}, 37(1/2), 17-23. \cr
#' Cambardella, C. A. et al. (1994). Field-scale variability of soil properties
#' in central Iowa soils. \emph{Soil Science Society of America Journal}, 58(5),
#' 1501-1511. \cr
#' Cressie, N. (1993). \emph{Statistics for Spatial Data}. Wiley, New York. \cr
#' Webster, R. & Oliver, M. A. (2007). \emph{Geostatistics for Environmental
#' Scientists}, 2nd ed. Wiley, Chichester.
#'
#' @seealso [calc_cv_shapes()], [calc_paranaiba()], [plot.trial_check()]
#' @examples
#' grid1 <- as.matrix(dados_ensaio_C1[dados_ensaio_C1$Rep == 1,
#'                                    paste0("C", 1:8)])
#' chk <- check_trial(grid1)
#' chk
#'
#' ## the fitted variogram, as numbers rather than as a picture
#' chk$checks[[1]]$variogram[c("model", "nugget", "sill", "range",
#'                             "nugget_ratio", "dependence")]
#'
#' \donttest{
#' ## the field map: kriged surface with the basic units drawn on top
#' plot(chk)
#'
#' ## several trials at once
#' grids <- lapply(split(dados_ensaio_C1, dados_ensaio_C1$Rep),
#'                 function(d) as.matrix(d[, paste0("C", 1:8)]))
#' names(grids) <- paste("Rep", names(grids))
#' check_trial(grids)$summary
#'
#' ## a grid whose sides are prime admits almost no plot shapes
#' check_trial(matrix(rnorm(77, 100, 10), nrow = 7))$checks[[1]]$issues
#' }
#' @export
check_trial <- function(.data, value = NULL, row_id = NULL, col_id = NULL,
                        trial = NULL, cell_size = c(1, 1), n_bins = 15,
                        max_dist = NULL, variogram = TRUE) {

  if (!is.numeric(cell_size) || length(cell_size) != 2 || any(cell_size <= 0))
    stop("`cell_size` must be two positive numbers, c(row, col).", call. = FALSE)
  if (!is.numeric(n_bins) || length(n_bins) != 1 || n_bins < 3)
    stop("`n_bins` must be a single number of at least 3.", call. = FALSE)

  mats <- .as_grid_list(.data, value, row_id, col_id, trial)
  for (nm in names(mats))
    if (!is.numeric(mats[[nm]]))
      stop(sprintf("Trial '%s' is not a numeric matrix.", nm), call. = FALSE)

  message(sprintf("Checking %d trial(s).", length(mats)))

  checks <- Map(function(mm, nm)
    .check_one(mm, nm, cell_size, n_bins, max_dist, variogram),
    mats, names(mats))

  summ <- do.call(rbind, lapply(checks, function(k) data.frame(
    trial = k$trial, rows = k$dim[["rows"]], cols = k$dim[["cols"]],
    n_missing = k$n_missing, n_shapes = k$n_shapes,
    mean = k$stats[["mean"]], cv = k$stats[["cv"]],
    n_outliers = if (is.null(k$outliers)) 0L else nrow(k$outliers),
    morans_i = k$spatial$morans_i,
    range = if (is.null(k$variogram)) NA_real_ else k$variogram$range,
    nugget_ratio = if (is.null(k$variogram)) NA_real_ else k$variogram$nugget_ratio,
    dependence = if (is.null(k$variogram)) NA_character_ else k$variogram$dependence,
    n_issues = length(k$issues),
    stringsAsFactors = FALSE)))
  rownames(summ) <- NULL

  structure(list(checks = checks, summary = summ,
                 meta = list(cell_size = cell_size, n_trials = length(mats))),
            class = "trial_check")
}

## Methods --------------------------------------------------------------------
#' Print a trial check
#'
#' @param x an object of class \code{"trial_check"}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.trial_check <- function(x, ...) {
  for (k in x$checks) {
    cat("Uniformity trial check --", k$trial, "\n")
    cat(sprintf("  Grid:      %d x %d = %d basic units, %d missing\n",
                k$dim[["rows"]], k$dim[["cols"]], k$n_cells, k$n_missing))
    cat(sprintf("  Shapes:    %d rectangular plot shapes available\n", k$n_shapes))
    cat(sprintf("  Values:    mean %.3f, sd %.3f, CV %.2f%%, range [%.3f, %.3f]\n",
                k$stats[["mean"]], k$stats[["sd"]], k$stats[["cv"]],
                k$stats[["min"]], k$stats[["max"]]))
    cat(sprintf("  Outliers:  %d by the boxplot rule\n",
                if (is.null(k$outliers)) 0L else nrow(k$outliers)))

    tp <- k$trend
    cat(sprintf("  Trend:     rows p = %s, columns p = %s\n",
                format.pval(tp$p_row, digits = 3),
                format.pval(tp$p_col, digits = 3)))
    sp <- k$spatial
    cat(sprintf("  Moran's I: %s (p = %s)   rho row %s, col %s\n",
                formatC(sp$morans_i, format = "f", digits = 3),
                format.pval(sp$morans_p, digits = 3),
                formatC(sp$rho_row, format = "f", digits = 3),
                formatC(sp$rho_col, format = "f", digits = 3)))

    v <- k$variogram
    if (is.null(v)) {
      cat("  Variogram: not fitted\n")
    } else {
      cat(sprintf("  Variogram: %s | nugget %.3f, sill %.3f, range %.2f\n",
                  v$model, v$nugget, v$sill, v$range))
      cat(sprintf("             nugget/sill %.2f -> %s spatial dependence\n",
                  v$nugget_ratio, v$dependence))
      if (identical(v$dependence, "weak"))
        cat("             little spatial structure: plot size will buy little\n",
            "             precision here, whatever method is used\n", sep = "")
      if (!is.null(k$outliers))
        cat("             note: squared differences drive the variogram, so the\n",
            "             outlier(s) above inflate it; check them first\n", sep = "")
    }

    if (length(k$issues)) {
      cat("  Issues:\n")
      for (i in k$issues) cat("    ! ", i, "\n", sep = "")
    } else {
      cat("  Issues:    none\n")
    }
    cat("\n")
  }
  invisible(x)
}

#' Summarize a trial check
#'
#' @param object an object of class \code{"trial_check"}.
#' @param ... ignored.
#' @return The summary data frame, invisibly.
#' @export
summary.trial_check <- function(object, ...) {
  print(object$summary, row.names = FALSE)
  invisible(object$summary)
}
