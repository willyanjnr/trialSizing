## ============================================================================
## trialSizing :: side-by-side comparison of the plot-size methods
## ----------------------------------------------------------------------------
## The methods disagree by construction, and the published articles report that
## disagreement as a finding. This runs them on one trial and tabulates what
## each recommends, so the spread is visible instead of inferred from separate
## calls.
## ============================================================================

## Comparable fit statistics (internal) ---------------------------------------
## Deliberately recomputed here rather than read from the fit objects: a
## weighted fit stores weighted statistics, and comparing those across methods
## would compare different quantities. On the original scale every method is
## predicting the same CV values, so these are apples to apples.
.comparable_stats <- function(f) {
  r  <- f$residuals
  y  <- f$data$cv
  c(R2 = 1 - sum(r^2) / sum((y - mean(y))^2),
    RMSE = sqrt(mean(r^2)))
}

## One trial, every method (internal) -----------------------------------------
.compare_one <- function(x, cv, w, grid_mat, trial_id, methods, step,
                         bootstrap, n_boot, conf_level) {
  rows <- list()

  add <- function(method, xo, cvxo, stats = c(R2 = NA, RMSE = NA),
                  ci = c(NA, NA), p = NA) {
    rows[[length(rows) + 1L]] <<- data.frame(
      trial = trial_id, method = method, Xo = xo, CVxo = cvxo,
      R2 = unname(stats["R2"]), RMSE = unname(stats["RMSE"]),
      Xo_lwr = ci[1], Xo_upr = ci[2], p_breakpoint = p,
      stringsAsFactors = FALSE)
  }

  if ("mcm" %in% methods) {
    ## the MCM weights through its Federer `df` argument, the same weighted
    ## least squares by another name
    f <- suppressWarnings(.mcm_fit_one(x, cv, "nls", w, bootstrap, n_boot,
                                       conf_level))
    add("MCM", unname(f$parameters["Breakpoint"]),
        unname(f$parameters["Breakpoint_Response"]), .comparable_stats(f),
        if (is.null(f$bootstrap)) c(NA, NA) else f$bootstrap$ci,
        NA)   # no existence test: the curve always has a curvature maximum
  }

  if ("lrp" %in% methods) {
    f <- suppressWarnings(.lrp_fit_one(x, cv, step, "segment", NULL, NULL, 0.10,
                                       bootstrap, n_boot, conf_level, w))
    add("LRP", unname(f$parameters["Breakpoint"]),
        unname(f$parameters["Breakpoint_Response"]), .comparable_stats(f),
        if (is.null(f$bootstrap)) c(NA, NA) else f$bootstrap$ci,
        if (is.null(f$bootstrap)) NA else f$bootstrap$p_value)
  }

  if ("qrp" %in% methods) {
    f <- suppressWarnings(.qrp_fit_one(x, cv, step, NULL, NULL, 0.10,
                                       bootstrap, n_boot, conf_level, w))
    add("QRP", unname(f$parameters["Breakpoint"]),
        unname(f$parameters["Breakpoint_Response"]), .comparable_stats(f),
        if (is.null(f$bootstrap)) c(NA, NA) else f$bootstrap$ci,
        if (is.null(f$bootstrap)) NA else f$bootstrap$p_value)
  }

  if ("paranaiba" %in% methods && !is.null(grid_mat)) {
    ## closed form on the raw grid: no residuals, so no R2 or RMSE
    s <- .paranaiba_one(grid_mat, "row", trial_id)
    add("Paranaiba", s$Xo, s$CVxo)
  }

  do.call(rbind, rows)
}

## Public function ------------------------------------------------------------
#' Compare the plot-size methods on the same trial
#'
#' Runs the CV-based methods on one set of data and tabulates what each
#' recommends: the optimal plot size \eqn{X_o}, the CV at that size, and
#' comparable fit statistics. Given the raw grid instead of a CV table, it also
#' builds the table with [calc_cv_shapes()] and adds [calc_paranaiba()], which
#' works on the basic units rather than on CV values.
#'
#' @section What the table is for:
#' The methods disagree systematically: the MCM optimum is the smallest, the
#' LRP intermediate and the QRP the largest, an ordering reported across many
#' crops. Seeing the three side by side, with intervals when
#' \code{bootstrap = TRUE}, shows whether that ordering is a real difference or
#' three readings of the same imprecise quantity.
#'
#' @section Which statistics are comparable:
#' \code{R2} and \code{RMSE} are recomputed here from the residuals on the
#' original CV scale, so they mean the same thing for every method even when
#' the fit itself was weighted. AIC and BIC are deliberately absent: the MCM has
#' two parameters against the plateau models' three plus a breakpoint, and the
#' breakpoint is not an ordinary parameter, so the information criteria are not
#' on a common footing. The Paranaiba estimate is a closed form with no fitted
#' residuals, so its \code{R2} and \code{RMSE} are \code{NA} by nature, not by
#' omission.
#'
#' @param .data one of: a data frame holding a CV table (with \code{x} and
#'   \code{cv} columns); a matrix, being the raw grid of basic experimental
#'   units of one trial; or a named list of such matrices. Grids are expanded
#'   into CV tables with [calc_cv_shapes()], and only then can the Paranaiba
#'   method take part.
#' @param x,cv column names in \code{.data} when it is a data frame (default
#'   \code{"x"} and \code{"cv"}, the names [calc_cv_shapes()] produces).
#' @param trial optional column name identifying the trial.
#' @param methods which methods to run; \code{NULL} (default) runs every
#'   applicable one. Any subset of \code{"mcm"}, \code{"lrp"}, \code{"qrp"} and
#'   \code{"paranaiba"}, the last requiring a grid.
#' @param step grid step for the breakpoint search of the LRP and QRP.
#' @param weights weighting for the fits, as in [fit_lrp()]. \code{TRUE} uses
#'   the \code{n} column, which a grid always provides. It reaches the MCM
#'   through that method's Federer \code{df} argument, the same weighted least
#'   squares by another name. The Paranaiba row is unaffected: it is a closed
#'   form over the basic units and never sees the CV table.
#' @param bootstrap logical; add a confidence interval for each \eqn{X_o}, and
#'   the breakpoint-existence p-value where the method has one.
#' @param n_boot,conf_level bootstrap size and confidence level.
#'
#' @return An object of class \code{"method_comparison"}: a list with
#'   \code{summary} (one row per method and trial: \code{trial},
#'   \code{method}, \code{Xo}, \code{CVxo}, \code{R2}, \code{RMSE}, and
#'   \code{Xo_lwr}, \code{Xo_upr}, \code{p_breakpoint} when
#'   \code{bootstrap = TRUE}) and \code{meta}.
#'
#' @references
#' Cargnelutti Filho, A. et al. (2025). \emph{Revista Vivencias}, 21(43),
#' 499-513, which compares the same three methods and reports 4.81, 7.19 and
#' 10.25 m2 for the MCM, LRP and QRP respectively.
#'
#' @seealso [fit_lrp()], [fit_qrp()], [fit_mcm()], [calc_paranaiba()],
#'   [calc_cv_shapes()]
#' @examples
#' ## The bundled simulated uniformity trial (see ?uniformity_trial): one
#' ## 8 x 12 grid of 1 m2 basic units per trial.
#' grid1 <- as.matrix(uniformity_trial[uniformity_trial$trial == "T1",
#'                                     grep("^col", names(uniformity_trial))])
#'
#' ## From a CV table (built here from the grid)
#' cv_tab <- calc_cv_shapes(list(T1 = grid1))
#' cmp <- compare_methods(data.frame(x = cv_tab$x, cv = cv_tab$cv), step = 0.05)
#' cmp
#'
#' \donttest{
#' ## From the raw grid: the CV table is built on the way, and the Paranaiba
#' ## method joins because it needs the basic units.
#' compare_methods(grid1, step = 0.05)
#'
#' ## Weighted by the number of plots per shape
#' compare_methods(grid1, step = 0.05, weights = TRUE)$summary
#'
#' ## With intervals, which is the point: the spread between methods is often
#' ## smaller than the uncertainty within each one.
#' set.seed(1)
#' compare_methods(grid1, step = 0.05, bootstrap = TRUE, n_boot = 200)
#'
#' ## Several trials at once
#' grids <- lapply(split(uniformity_trial, uniformity_trial$trial),
#'                 function(d) as.matrix(d[, grep("^col", names(d))]))
#' compare_methods(grids, step = 0.05)$summary
#' }
#' @export
compare_methods <- function(.data, x = NULL, cv = NULL, trial = NULL,
                            methods = NULL, step = 0.001, weights = FALSE,
                            bootstrap = FALSE, n_boot = 1000,
                            conf_level = 0.95) {

  all_methods <- c("mcm", "lrp", "qrp", "paranaiba")

  ## ---- grid or CV table? ----
  from_grid <- is.matrix(.data) || (is.list(.data) && !is.data.frame(.data))
  grids <- NULL

  if (from_grid) {
    grids <- .as_grid_list(.data)
    tab   <- suppressMessages(calc_cv_shapes(.data))
    xcol <- "x"; cvcol <- "cv"; trialcol <- "trial"
  } else {
    if (!is.data.frame(.data))
      stop("`.data` must be a data frame (a CV table), a matrix, or a list of ",
           "matrices (raw grids).", call. = FALSE)
    tab      <- .data
    xcol     <- if (is.null(x))  "x"  else x
    cvcol    <- if (is.null(cv)) "cv" else cv
    trialcol <- trial
    need <- c(xcol, cvcol, if (!is.null(trialcol)) trialcol)
    missing_cols <- setdiff(need, names(tab))
    if (length(missing_cols))
      stop(sprintf("Column(s) not found in `.data`: %s.\n  Available: %s.",
                   paste(missing_cols, collapse = ", "),
                   paste(names(tab), collapse = ", ")), call. = FALSE)
  }

  ## ---- which methods ----
  if (is.null(methods)) {
    methods <- if (from_grid) all_methods else setdiff(all_methods, "paranaiba")
  } else {
    methods <- tolower(methods)
    bad <- setdiff(methods, all_methods)
    if (length(bad))
      stop(sprintf("Unknown method(s): %s. Choose from %s.",
                   paste(bad, collapse = ", "),
                   paste(all_methods, collapse = ", ")), call. = FALSE)
    if ("paranaiba" %in% methods && !from_grid)
      stop("The Paranaiba method works on the raw grid of basic units, so it ",
           "needs a matrix or a list of matrices, not a CV table.",
           call. = FALSE)
  }

  wvec <- .resolve_weights(weights, tab, nrow(tab))

  ## ---- split into trials ----
  if (is.null(trialcol) || !trialcol %in% names(tab)) {
    idx <- list(`Trial 1` = seq_len(nrow(tab)))
  } else {
    idx <- split(seq_len(nrow(tab)), as.character(tab[[trialcol]]))
  }

  message(sprintf("Comparing %d method(s) on %d trial(s).",
                  length(methods), length(idx)))

  out <- do.call(rbind, Map(function(rows, nm)
    .compare_one(tab[[xcol]][rows], tab[[cvcol]][rows],
                 if (is.null(wvec)) NULL else wvec[rows],
                 if (is.null(grids)) NULL else grids[[nm]],
                 nm, methods, step, bootstrap, n_boot, conf_level),
    idx, names(idx)))
  rownames(out) <- NULL

  ## the interval columns only earn their place when they hold something
  if (!isTRUE(bootstrap))
    out <- out[, setdiff(names(out), c("Xo_lwr", "Xo_upr", "p_breakpoint"))]

  structure(list(summary = out,
                 meta = list(methods = methods, from_grid = from_grid,
                             weighted = !is.null(wvec), step = step,
                             bootstrap = isTRUE(bootstrap),
                             n_trials = length(idx))),
            class = "method_comparison")
}

## Methods --------------------------------------------------------------------
#' Print a method comparison
#'
#' @param x an object of class \code{"method_comparison"}.
#' @param digits decimals for the printed table.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.method_comparison <- function(x, digits = 3, ...) {
  m <- x$meta
  cat("Plot-size methods compared\n")
  cat("Trials:", m$n_trials, " | source:",
      if (m$from_grid) "raw grid" else "CV table",
      " | weighted:", m$weighted, "\n\n")

  tab <- x$summary
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], round, digits)
  print(tab, row.names = FALSE)

  ## the headline: how far apart the recommendations are
  xo <- x$summary$Xo
  if (sum(is.finite(xo)) > 1) {
    cat(sprintf("\nXo ranges from %.*f to %.*f (a factor of %.1f).\n",
                digits, min(xo, na.rm = TRUE), digits, max(xo, na.rm = TRUE),
                max(xo, na.rm = TRUE) / min(xo, na.rm = TRUE)))
    if (m$bootstrap) {
      ## overlapping intervals mean the methods are not really disagreeing
      lo <- max(x$summary$Xo_lwr, na.rm = TRUE)
      hi <- min(x$summary$Xo_upr, na.rm = TRUE)
      cat(if (is.finite(lo) && is.finite(hi) && lo <= hi)
        sprintf("All intervals share [%.*f, %.*f]: the methods do not disagree\n  beyond their own uncertainty.\n",
                digits, lo, digits, hi)
        else "The intervals do not all overlap: the methods differ beyond their\n  own uncertainty.\n")
    }
  }
  invisible(x)
}

#' Summarize a method comparison
#'
#' @param object an object of class \code{"method_comparison"}.
#' @param ... ignored.
#' @return The summary data frame, invisibly.
#' @export
summary.method_comparison <- function(object, ...) {
  print(object$summary, row.names = FALSE)
  invisible(object$summary)
}
