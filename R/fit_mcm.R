## ============================================================================
## trialSize :: Modified Maximum Curvature (MCM) plot-size model
## ----------------------------------------------------------------------------
## Meier & Lessman (1971), after Lessman & Atkins (1963). The CV vs plot-size
## relationship CV = a * X^(-b) is fit and the optimum plot size is the point of
## maximum curvature of that curve:
##   Xc = [ a^2 b^2 (2b + 1) / (b + 2) ] ^ ( 1 / (2b + 2) ).
##
## Estimation (method):
##   "nls"       nonlinear least squares of CV = a X^(-b) on the original scale,
##               seeded from the log-log fit. Reproduces Cargnelutti Filho et al.
##               (2025) and the modern Brazilian plot-size papers. DEFAULT.
##   "loglinear" linear regression of log(CV) on log(X) (the classic Meier &
##               Lessman route). Pass `df` for the Federer (1955) weighting.
## Degrees-of-freedom correction (Federer, 1955): pass `df`, the degrees of
## freedom per point, to weight the fit. Some authors apply it, others do not;
## `df = NULL` (default) is the unweighted fit.
##
## The cost-factor modification (K1, K2) of the classic method is NOT applied
## here; Xc is reported in the units of `x` (basic units in the original method).
## Uses the shared internal saver .save_lrp() defined in fit_lrp.R.
## ============================================================================

## Bootstrap of the optimum (internal) ----------------------------------------
## Only half of what the plateau models get. The interval carries over: resample
## the shapes, refit, take the quantiles of Xc. The test of existence does not:
## there is no breakpoint whose presence is in doubt, since the curve
## CV = a X^(-b) has a maximum-curvature point whenever b > 0. What can be
## questioned is b itself, so the interval for b is reported alongside; if it
## covers zero, the CV does not demonstrably fall with plot size and Xc means
## nothing.
.mcm_boot <- function(x, cv, method, df, n_boot, conf_level) {
  n <- length(x)
  xc <- rep(NA_real_, n_boot)
  bb <- rep(NA_real_, n_boot)

  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(x[idx])) < 3) next
    f <- tryCatch(
      suppressWarnings(.mcm_fit_one(x[idx], cv[idx], method, df[idx])),
      error = function(e) NULL)
    if (is.null(f)) next
    xc[i] <- unname(f$parameters["Breakpoint"])
    bb[i] <- unname(f$coefficients["b"])
  }
  ok <- !is.na(xc)

  alpha <- (1 - conf_level) / 2
  qs <- function(v) if (sum(ok) >= 2)
    unname(stats::quantile(v[ok], c(alpha, 1 - alpha))) else c(NA, NA)

  list(
    ci         = qs(xc),
    se         = if (sum(ok) >= 2) stats::sd(xc[ok]) else NA_real_,
    ci_b       = qs(bb),
    replicates = xc[ok],
    n_valid    = sum(ok),
    n_boot     = n_boot,
    conf_level = conf_level
  )
}

## Internal single-series engine ---------------------------------------------
.mcm_fit_one <- function(x, cv, method = "nls", df = NULL, bootstrap = FALSE,
                         n_boot = 1000, conf_level = 0.95) {

  if (!is.numeric(x) || !is.numeric(cv))
    stop("`x` and `cv` must be numeric.", call. = FALSE)
  if (length(x) != length(cv))
    stop("`x` and `cv` must have the same length.", call. = FALSE)
  if (anyNA(x) || anyNA(cv))
    stop("`x` and `cv` cannot contain missing values (NA).", call. = FALSE)
  if (any(x <= 0) || any(cv <= 0))
    stop("`x` and `cv` must be strictly positive.", call. = FALSE)
  if (length(x) < 3 || length(unique(x)) < 3)
    stop("At least 3 distinct `x` values are required.", call. = FALSE)
  if (!is.null(df)) {
    if (!is.numeric(df) || length(df) != length(x) || anyNA(df) || any(df <= 0))
      stop("`df` must be a positive numeric vector matching `x` in length.",
           call. = FALSE)
  }
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

  ## log-log fit (also the starting values for nls)
  ll <- if (is.null(df)) stats::lm(log(cv) ~ log(x))
  else             stats::lm(log(cv) ~ log(x), weights = df)
  a0 <- unname(exp(stats::coef(ll)[1]))
  b0 <- unname(-stats::coef(ll)[2])
  used <- method

  if (method == "loglinear") {
    a <- a0; b <- b0
  } else {                                   # method == "nls"
    nlsfit <- tryCatch(
      if (is.null(df))
        stats::nls(cv ~ a * x^(-b), start = list(a = a0, b = b0))
      else
        stats::nls(cv ~ a * x^(-b), start = list(a = a0, b = b0), weights = df),
      error = function(e) NULL)
    if (is.null(nlsfit)) {
      warning("nls did not converge; using the log-linear estimate instead.",
              call. = FALSE)
      a <- a0; b <- b0; used <- "loglinear (nls fallback)"
    } else {
      cf <- stats::coef(nlsfit); a <- unname(cf["a"]); b <- unname(cf["b"])
    }
  }

  if (b <= 0)
    warning("Estimated exponent b <= 0; CV does not decrease with plot size, ",
            "so the maximum-curvature point is not meaningful.", call. = FALSE)

  ## point of maximum curvature (in the units of x)
  breakpoint          <- ((a^2 * b^2 * (2 * b + 1)) / (b + 2))^(1 / (2 * b + 2))
  breakpoint_response <- a * breakpoint^(-b)

  pred      <- a * x^(-b)
  residuals <- cv - pred
  r2   <- 1 - sum(residuals^2) / sum((cv - mean(cv))^2)
  rmse <- sqrt(mean(residuals^2))

  boot <- NULL
  if (isTRUE(bootstrap))
    boot <- .mcm_boot(x, cv, method, df, n_boot, conf_level)

  structure(
    list(
      coefficients = c(a = a, b = b),
      parameters   = c(Breakpoint = breakpoint,
                       Breakpoint_Response = breakpoint_response,
                       R2 = r2, RMSE = rmse),
      fitted = pred, residuals = residuals,
      data = data.frame(x = x, cv = cv),
      method = used, weighted = !is.null(df), bootstrap = boot
    ),
    class = "mcm_fit"
  )
}

## Public fitter --------------------------------------------------------------
#' Fit the Modified Maximum Curvature (MCM) plot-size model
#'
#' Estimates the optimal plot size by the modified maximum curvature method of
#' Meier & Lessman (1971). The relationship \eqn{CV = a\,X^{-b}} is fit by linear
#' regression of \eqn{\log CV} on \eqn{\log X}, and the optimum is the point of
#' maximum curvature
#' \deqn{X_c = \left[ a^2 b^2 (2b + 1) / (b + 2) \right]^{1/(2b + 2)}.}
#'
#' @section Estimation method:
#' \code{method = "nls"} (default) fits \eqn{CV = a X^{-b}} by nonlinear least
#' squares on the original scale (seeded from the log-log fit), reproducing
#' Cargnelutti Filho et al. (2025) and the modern Brazilian plot-size papers.
#' \code{method = "loglinear"} fits the line \eqn{\log CV = \log a - b \log X}
#' (the classic Meier & Lessman route). If \code{"nls"} fails to converge it
#' falls back to the log-linear estimate with a warning.
#'
#' @section Degrees-of-freedom correction:
#' Pass \code{df}, the degrees of freedom of each point, to weight the fit as
#' proposed by Federer (1955); this down-weights the CV of larger plot sizes,
#' estimated from fewer plots. Some authors apply it, others do not;
#' \code{df = NULL} (default) is unweighted. The cost-factor modification
#' (K1, K2) of the classic method is not applied; \eqn{X_c} is in the units of
#' \code{x}.
#'
#' @section Uncertainty of the optimum:
#' \code{bootstrap = TRUE} resamples the shapes with replacement, refits, and
#' returns a percentile interval for \eqn{X_c} in \code{$bootstrap$ci}.
#'
#' Unlike [fit_lrp()] and [fit_qrp()], there is no test for the existence of the
#' optimum, because there is no breakpoint whose presence is in doubt: the curve
#' \eqn{CV = a X^{-b}} has a point of maximum curvature whenever \eqn{b > 0}.
#' What can be questioned is \eqn{b}. Its bootstrap interval is therefore
#' reported as \code{$bootstrap$ci_b}; if that interval covers zero, the CV does
#' not demonstrably fall with plot size and \eqn{X_c} carries no meaning.
#'
#' Set the random seed before calling to make the result reproducible.
#'
#' @section Two ways to call:
#' \describe{
#'   \item{Vectors}{\code{fit_mcm(x, cv, df = NULL)} returns an \code{"mcm_fit"}.}
#'   \item{Data frame}{\code{fit_mcm(.data, x = "x", cv = "cv", df = "df",
#'     trial = "trial")}; with \code{trial}, one model per trial is fit and an
#'     \code{"mcm_multi"} object is returned.}
#' }
#'
#' @param .data optional data frame; when supplied the other arguments are
#'   column names.
#' @param x,cv numeric vectors (plot size and CV in percent), or column names.
#' @param df optional degrees of freedom per point for the Federer (1955)
#'   weighting, as a vector or a column name. \code{NULL} = unweighted.
#' @param trial optional column name identifying the trial.
#' @param method estimation method: \code{"nls"} (default, original scale;
#'   reproduces the chickpea article) or \code{"loglinear"} (Meier & Lessman
#'   regression).
#' @param bootstrap logical; if \code{TRUE}, also estimate the uncertainty of
#'   \eqn{X_c} by resampling the shapes (default \code{FALSE}). See the section
#'   "Uncertainty of the optimum".
#' @param n_boot number of bootstrap resamples (default 1000), used only when
#'   \code{bootstrap = TRUE}.
#' @param conf_level confidence level of the percentile interval (default 0.95),
#'   used only when \code{bootstrap = TRUE}.
#' @return An \code{"mcm_fit"} (single series) or \code{"mcm_multi"} (per trial).
#'   With \code{bootstrap = TRUE} the fit also carries \code{bootstrap}, a list
#'   with \code{ci}, \code{se}, \code{ci_b}, \code{replicates}, \code{n_valid}
#'   and \code{conf_level}.
#' @references
#' Meier, V. D. & Lessman, K. J. (1971). Estimation of optimum field plot shape
#' and size for testing yield in \emph{Crambe abyssinica} Hochst.
#' \emph{Crop Science}, 11, 648-650. \cr
#' Federer, W. T. (1955). \emph{Experimental Design}. Macmillan, New York.
#' @seealso [fit_lrp()], [fit_qrp()], [plot.mcm_fit()]
#' @examples
#' ## Chickpea uniformity trial, trial 1 (Cargnelutti Filho et al., 2025)
#' X   <- c(1, 2, 2, 3, 3, 4, 6, 6, 6, 6, 9, 12, 12, 18, 18)
#' CV1 <- c(30.40, 19.51, 23.72, 12.89, 21.32, 16.69, 6.71, 10.75,
#'          17.58, 14.94, 11.93, 3.18, 8.63, 4.25, 11.41)
#'
#' fit <- fit_mcm(X, CV1)
#' fit
#'
#' ## The fitted decay CV = a * X^-b, and the CV expected at unobserved sizes
#' coef(fit)
#' predict(fit, newx = c(2, 5, 7.5, 15))
#'
#' ## "loglinear" is the classic Meier & Lessman route (regression of log CV on
#' ## log X); "nls" (default) fits on the original scale and is what the modern
#' ## plot-size articles report. They rarely agree exactly.
#' c(nls       = unname(fit$parameters["Breakpoint"]),
#'   loglinear = unname(fit_mcm(X, CV1, method = "loglinear")$parameters["Breakpoint"]))
#'
#' ## Federer (1955) weighting: the CV of a large plot size rests on fewer
#' ## plots, so pass the degrees of freedom of each point to down-weight it.
#' n_plots <- 72 / X
#' fit_mcm(X, CV1, df = n_plots - 1)$parameters["Breakpoint"]
#'
#' ## Uncertainty of Xc, off by default. There is no existence test here: the
#' ## curve always has a maximum-curvature point when b > 0, so what gets an
#' ## interval is b itself.
#' set.seed(1)
#' unc <- fit_mcm(X, CV1, bootstrap = TRUE, n_boot = 500)
#' unc
#' unc$bootstrap$ci_b
#'
#' ## One model per trial
#' trials <- rbind(
#'   data.frame(x = X, cv = CV1,        trial = "T1"),
#'   data.frame(x = X, cv = CV1 * 0.85, trial = "T2")
#' )
#' fit_mcm(trials, x = "x", cv = "cv", trial = "trial")$summary
#'
#' \donttest{
#' plot(fit, title = "Chickpea, trial 1")
#'
#' ## The three CV-based methods on the same data. MCM is the most
#' ## conservative and QRP the most generous; the ordering is systematic.
#' c(MCM = unname(fit$parameters["Breakpoint"]),
#'   LRP = unname(fit_lrp(X, CV1, step = 0.01)$parameters["Breakpoint"]),
#'   QRP = unname(fit_qrp(X, CV1, step = 0.01)$parameters["Breakpoint"]))
#' }
#' @export
fit_mcm <- function(.data = NULL, x = NULL, cv = NULL, df = NULL, trial = NULL,
                    method = c("nls", "loglinear"), bootstrap = FALSE,
                    n_boot = 1000, conf_level = 0.95) {

  method <- match.arg(method)

  if (!is.null(.data) && !is.data.frame(.data)) {
    cv <- x; x <- .data; .data <- NULL
  }

  if (is.null(.data)) {
    if (is.null(x) || is.null(cv))
      stop("Provide numeric `x` and `cv`, or a data frame as `.data`.",
           call. = FALSE)
    return(.mcm_fit_one(x, cv, method, df, bootstrap, n_boot, conf_level))
  }

  if (!is.data.frame(.data)) stop("`.data` must be a data frame.", call. = FALSE)

  xcol  <- if (is.null(x))  "x"  else x
  cvcol <- if (is.null(cv)) "cv" else cv
  dfcol <- df  # column name or NULL
  need  <- c(xcol, cvcol, if (!is.null(dfcol)) dfcol, if (!is.null(trial)) trial)
  missing_cols <- setdiff(need, names(.data))
  if (length(missing_cols))
    stop(sprintf("Column(s) not found in `.data`: %s.\n  Available columns: %s.",
                 paste(missing_cols, collapse = ", "),
                 paste(names(.data), collapse = ", ")), call. = FALSE)

  get_df <- function(g) if (is.null(dfcol)) NULL else g[[dfcol]]

  if (is.null(trial)) {
    message(sprintf("Using x = '%s', cv = '%s'%s (single series).",
                    xcol, cvcol,
                    if (is.null(dfcol)) "" else sprintf(", df = '%s'", dfcol)))
    return(.mcm_fit_one(.data[[xcol]], .data[[cvcol]], method, get_df(.data),
                        bootstrap, n_boot, conf_level))
  }

  groups <- split(.data, .data[[trial]])
  message(sprintf("Using x = '%s', cv = '%s'%s, trial = '%s' -> %d trials.",
                  xcol, cvcol,
                  if (is.null(dfcol)) "" else sprintf(", df = '%s'", dfcol),
                  trial, length(groups)))

  fits <- lapply(groups, function(g) .mcm_fit_one(g[[xcol]], g[[cvcol]], method,
                                                  get_df(g), bootstrap, n_boot,
                                                  conf_level))
  summ <- do.call(rbind, Map(function(f, nm) data.frame(
    trial      = nm,
    a          = round(unname(f$coefficients["a"]), 4),
    b          = round(unname(f$coefficients["b"]), 4),
    breakpoint = round(unname(f$parameters["Breakpoint"]), 4),
    plateau    = round(unname(f$parameters["Breakpoint_Response"]), 4),
    R2         = round(unname(f$parameters["R2"]), 4),
    RMSE       = round(unname(f$parameters["RMSE"]), 4),
    method     = f$method,
    weighted   = f$weighted,
    stringsAsFactors = FALSE), fits, names(fits)))
  rownames(summ) <- NULL

  if (isTRUE(bootstrap)) {
    summ$Xc_lwr <- round(vapply(fits, function(f) f$bootstrap$ci[1], 0), 4)
    summ$Xc_upr <- round(vapply(fits, function(f) f$bootstrap$ci[2], 0), 4)
  }

  structure(list(fits = fits, summary = summ), class = "mcm_multi")
}

## Methods --------------------------------------------------------------------
#' Predictions from an MCM fit
#' @param object an \code{"mcm_fit"} object.
#' @param newx numeric predictor values; defaults to the fitted data.
#' @param ... ignored.
#' @return numeric vector of predicted CV values.
#' @export
predict.mcm_fit <- function(object, newx = NULL, ...) {
  a <- object$coefficients["a"]; b <- object$coefficients["b"]
  if (is.null(newx)) newx <- object$data$x
  unname(a * newx^(-b))
}

#' @export
print.mcm_fit <- function(x, ...) {
  cat("Modified Maximum Curvature (MCM) fit\n")
  cat("Method:                 ", x$method, "\n")
  cat("Weighted by df:         ", x$weighted, "\n")
  cat("Breakpoint (Xo):        ", sprintf("%.3f", x$parameters["Breakpoint"]), "\n")
  if (!is.null(x$bootstrap)) {
    b <- x$bootstrap
    cat(sprintf("  %.0f%% CI (percentile):  [%.3f, %.3f]   SE %.3f\n",
                100 * b$conf_level, b$ci[1], b$ci[2], b$se))
    cat(sprintf("  exponent b:            [%.3f, %.3f]  (%d resamples)\n",
                b$ci_b[1], b$ci_b[2], b$n_valid))
  }
  cat("CV at breakpoint:       ", sprintf("%.3f", x$parameters["Breakpoint_Response"]), "\n")
  cat("R2:", sprintf("%.3f", x$parameters["R2"]),
      " RMSE:", sprintf("%.3f", x$parameters["RMSE"]), "\n")
  invisible(x)
}

#' @export
summary.mcm_fit <- function(object, ...) {
  cat("Model coefficients (CV = a * X^(-b)):\n"); print(object$coefficients)
  cat("\nGoodness of fit:\n"); print(object$parameters)
  invisible(object)
}

#' Plot an MCM fit (publication style)
#'
#' Draws the observed points, the fitted power curve \eqn{a X^{-b}}, dotted
#' guides to the maximum-curvature point and the model annotation. Unlike the
#' plateau models there is no flat segment: the curve keeps decreasing.
#'
#' @param x an \code{"mcm_fit"} object.
#' @param title,annotate_model,xlab,ylab,decimal_mark,digits_coef,digits_stat
#'   as in [plot.lrp_fit()].
#' @param point_size,line_size,point_colour,line_colour,bp_colour aesthetics.
#' @param base_size,label_size,title_size,family,theme theme and text controls.
#' @param save,file,format,dpi,width,height,units,compression saving controls.
#' @param ... ignored.
#' @return A \code{ggplot} object (invisibly when saved).
#' @export
plot.mcm_fit <- function(x, title = "Modified Maximum Curvature",
                         annotate_model = TRUE, xlab = NULL, ylab = "CV (%)",
                         decimal_mark = c(".", ","),
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
  curve <- data.frame(x = seq(min(d$x), max(d$x), length.out = 400))
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
    a_s  <- fmtn(a, digits_coef);  b_s <- fmtn(b, digits_coef)
    xo_s <- fmtn(xo, digits_stat); p_s <- fmtn(p, digits_stat)
    r2_s <- fmtn(r2, digits_stat)

    eq1 <- paste0("CV[(x)]=='", a_s, "'*X^{-'", b_s, "'}")
    eq2 <- paste0("R^2=='", r2_s, "'")
    lXo <- paste0("X[o]=='", xo_s, "'")
    lCV <- paste0("CV[Xo]=='", p_s, "'")

    xm <- 0.55 * xmax
    ym <- ymax * c(0.97, 0.88)
    for (i in 1:2)
      g <- g + ggplot2::annotate("text", x = xm, y = ym[i],
                                 label = c(eq1, eq2)[i], parse = TRUE,
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

#' @export
print.mcm_multi <- function(x, ...) {
  cat(sprintf("MCM fits for %d trials\n\n", length(x$fits)))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
summary.mcm_multi <- function(object, ...) {
  print(object$summary, row.names = FALSE)
  invisible(object)
}

#' Plot MCM fits for several trials (paginated grid)
#'
#' @param x an \code{"mcm_multi"} object.
#' @param save,file,format,dpi,width,height,units,compression saving controls.
#' @param ... styling arguments forwarded to [plot.mcm_fit()].
#' @return A list of page objects, invisibly.
#' @export
plot.mcm_multi <- function(x, save = FALSE, file = NULL,
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
    base <- if (is.null(file)) "mcm_trials" else sub("\\.[^.]*$", "", file)
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
