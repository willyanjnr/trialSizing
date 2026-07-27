## ============================================================================
## trialSizing :: Optimal number of replications (Cargnelutti Filho et al., 2014)
## ----------------------------------------------------------------------------
## For each (treatments, LSD) the required replications solve, iteratively,
##   r = ( q_alpha * CV / LSD )^2,   q_alpha = studentized range (Tukey),
## where q_alpha depends on the error df, which depends on r. Two readings of the
## same calculation are reported:
##   r_continuous  the continuous fixed point (reproduces the article's tables);
##   r_optimal     the practical value, ceiling(r_continuous) floored at 2.
## Error df: t(r-1) for CRD, (t-1)(r-1) for RCBD, with t = number of treatments.
## Validated against Cargnelutti Filho et al. (2014), Tables 2 (CRD) and 3 (RCBD),
## CVxo = 9.25%. Uses the shared internal saver .save_lrp() from fit_lrp.R.
## ============================================================================

#' Estimate the optimal number of replications
#'
#' Estimates the number of replications for CRD or RCBD experiments following
#' Cargnelutti Filho et al. (2014), from the number of treatments, the
#' experimental coefficient of variation, the least significant difference (LSD,
#' as a percent of the mean) and the significance level, using the Tukey
#' (studentized range) critical value.
#'
#' @details
#' The required replications solve \eqn{r = (q_\alpha CV / LSD)^2}, where
#' \eqn{q_\alpha} is the Tukey critical value. Since \eqn{q_\alpha} depends on the
#' error df, which depends on \eqn{r}, the problem is solved iteratively. Two
#' readings are returned: \code{r_continuous}, the continuous fixed point (this
#' reproduces the published tables), and \code{r_optimal}, the practical
#' \code{ceiling(r_continuous)} floored at 2 (a design needs at least two
#' replications). Error df: \eqn{t(r-1)} for CRD and \eqn{(t-1)(r-1)} for RCBD.
#'
#' @param treatments numeric vector; numbers of treatments (>= 2).
#' @param cv_percent single numeric; experimental CV (%), e.g. the CVxo from
#'   [fit_lrp()]/[fit_qrp()]/[fit_mcm()].
#' @param lsd_percent numeric vector; least significant differences (% of mean).
#' @param alpha significance level for the Tukey test (default 0.05).
#' @param design "CRD" or "RCBD".
#' @param tol,max_iter convergence tolerance and iteration cap for the fixed
#'   point (defaults 1e-9 and 1000).
#'
#' @return An object of class \code{"replicates_fit"}: \code{data} (Treatments,
#'   CV_percent, LSD_percent, Alpha, Design, r_continuous, r_optimal, df_error,
#'   q_tukey, converged) and \code{meta}.
#'
#' @references
#' Cargnelutti Filho, A. et al. (2014). Tamanho de parcela e numero de
#' repeticoes em aveia preta. \emph{Ciencia Rural}, 44(10), 1732-1739.
#'
#' @examples
#' ## CV = 9.25% is the CVxo of the black oat trial (Cargnelutti Filho et al.,
#' ## 2014). How many replications to detect a difference of 10% or 20% of the
#' ## mean, for a few treatment numbers?
#' fit <- calc_replicates(treatments = c(5, 10, 20, 30), cv_percent = 9.25,
#'                        lsd_percent = c(10, 20), design = "CRD")
#' fit
#'
#' ## r_continuous is the fixed point the published tables report; r_optimal is
#' ## the practical ceiling, floored at 2.
#' fit$data[fit$data$LSD_percent == 10, ]
#'
#' ## A randomized complete block design spends one df per block, so it needs
#' ## slightly more replications than a completely randomized design.
#' rcbd <- calc_replicates(treatments = c(5, 10, 20, 30), cv_percent = 9.25,
#'                         lsd_percent = c(10, 20), design = "RCBD")
#' cbind(CRD = fit$data$r_optimal, RCBD = rcbd$data$r_optimal)
#'
#' \donttest{
#' ## The full published table: every treatment number from 3 to 50, at three
#' ## precision levels.
#' full <- calc_replicates(treatments = 3:50, cv_percent = 9.25,
#'                         lsd_percent = c(10, 20, 30), design = "CRD")
#' head(full$data)
#' plot(full)
#'
#' ## A stricter test costs replications
#' calc_replicates(treatments = 10, cv_percent = 9.25, lsd_percent = 10,
#'                 alpha = 0.01)$data[, c("Alpha", "r_continuous", "r_optimal")]
#' }
#' @export
calc_replicates <- function(treatments, cv_percent, lsd_percent, alpha = 0.05,
                            design = c("CRD", "RCBD"),
                            tol = 1e-9, max_iter = 1000) {

  design <- match.arg(design)

  if (any(treatments < 2))
    stop("`treatments` must be 2 or greater.", call. = FALSE)
  if (length(cv_percent) != 1 || cv_percent <= 0)
    stop("`cv_percent` must be a single positive number.", call. = FALSE)
  if (any(lsd_percent <= 0))
    stop("`lsd_percent` must be strictly positive.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be between 0 and 1.", call. = FALSE)

  df_of <- function(t, r) if (design == "CRD") t * (r - 1) else (t - 1) * (r - 1)

  ## continuous fixed point r = (q(df(r)) CV / LSD)^2
  ##
  ## The fixed point only exists while the implied error df stay positive,
  ## i.e. while r > 1. When the required precision is so loose that r would
  ## fall below the two replications any design needs, the iteration has no
  ## meaningful solution: it is reported at the r = 2 floor, flagged by
  ## `at_floor`, and r_optimal is 2.
  r_fixed <- function(t, lsd) {
    ## requirement evaluated at the floor: if two replications already
    ## satisfy it, the continuous fixed point lies below the floor
    q2 <- suppressWarnings(stats::qtukey(1 - alpha, nmeans = t,
                                         df = df_of(t, 2)))
    if (!is.nan(q2) && (q2 * cv_percent / lsd)^2 <= 2)
      return(list(r = 2, converged = TRUE, at_floor = TRUE))

    r <- 4; converged <- FALSE
    for (it in seq_len(max_iter)) {
      df <- df_of(t, r)
      if (df < 1) { r <- r + 0.5; next }
      q  <- suppressWarnings(stats::qtukey(1 - alpha, nmeans = t, df = df))
      if (is.nan(q)) { r <- r + 0.5; next }
      rn <- (q * cv_percent / lsd)^2
      if (abs(rn - r) < tol) { r <- rn; converged <- TRUE; break }
      r <- (r + rn) / 2
    }
    list(r = r, converged = converged, at_floor = FALSE)
  }

  one <- function(t, lsd) {
    fp <- r_fixed(t, lsd)
    r_opt <- max(2, ceiling(fp$r))
    df    <- df_of(t, r_opt)
    q     <- suppressWarnings(stats::qtukey(1 - alpha, nmeans = t, df = df))
    data.frame(Treatments = t, CV_percent = cv_percent, LSD_percent = lsd,
               Alpha = alpha, Design = design,
               r_continuous = fp$r, r_optimal = r_opt,
               df_error = df, q_tukey = q, converged = fp$converged,
               at_floor = fp$at_floor)
  }

  result_data <- do.call(rbind, lapply(lsd_percent, function(lsd)
    do.call(rbind, lapply(treatments, one, lsd = lsd))))
  rownames(result_data) <- NULL

  if (any(!result_data$converged))
    warning(sprintf("%d combination(s) did not converge within max_iter = %d.",
                    sum(!result_data$converged), max_iter), call. = FALSE)

  structure(list(data = result_data,
                 meta = list(cv_percent = cv_percent, alpha = alpha,
                             design = design)),
            class = "replicates_fit")
}

#' @export
print.replicates_fit <- function(x, ...) {
  cat("Optimal number of replications\n")
  cat("Design:", x$meta$design, " CV:", sprintf("%.2f%%", x$meta$cv_percent),
      " alpha:", x$meta$alpha, "\n")
  cat("Rows:", nrow(x$data), " | Non-converged:", sum(!x$data$converged), "\n\n")
  d <- x$data
  d$r_continuous <- round(d$r_continuous, 2); d$q_tukey <- round(d$q_tukey, 3)
  print(utils::head(d), row.names = FALSE)
  invisible(x)
}

#' @export
summary.replicates_fit <- function(object, ...) {
  cat("Optimal replications (integer) by LSD (%):\n")
  print(stats::aggregate(r_optimal ~ LSD_percent, data = object$data,
                         FUN = function(v) c(min = min(v), max = max(v))))
  invisible(object)
}

#' Plot the number of replications
#'
#' Draws replications against the number of treatments, one line per LSD level.
#'
#' @param x a \code{"replicates_fit"} object.
#' @param y_var which value to plot: \code{"r_optimal"} (integer, default) or
#'   \code{"r_continuous"} (the article's tabulated value).
#' @param title,xlab,ylab,colour_lab labels.
#' @param line_size,point_size line and point sizes.
#' @param base_size,title_size,family,theme theme controls.
#' @param save,file,format,dpi,width,height,units,compression saving controls.
#' @param ... ignored.
#' @return A \code{ggplot} object (invisibly when saved).
#' @export
plot.replicates_fit <- function(x, y_var = c("r_optimal", "r_continuous"),
                                title = "Number of replications",
                                xlab = "Number of treatments",
                                ylab = NULL, colour_lab = "LSD (%)",
                                line_size = 0.8, point_size = 1.8,
                                base_size = 12, title_size = NULL,
                                family = "sans", theme = NULL,
                                save = FALSE, file = NULL,
                                format = c("tiff", "png", "jpeg", "pdf", "eps"),
                                dpi = 300, width = 18, height = 12, units = "cm",
                                compression = "lzw", ...) {
  y_var  <- match.arg(y_var)
  format <- match.arg(format)
  if (is.null(title_size)) title_size <- base_size
  if (is.null(ylab))
    ylab <- if (y_var == "r_optimal") "Optimal number of replications"
  else "Replications (continuous)"

  d <- x$data
  d$LSD <- factor(d$LSD_percent)
  d$.y  <- d[[y_var]]

  g <- ggplot2::ggplot(d, ggplot2::aes(Treatments, .y, colour = LSD)) +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(title = title, x = xlab, y = ylab, colour = colour_lab) +
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
