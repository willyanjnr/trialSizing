## ============================================================================
## trialSize :: Paranaiba method (maximum curvature of the CV model)
## ----------------------------------------------------------------------------
## Paranaiba, Ferreira & Morais (2009). From the raw uniformity-trial grid, the
## first-order spatial autocorrelation (rho) is estimated along a serpentine
## walk; the optimal plot size and its CV follow in closed form:
##   Xo   = 10 * ( 2 (1 - rho^2) s^2 m )^(1/3) / m
##   CVxo = 100 * sqrt( (1 - rho^2) s^2 / m^2 ) / sqrt(Xo)
## with s^2 the variance and m the mean of the basic units.
##
## rho direction: the original method walks in the direction of the ROWS
## ("no sentido das linhas"), which is the default here. "col" and "mean" are
## offered for studies that walk down columns or average both directions.
## Validated against Cargnelutti Filho et al. (2014), black oat, Table 1.
## Uses the shared internal saver .save_lrp() defined in fit_lrp.R.
## ============================================================================

## Internal: serpentine first-order autocorrelation ---------------------------
.rho_serpentine <- function(mat, direction = "row") {
  err <- mat - mean(mat)
  ss  <- sum(err^2)
  if (ss == 0) return(NA_real_)
  walk <- c()
  if (direction == "row") {
    for (i in seq_len(nrow(mat)))
      walk <- c(walk, if (i %% 2 == 1) err[i, ] else rev(err[i, ]))
  } else {
    for (j in seq_len(ncol(mat)))
      walk <- c(walk, if (j %% 2 == 1) err[, j] else rev(err[, j]))
  }
  sum(walk[-1] * walk[-length(walk)]) / ss
}

## Internal: one trial --------------------------------------------------------
.paranaiba_one <- function(mat, rho_direction, trial_id) {

  m  <- mean(mat)
  s2 <- stats::var(as.vector(mat))
  cv <- 100 * sqrt(s2) / m

  if (s2 == 0 || m == 0) {
    warning(sprintf("Trial %s has zero variance or zero mean; results are NA.",
                    trial_id), call. = FALSE)
    return(data.frame(trial = trial_id, mean = m, variance = s2, CV = NA_real_,
                      rho_row = NA_real_, rho_col = NA_real_, rho = NA_real_,
                      Xo = NA_real_, CVxo = NA_real_, valid = FALSE,
                      stringsAsFactors = FALSE))
  }

  rho_row <- .rho_serpentine(mat, "row")
  rho_col <- .rho_serpentine(mat, "col")
  rho <- switch(rho_direction,
                row  = rho_row,
                col  = rho_col,
                mean = mean(c(rho_row, rho_col)))

  Xo   <- 10 * (2 * (1 - rho^2) * s2 * m)^(1 / 3) / m
  CVxo <- 100 * sqrt((1 - rho^2) * s2 / m^2) / sqrt(Xo)

  data.frame(trial = trial_id, mean = m, variance = s2, CV = cv,
             rho_row = rho_row, rho_col = rho_col, rho = rho,
             Xo = Xo, CVxo = CVxo, valid = TRUE, stringsAsFactors = FALSE)
}

## Public function ------------------------------------------------------------
#' Optimal plot size by the Paranaiba method
#'
#' Estimates the optimal plot size from the raw uniformity-trial grid using the
#' maximum curvature of the coefficient of variation model (Paranaiba, Ferreira
#' & Morais, 2009). Unlike [fit_lrp()], [fit_qrp()] and [fit_mcm()], which take
#' CV values already computed for several plot sizes, this method works directly
#' on the basic experimental units (BEU) and returns a closed-form estimate:
#' \deqn{X_o = \frac{10 \sqrt[3]{2 (1 - \rho^2) s^2 m}}{m}, \qquad
#'       CV_{Xo} = \frac{100 \sqrt{(1 - \rho^2) s^2 / m^2}}{\sqrt{X_o}}}
#' where \eqn{m} and \eqn{s^2} are the mean and variance of the BEU values and
#' \eqn{\rho} is the first-order spatial autocorrelation.
#'
#' @section Direction of the autocorrelation:
#' \eqn{\rho} is estimated along a serpentine walk through the grid. The original
#' method walks in the direction of the rows (\code{rho_direction = "row"},
#' the default). \code{"col"} walks down the columns and \code{"mean"} averages
#' both directions; these can give visibly different \eqn{\rho} and so a
#' different \eqn{X_o}.
#'
#' @section Input:
#' Supply either a matrix (one trial), a list of matrices (several trials), or a
#' data frame in long format with the value column plus row/column indices, or a
#' data frame whose \code{n_col} value columns hold the grid (see \code{value},
#' \code{row_id}, \code{col_id}, \code{trial}).
#'
#' @param .data a matrix, a list of matrices, or a data frame.
#' @param value name of the column holding the BEU measurement (long format).
#' @param row_id,col_id names of the row and column index columns (long format).
#' @param trial optional column name identifying the trial.
#' @param n_row,n_col grid dimensions; required only for a matrix supplied as a
#'   plain vector, or to check the grid of a long data frame.
#' @param rho_direction \code{"row"} (default), \code{"col"} or \code{"mean"}.
#'
#' @return An object of class \code{"paranaiba_fit"}: a list with \code{summary}
#'   (one row per trial: mean, variance, CV, rho_row, rho_col, rho, Xo, CVxo,
#'   valid) and \code{meta}.
#'
#' @references
#' Paranaiba, P. F., Ferreira, D. F. & Morais, A. R. (2009). Tamanho otimo de
#' parcelas experimentais: proposicao de metodos de estimacao. \emph{Revista
#' Brasileira de Biometria}, 27(2), 255-268. \cr
#' Cargnelutti Filho, A. et al. (2014). Tamanho de parcela e numero de
#' repeticoes em aveia preta. \emph{Ciencia Rural}, 44(10), 1732-1739.
#'
#' @seealso [fit_lrp()], [fit_mcm()], [calc_replicates()]
#' @examples
#' ## The packaged uniformity trial holds three replications, each a 6 x 8 grid
#' ## of basic experimental units (columns C1-C8).
#' rep1 <- as.matrix(dados_ensaio_C1[dados_ensaio_C1$Rep == 1, paste0("C", 1:8)])
#' dim(rep1)
#'
#' ## One trial, from a matrix
#' calc_paranaiba(rep1)
#'
#' ## Several trials, from a named list of matrices
#' grids <- lapply(split(dados_ensaio_C1, dados_ensaio_C1$Rep),
#'                 function(d) as.matrix(d[, paste0("C", 1:8)]))
#' names(grids) <- paste("Rep", names(grids))
#' par_fit <- calc_paranaiba(grids)
#' par_fit$summary
#'
#' ## rho is estimated along a serpentine walk. The original method walks the
#' ## rows; walking the columns instead can give a visibly different Xo.
#' cbind(
#'   row  = calc_paranaiba(grids)$summary$Xo,
#'   col  = calc_paranaiba(grids, rho_direction = "col")$summary$Xo,
#'   mean = calc_paranaiba(grids, rho_direction = "mean")$summary$Xo
#' )
#'
#' ## Long format: one row per basic unit, with row and column indices
#' long <- data.frame(
#'   rep   = rep(dados_ensaio_C1$Rep, times = 8),
#'   linha = rep(dados_ensaio_C1$Linha, times = 8),
#'   coluna = rep(paste0("C", 1:8), each = nrow(dados_ensaio_C1)),
#'   valor = unlist(dados_ensaio_C1[, paste0("C", 1:8)], use.names = FALSE)
#' )
#' calc_paranaiba(long, value = "valor", row_id = "linha", col_id = "coluna",
#'                trial = "rep")$summary
#'
#' \donttest{
#' plot(par_fit)
#'
#' ## CVxo of each trial feeds the number of replications
#' calc_replicates(treatments = c(5, 10, 20),
#'                 cv_percent = mean(par_fit$summary$CVxo),
#'                 lsd_percent = c(10, 20))
#' }
#' @export
calc_paranaiba <- function(.data, value = NULL, row_id = NULL, col_id = NULL,
                           trial = NULL, n_row = NULL, n_col = NULL,
                           rho_direction = c("row", "col", "mean")) {

  rho_direction <- match.arg(rho_direction)

  ## ---- build a named list of matrices (shared with calc_cv_shapes) ----
  mats <- .as_grid_list(.data, value, row_id, col_id, trial)

  ## optional dimension check
  if (!is.null(n_row) || !is.null(n_col)) {
    for (nm in names(mats)) {
      d <- dim(mats[[nm]])
      if ((!is.null(n_row) && d[1] != n_row) ||
          (!is.null(n_col) && d[2] != n_col))
        stop(sprintf("Trial '%s' is %d x %d, which does not match n_row/n_col.",
                     nm, d[1], d[2]), call. = FALSE)
    }
  }

  message(sprintf("Paranaiba method on %d trial(s); rho direction = '%s'.",
                  length(mats), rho_direction))

  summ <- do.call(rbind, Map(function(mm, nm)
    .paranaiba_one(mm, rho_direction, nm), mats, names(mats)))
  rownames(summ) <- NULL

  structure(list(summary = summ, matrices = mats,
                 meta = list(rho_direction = rho_direction,
                             n_trials = length(mats))),
            class = "paranaiba_fit")
}

## Methods --------------------------------------------------------------------
#' @export
print.paranaiba_fit <- function(x, ...) {
  cat("Paranaiba optimal plot size\n")
  cat("Trials:", x$meta$n_trials,
      " | rho direction:", x$meta$rho_direction,
      " | invalid:", sum(!x$summary$valid), "\n\n")
  d <- x$summary
  num <- c("mean", "variance", "CV", "rho_row", "rho_col", "rho", "Xo", "CVxo")
  d[num] <- lapply(d[num], round, 3)
  print(d, row.names = FALSE)
  if (x$meta$n_trials > 1)
    cat(sprintf("\nMean Xo: %.2f  |  Mean CVxo: %.2f\n",
                mean(x$summary$Xo, na.rm = TRUE),
                mean(x$summary$CVxo, na.rm = TRUE)))
  invisible(x)
}

#' @export
summary.paranaiba_fit <- function(object, ...) {
  d <- object$summary[object$summary$valid, ]
  cat("Optimal plot size (Xo) across trials:\n")
  print(summary(d$Xo))
  cat("\nCV at optimal plot size (CVxo):\n")
  print(summary(d$CVxo))
  invisible(object)
}

#' Plot Paranaiba estimates by trial
#'
#' Shows the optimal plot size per trial, with the mean across trials as a
#' dashed reference line.
#'
#' @param x a \code{"paranaiba_fit"} object.
#' @param y_var \code{"Xo"} (default) or \code{"CVxo"}.
#' @param title,xlab,ylab labels; \code{ylab = NULL} is chosen from \code{y_var}.
#' @param show_mean draw the across-trial mean as a dashed line.
#' @param fill_colour bar fill colour.
#' @param base_size,title_size,family,theme theme controls.
#' @param save,file,format,dpi,width,height,units,compression saving controls.
#' @param ... ignored.
#' @return A \code{ggplot} object (invisibly when saved).
#' @export
plot.paranaiba_fit <- function(x, y_var = c("Xo", "CVxo"),
                               title = "Paranaiba method",
                               xlab = "Trial", ylab = NULL, show_mean = TRUE,
                               fill_colour = "grey35",
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
    ylab <- if (y_var == "Xo") expression("Optimal plot size (" * m^2 * ")")
  else "CV at optimal plot size (%)"

  d <- x$summary
  d$.y <- d[[y_var]]

  g <- ggplot2::ggplot(d, ggplot2::aes(factor(trial), .y)) +
    ggplot2::geom_col(fill = fill_colour, width = 0.65, na.rm = TRUE)

  if (show_mean && sum(d$valid) > 1)
    g <- g + ggplot2::geom_hline(yintercept = mean(d$.y, na.rm = TRUE),
                                 linetype = 2, linewidth = 0.6)

  g <- g +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    (if (is.null(theme))
      .theme_trialsize(base_size = base_size, family = family)
     else theme) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = title_size))

  if (save) {
    .save_lrp(g, file, title, format, dpi, width, height, units, compression)
    return(invisible(g))
  }
  g
}
