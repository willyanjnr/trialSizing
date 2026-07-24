## ============================================================================
## trialSize :: CV by plot shape, from the raw uniformity-trial grid
## ----------------------------------------------------------------------------
## The step before any of the CV-based models. A uniformity trial is harvested
## as a grid of basic experimental units (BEU); the CV table that fit_lrp(),
## fit_qrp() and fit_mcm() consume is built from it by grouping adjacent BEU
## into every rectangular plot shape the grid admits, and computing the mean,
## standard deviation and CV of the plot totals for each shape.
## ============================================================================

## Shared grid builder (internal) ---------------------------------------------
## Accepts a matrix, a named list of matrices, or a long data frame, and returns
## a named list of matrices. Used by calc_cv_shapes() and calc_paranaiba() so
## both take exactly the same inputs.
.as_grid_list <- function(.data, value = NULL, row_id = NULL, col_id = NULL,
                          trial = NULL) {
  if (is.matrix(.data))
    return(list(`Trial 1` = .data))

  if (is.list(.data) && !is.data.frame(.data)) {
    if (!all(vapply(.data, is.matrix, logical(1))))
      stop("When `.data` is a list, every element must be a matrix.",
           call. = FALSE)
    nms <- names(.data)
    if (is.null(nms)) nms <- paste("Trial", seq_along(.data))
    return(stats::setNames(.data, nms))
  }

  if (is.data.frame(.data)) {
    if (is.null(value))
      stop("For a data frame, give `value` (the measurement column), plus ",
           "`row_id` and `col_id`.", call. = FALSE)
    need <- c(value, row_id, col_id, trial)
    missing_cols <- setdiff(need, names(.data))
    if (length(missing_cols))
      stop(sprintf("Column(s) not found in `.data`: %s.\n  Available: %s.",
                   paste(missing_cols, collapse = ", "),
                   paste(names(.data), collapse = ", ")), call. = FALSE)
    if (is.null(row_id) || is.null(col_id))
      stop("`row_id` and `col_id` are required for a long-format data frame.",
           call. = FALSE)

    split_by <- if (is.null(trial)) rep("Trial 1", nrow(.data))
    else as.character(.data[[trial]])
    groups <- split(.data, split_by)

    return(lapply(groups, function(g) {
      ri <- as.integer(factor(g[[row_id]]))
      ci <- as.integer(factor(g[[col_id]]))
      nr <- max(ri); nc <- max(ci)
      if (nrow(g) != nr * nc)
        stop(sprintf(paste0("Trial has %d rows but the grid implied by ",
                            "`row_id`/`col_id` is %d x %d = %d cells."),
                     nrow(g), nr, nc, nr * nc), call. = FALSE)
      mm <- matrix(NA_real_, nrow = nr, ncol = nc)
      mm[cbind(ri, ci)] <- as.numeric(g[[value]])
      if (anyNA(mm))
        stop("The grid has missing cells; every row/column combination must ",
             "be present.", call. = FALSE)
      mm
    }))
  }

  stop("`.data` must be a matrix, a list of matrices, or a data frame.",
       call. = FALSE)
}

## CV table for one grid (internal) -------------------------------------------
.cv_shapes_one <- function(m, trial_id, min_plots) {
  nr <- nrow(m); nc <- ncol(m)
  div <- function(k) seq_len(k)[k %% seq_len(k) == 0]

  out <- list()
  for (a in div(nr)) for (b in div(nc)) {
    n_plots <- (nr / a) * (nc / b)
    ## a single plot has no variance, and the article tables stop before it
    if (n_plots < min_plots) next

    ## total of every a x b block of adjacent basic units
    totals <- numeric(n_plots)
    k <- 0L
    for (i in seq_len(nr / a)) for (j in seq_len(nc / b)) {
      k <- k + 1L
      totals[k] <- sum(m[((i - 1) * a + 1):(i * a),
                         ((j - 1) * b + 1):(j * b)])
    }

    mu <- mean(totals); s <- stats::sd(totals)
    out[[length(out) + 1L]] <- data.frame(
      trial = trial_id, X_L = a, X_C = b, x = a * b, n = n_plots,
      mean = mu, sd = s, cv = 100 * s / mu,
      stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, out)
  res[order(res$x, res$X_L), , drop = FALSE]
}

## Public function ------------------------------------------------------------
#' CV by plot shape from a uniformity-trial grid
#'
#' Builds the coefficient-of-variation table that [fit_lrp()], [fit_qrp()] and
#' [fit_mcm()] consume, starting from the raw grid of basic experimental units
#' (BEU) as the trial was harvested. Adjacent BEU are grouped into every
#' rectangular plot shape the grid admits: for a grid of \eqn{L} rows by \eqn{C}
#' columns, each shape is \eqn{X_L \times X_C} basic units with \eqn{X_L} a
#' divisor of \eqn{L} and \eqn{X_C} a divisor of \eqn{C}. For each shape the
#' plot totals are formed and their mean, standard deviation and CV reported.
#'
#' This closes the pipeline: the grid goes in, the CV table comes out, and it
#' feeds the fitters directly, since the returned columns are already named
#' \code{x}, \code{cv} and \code{trial}.
#'
#' @section Shapes and plot counts:
#' A shape of area \eqn{X = X_L X_C} yields \eqn{n = LC/X} plots, so the CV of a
#' large plot size rests on few plots and is estimated with little information:
#' the last rows of the table are the least reliable. The \code{n} column is
#' returned for exactly this reason and is the natural weight for a weighted
#' fit. Shapes leaving fewer than \code{min_plots} plots are dropped, which by
#' default removes only the whole grid (a single plot, no variance).
#'
#' Shapes are not interchangeable at equal area: \eqn{1 \times 2} and
#' \eqn{2 \times 1} cover the same two basic units but run along different
#' directions of the field, and their CVs differ whenever fertility is not
#' isotropic. Both are reported, which is why the CV table has repeated
#' \code{x} values.
#'
#' @param .data a matrix (one trial), a named list of matrices (several
#'   trials), or a long data frame with one row per basic unit.
#' @param value name of the column holding the BEU measurement (long format).
#' @param row_id,col_id names of the row and column index columns (long format).
#' @param trial optional column name identifying the trial (long format).
#' @param min_plots minimum number of plots a shape must yield to be kept
#'   (default 2).
#'
#' @return A data frame with one row per shape and trial: \code{trial},
#'   \code{X_L}, \code{X_C} (shape in basic units), \code{x} (plot size),
#'   \code{n} (number of plots), \code{mean}, \code{sd} and \code{cv} (percent),
#'   ordered by plot size.
#'
#' @references
#' Cargnelutti Filho, A. et al. (2025). Determinacao do tamanho de parcela para
#' avaliar a massa de parte aerea de grao-de-bico. \emph{Revista Vivencias},
#' 21(43), 499-513. \cr
#' Paranaiba, P. F., Ferreira, D. F. & Morais, A. R. (2009). \emph{Revista
#' Brasileira de Biometria}, 27(2), 255-268.
#'
#' @seealso [fit_lrp()], [fit_qrp()], [fit_mcm()], [calc_paranaiba()]
#' @examples
#' ## Three replications of a uniformity trial, each a 6 x 8 grid of basic units
#' grids <- lapply(split(dados_ensaio_C1, dados_ensaio_C1$Rep),
#'                 function(d) as.matrix(d[, paste0("C", 1:8)]))
#' names(grids) <- paste("Rep", names(grids))
#'
#' tab <- calc_cv_shapes(grids)
#' head(tab, 8)
#'
#' ## Same area, different orientation, different CV
#' tab[tab$trial == "Rep 1" & tab$x == 2, ]
#'
#' ## The table feeds the fitters unchanged
#' fit_lrp(tab[tab$trial == "Rep 1", ], x = "x", cv = "cv", step = 0.05)
#'
#' \donttest{
#' ## One model per replication, straight from the grids
#' fit_lrp(tab, x = "x", cv = "cv", trial = "trial", step = 0.01)$summary
#'
#' ## n falls as the plot grows: the last shapes rest on very few plots
#' unique(tab[, c("x", "n")])
#' }
#' @export
calc_cv_shapes <- function(.data, value = NULL, row_id = NULL, col_id = NULL,
                           trial = NULL, min_plots = 2) {

  if (!is.numeric(min_plots) || length(min_plots) != 1 || is.na(min_plots) ||
      min_plots < 2)
    stop("`min_plots` must be a single number of at least 2.", call. = FALSE)

  mats <- .as_grid_list(.data, value, row_id, col_id, trial)

  for (nm in names(mats)) {
    if (!is.numeric(mats[[nm]]))
      stop(sprintf("Trial '%s' is not a numeric matrix.", nm), call. = FALSE)
    if (anyNA(mats[[nm]]))
      stop(sprintf("Trial '%s' has missing values.", nm), call. = FALSE)
  }

  message(sprintf("CV by shape for %d trial(s).", length(mats)))

  res <- do.call(rbind, Map(function(mm, nm)
    .cv_shapes_one(mm, nm, min_plots), mats, names(mats)))
  rownames(res) <- NULL
  res
}
