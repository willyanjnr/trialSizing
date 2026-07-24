## ============================================================================
## trialSize :: field map of a uniformity trial
## ----------------------------------------------------------------------------
## The surface is ordinary kriging with the variogram fitted by check_trial();
## the points are the basic units themselves, filled on the same scale, so the
## reader sees the data and the interpolation in one figure and can tell them
## apart. Styling aims at a journal figure: no decoration that does not carry
## information, square cells, one colourbar.
## ============================================================================

## Theme (internal) -----------------------------------------------------------
## Kept separate so the other plot methods can adopt it later.
.theme_trialsize <- function(base_size = 12, family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = family) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "grey20",
                                               linewidth = 0.4),
      axis.ticks       = ggplot2::element_line(colour = "grey20",
                                               linewidth = 0.4),
      axis.ticks.length = ggplot2::unit(3, "pt"),
      axis.text        = ggplot2::element_text(colour = "grey20",
                                               size = base_size * 0.85),
      axis.title       = ggplot2::element_text(colour = "grey10"),
      plot.title       = ggplot2::element_text(colour = "grey10",
                                               size = base_size * 1.1,
                                               hjust = 0),
      plot.subtitle    = ggplot2::element_text(colour = "grey40",
                                               size = base_size * 0.85,
                                               hjust = 0),
      plot.caption     = ggplot2::element_text(colour = "grey45",
                                               size = base_size * 0.75,
                                               hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.title     = ggplot2::element_text(colour = "grey10",
                                               size = base_size * 0.85),
      legend.text      = ggplot2::element_text(colour = "grey30",
                                               size = base_size * 0.75),
      legend.background = ggplot2::element_blank(),
      legend.key       = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(colour = "grey10",
                                               size = base_size * 0.9,
                                               hjust = 0),
      plot.margin      = ggplot2::margin(8, 8, 8, 8))
}

## Axis breaks (internal) -----------------------------------------------------
## The axes index rows and columns, so fractional breaks are meaningless. When
## the units are one apart, only whole numbers are kept.
.grid_breaks <- function(spacing) {
  force(spacing)
  function(lims) {
    b <- pretty(lims, n = 6)
    b <- b[b >= lims[1] - 1e-9 & b <= lims[2] + 1e-9]
    if (isTRUE(all.equal(unname(spacing), 1)))
      b <- b[abs(b - round(b)) < 1e-9]
    b
  }
}

## Colour scale (internal) ----------------------------------------------------
.fill_scale <- function(palette, name, aes_type = "fill") {
  fun <- if (aes_type == "fill") ggplot2::scale_fill_gradientn
  else ggplot2::scale_colour_gradientn
  guide <- ggplot2::guide_colourbar(
    barwidth = ggplot2::unit(0.5, "lines"),
    barheight = ggplot2::unit(7, "lines"),
    ticks.colour = "white", frame.colour = NA)

  if (palette == "viridis") {
    if (aes_type == "fill")
      return(ggplot2::scale_fill_viridis_c(name = name, guide = guide))
    return(ggplot2::scale_colour_viridis_c(name = name, guide = guide))
  }
  cols <- switch(palette,
    blues = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    greys = c("#FFFFFF", "#D9D9D9", "#969696", "#525252", "#000000"),
    terrain = c("#FFFFCC", "#C7E9B4", "#7FCDBB", "#41B6C4", "#225EA8"),
    stop("Unknown palette.", call. = FALSE))
  fun(name = name, colours = cols, guide = guide)
}

#' Field map of a checked uniformity trial
#'
#' Draws the trial as a map: an ordinary-kriging surface built from the
#' variogram fitted by [check_trial()], with the basic experimental units drawn
#' on top and filled on the same colour scale. The surface shows where the
#' field is systematically better or worse; the points show the data the
#' surface came from, so an interpolation artefact cannot be mistaken for a
#' measurement.
#'
#' With several trials the maps are faceted and share one colour scale, which
#' is what makes them comparable.
#'
#' @param x an object of class \code{"trial_check"}.
#' @param resolution number of interpolation cells along the longer side of the
#'   field (default 120). The surface costs one linear solve, so raising this
#'   is cheap.
#' @param points logical; draw the basic units (default \code{TRUE}).
#' @param point_values logical; fill the points with their own value
#'   (default \code{TRUE}). \code{FALSE} draws them as plain markers, showing
#'   only the sampling positions.
#' @param point_size,point_stroke,point_colour size of the unit markers, the
#'   width of their outline, and its colour. A dark rim keeps the markers
#'   visible over the pale end of any palette.
#' @param surface logical; draw the kriged surface (default \code{TRUE}). With
#'   \code{FALSE} only the units are drawn, which is the honest picture when the
#'   variogram shows weak spatial dependence.
#' @param palette one of \code{"viridis"} (default), \code{"blues"},
#'   \code{"greys"} or \code{"terrain"}. The default is perceptually uniform and
#'   readable in greyscale and to colour-blind readers, which several journals
#'   now require; \code{"blues"} gives the classic look.
#' @param title,subtitle,caption,legend_title plot labels. \code{NULL} leaves a
#'   sensible default; \code{NA} removes the element.
#' @param xlab,ylab axis titles.
#' @param base_size,family base font size and family.
#' @param ... ignored.
#'
#' @return A \code{ggplot} object.
#' @seealso [check_trial()]
#' @examples
#' grid1 <- as.matrix(dados_ensaio_C1[dados_ensaio_C1$Rep == 1,
#'                                    paste0("C", 1:8)])
#' chk <- check_trial(grid1)
#'
#' \donttest{
#' plot(chk)
#'
#' ## the classic look, and points as plain position markers
#' plot(chk, palette = "blues", point_values = FALSE)
#'
#' ## data only, no interpolation
#' plot(chk, surface = FALSE)
#' }
#' @export
plot.trial_check <- function(x, resolution = 120, points = TRUE,
                             point_values = TRUE, point_size = 2.4,
                             point_stroke = 0.4, point_colour = "grey20",
                             surface = TRUE,
                             palette = c("viridis", "blues", "greys", "terrain"),
                             title = NULL, subtitle = NULL, caption = NULL,
                             legend_title = NULL, xlab = "Row", ylab = "Column",
                             base_size = 12, family = "sans", ...) {

  palette <- match.arg(palette)
  if (!is.numeric(resolution) || length(resolution) != 1 || resolution < 10)
    stop("`resolution` must be a single number of at least 10.", call. = FALSE)

  surf <- list(); pts <- list()

  for (k in x$checks) {
    m  <- k$matrix
    xy <- k$coords
    v  <- as.numeric(m)
    ok <- !is.na(v)

    pts[[length(pts) + 1L]] <- data.frame(
      trial = k$trial, xpos = xy[ok, 1], ypos = xy[ok, 2], value = v[ok],
      stringsAsFactors = FALSE)

    if (isTRUE(surface) && !is.null(k$variogram) && sum(ok) >= 4) {
      ## pad by half a unit spacing so the surface reaches the panel edge and
      ## the outermost units are not sliced in half by the frame
      pad <- x$meta$cell_size / 2
      rx <- range(xy[, 1]) + c(-1, 1) * pad[1]
      ry <- range(xy[, 2]) + c(-1, 1) * pad[2]
      ## keep the cells square: the denser axis gets `resolution`
      span <- c(diff(rx), diff(ry))
      nx <- max(10, round(resolution * span[1] / max(span)))
      ny <- max(10, round(resolution * span[2] / max(span)))
      gx <- seq(rx[1], rx[2], length.out = nx)
      gy <- seq(ry[1], ry[2], length.out = ny)
      new <- as.matrix(expand.grid(x = gx, y = gy))

      z <- .krige(v[ok], xy[ok, , drop = FALSE], new, k$variogram)
      surf[[length(surf) + 1L]] <- data.frame(
        trial = k$trial, xpos = new[, 1], ypos = new[, 2], value = z,
        stringsAsFactors = FALSE)
    }
  }

  pts  <- do.call(rbind, pts)
  surf <- if (length(surf)) do.call(rbind, surf) else NULL

  if (is.null(legend_title)) legend_title <- "Value"
  if (is.null(title)) title <- "Uniformity trial"
  if (is.null(subtitle) && length(x$checks) == 1) {
    v <- x$checks[[1]]$variogram
    subtitle <- if (is.null(v))
      sprintf("%d x %d basic units", x$checks[[1]]$dim[["rows"]],
              x$checks[[1]]$dim[["cols"]])
    else sprintf("%d x %d basic units | %s variogram, range %.2f, %s spatial dependence",
                 x$checks[[1]]$dim[["rows"]], x$checks[[1]]$dim[["cols"]],
                 v$model, v$range, v$dependence)
  }
  if (is.null(caption) && !is.null(surf))
    caption <- "Surface: ordinary kriging. Points: basic experimental units."

  g <- ggplot2::ggplot()

  if (!is.null(surf))
    g <- g + ggplot2::geom_raster(
      data = surf, ggplot2::aes(xpos, ypos, fill = value), interpolate = TRUE)

  if (isTRUE(points)) {
    g <- if (isTRUE(point_values))
      ## a dark rim reads against every palette; a white one disappears
      ## wherever the scale itself is pale
      g + ggplot2::geom_point(
        data = pts, ggplot2::aes(xpos, ypos, fill = value),
        shape = 21, size = point_size, stroke = point_stroke,
        colour = point_colour)
    else
      g + ggplot2::geom_point(
        data = pts, ggplot2::aes(xpos, ypos),
        size = point_size, colour = point_colour)
  }

  g <- g + .fill_scale(palette, legend_title, "fill")

  if (length(x$checks) > 1)
    g <- g + ggplot2::facet_wrap(~ trial)

  cs <- x$meta$cell_size
  g +
    ggplot2::scale_x_continuous(breaks = .grid_breaks(cs[1])) +
    ggplot2::scale_y_continuous(breaks = .grid_breaks(cs[2])) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(title = if (is.na(title)) NULL else title,
                  subtitle = if (length(subtitle) && is.na(subtitle)) NULL else subtitle,
                  caption = if (length(caption) && is.na(caption)) NULL else caption,
                  x = xlab, y = ylab) +
    .theme_trialsize(base_size = base_size, family = family)
}
