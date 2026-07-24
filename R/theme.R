## ============================================================================
## trialSize :: shared plot theme
## ----------------------------------------------------------------------------
## The package theme. It follows the look of the tidyplots figures: a minimal
## frame with axis lines and ticks but no panel grid, grey type on white,
## left-aligned title/subtitle/caption, and a legend without a box. Every plot
## method in the package builds on it so the figures read as one family. It is
## first used by the uniformity-trial field map (see plot.trial_check()) and
## then adopted by the model plots (LRP, QRP, MCM, Paranaiba, replicates).
## ============================================================================

## Theme (internal) -----------------------------------------------------------
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
