# Plot an LRP fit (publication style)

Draws the observed points, the fitted broken line, dotted guides to the
breakpoint and the plateau-model annotations, in the layout used in
plot-size articles: the model block (equations and \\R^2\\) centred at
the top and `Xo`/`CVxo` next to the breakpoint. Axis limits come from
the data.

## Usage

``` r
# S3 method for class 'lrp_fit'
plot(
  x,
  title = "Linear Plateau",
  annotate_model = TRUE,
  xlab = NULL,
  ylab = "CV (%)",
  decimal_mark = c(".", ","),
  cond_word = "if",
  digits_coef = 3,
  digits_stat = 2,
  point_size = 2.3,
  line_size = 0.8,
  point_colour = "black",
  line_colour = "black",
  bp_colour = "red",
  base_size = 12,
  label_size = 4.6,
  title_size = NULL,
  family = "sans",
  theme = NULL,
  save = FALSE,
  file = NULL,
  format = c("tiff", "png", "jpeg", "pdf", "eps"),
  dpi = 300,
  width = 18,
  height = 12,
  units = "cm",
  compression = "lzw",
  ...
)
```

## Arguments

- x:

  an object of class `"lrp_fit"`.

- title:

  plot title.

- annotate_model:

  logical; draw the model equations and statistics.

- xlab, ylab:

  axis titles. `xlab = NULL` uses "Plot size (m^2)".

- decimal_mark:

  decimal separator for the annotations, "." or ",".

- cond_word:

  conditional word in the equations (e.g. "if" or "se").

- digits_coef, digits_stat:

  decimals for the coefficients (a, b) and for the statistics (Xo, CVxo,
  R2).

- point_size, line_size:

  sizes of the points and the fitted line.

- point_colour, line_colour, bp_colour:

  colours of the points, the fitted line and the breakpoint marker.

- base_size:

  base font size for the theme; axis titles and text scale from it.

- label_size:

  size of the annotation text (equations, Xo, CVxo).

- title_size:

  size of the plot title. `NULL` uses `base_size`.

- family:

  font family for the theme and annotations.

- theme:

  optional ggplot2 theme object used instead of the package default (the
  shared `trialSizing` theme); the `title_size` tweak is applied on top.

- save:

  logical; if `TRUE`, write the figure to disk (default FALSE).

- file:

  output file name; if `NULL`, derived from `title`.

- format:

  one of "tiff", "png", "jpeg", "pdf", "eps". TIFF/PDF/EPS are preferred
  for journals; TIFF is written with LZW compression.

- dpi:

  resolution for raster formats.

- width, height, units:

  figure size.

- compression:

  TIFF compression (default "lzw").

- ...:

  ignored.

## Value

A `ggplot` object (invisibly when saved).
