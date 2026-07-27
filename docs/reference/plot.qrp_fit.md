# Plot a QRP fit (publication style)

Same layout and options as
[`plot.lrp_fit()`](https://willyanjnr.github.io/trialSizing/reference/plot.lrp_fit.md),
with the quadratic descending arm and the quadratic equation in the
annotation.

## Usage

``` r
# S3 method for class 'qrp_fit'
plot(
  x,
  title = "Quadratic Plateau",
  annotate_model = TRUE,
  xlab = NULL,
  ylab = "CV (%)",
  decimal_mark = c(".", ","),
  cond_word = "if",
  digits_coef = 3,
  digits_c = 4,
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

  a `"qrp_fit"` object.

- title:

  plot title.

- annotate_model:

  draw the model equations and statistics.

- xlab, ylab:

  axis titles. `xlab = NULL` uses "Plot size (m^2)".

- decimal_mark:

  decimal separator, "." or ",".

- cond_word:

  conditional word in the equations (e.g. "if" or "se").

- digits_coef, digits_c, digits_stat:

  decimals for a/b, for c, and for the statistics.

- point_size, line_size:

  sizes of points and fitted line.

- point_colour, line_colour, bp_colour:

  point, line and breakpoint colours.

- base_size, label_size, title_size, family:

  theme and annotation sizes and font family.

- theme:

  optional ggplot2 theme used instead of the default.

- save, file, format, dpi, width, height, units, compression:

  saving controls; see
  [`plot.lrp_fit()`](https://willyanjnr.github.io/trialSizing/reference/plot.lrp_fit.md).

- ...:

  ignored.

## Value

A `ggplot` object (invisibly when saved).
