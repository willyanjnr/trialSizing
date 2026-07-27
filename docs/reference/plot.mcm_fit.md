# Plot an MCM fit (publication style)

Draws the observed points, the fitted power curve \\a X^{-b}\\, dotted
guides to the maximum-curvature point and the model annotation. Unlike
the plateau models there is no flat segment: the curve keeps decreasing.

## Usage

``` r
# S3 method for class 'mcm_fit'
plot(
  x,
  title = "Modified Maximum Curvature",
  annotate_model = TRUE,
  xlab = NULL,
  ylab = "CV (%)",
  decimal_mark = c(".", ","),
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

  an `"mcm_fit"` object.

- title, annotate_model, xlab, ylab, decimal_mark, digits_coef,
  digits_stat:

  as in
  [`plot.lrp_fit()`](https://willyanjnr.github.io/trialSizing/reference/plot.lrp_fit.md).

- point_size, line_size, point_colour, line_colour, bp_colour:

  aesthetics.

- base_size, label_size, title_size, family, theme:

  theme and text controls.

- save, file, format, dpi, width, height, units, compression:

  saving controls.

- ...:

  ignored.

## Value

A `ggplot` object (invisibly when saved).
