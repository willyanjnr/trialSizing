# Plot Paranaiba estimates by trial

Shows the optimal plot size per trial, with the mean across trials as a
dashed reference line.

## Usage

``` r
# S3 method for class 'paranaiba_fit'
plot(
  x,
  y_var = c("Xo", "CVxo"),
  title = "Paranaiba method",
  xlab = "Trial",
  ylab = NULL,
  show_mean = TRUE,
  fill_colour = "grey35",
  base_size = 12,
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

  a `"paranaiba_fit"` object.

- y_var:

  `"Xo"` (default) or `"CVxo"`.

- title, xlab, ylab:

  labels; `ylab = NULL` is chosen from `y_var`.

- show_mean:

  draw the across-trial mean as a dashed line.

- fill_colour:

  bar fill colour.

- base_size, title_size, family, theme:

  theme controls.

- save, file, format, dpi, width, height, units, compression:

  saving controls.

- ...:

  ignored.

## Value

A `ggplot` object (invisibly when saved).
