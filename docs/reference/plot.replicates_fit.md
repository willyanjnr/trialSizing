# Plot the number of replications

Draws replications against the number of treatments, one line per LSD
level.

## Usage

``` r
# S3 method for class 'replicates_fit'
plot(
  x,
  y_var = c("r_optimal", "r_continuous"),
  title = "Number of replications",
  xlab = "Number of treatments",
  ylab = NULL,
  colour_lab = "LSD (%)",
  line_size = 0.8,
  point_size = 1.8,
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

  a `"replicates_fit"` object.

- y_var:

  which value to plot: `"r_optimal"` (integer, default) or
  `"r_continuous"` (the article's tabulated value).

- title, xlab, ylab, colour_lab:

  labels.

- line_size, point_size:

  line and point sizes.

- base_size, title_size, family, theme:

  theme controls.

- save, file, format, dpi, width, height, units, compression:

  saving controls.

- ...:

  ignored.

## Value

A `ggplot` object (invisibly when saved).
