# Plot QRP fits for several trials (paginated grid)

Arranges the per-trial plots in a grid of at most 6 panels (3 x 2),
paginating when there are more. Requires patchwork. See
[`plot.lrp_multi()`](https://willyanjnr.github.io/trialSizing/reference/plot.lrp_multi.md).

## Usage

``` r
# S3 method for class 'qrp_multi'
plot(
  x,
  save = FALSE,
  file = NULL,
  format = c("tiff", "png", "jpeg", "pdf", "eps"),
  dpi = 300,
  width = 18,
  height = 24,
  units = "cm",
  compression = "lzw",
  ...
)
```

## Arguments

- x:

  a `"qrp_multi"` object.

- save, file, format, dpi, width, height, units, compression:

  saving controls.

- ...:

  styling arguments forwarded to
  [`plot.qrp_fit()`](https://willyanjnr.github.io/trialSizing/reference/plot.qrp_fit.md).

## Value

A list of page objects, invisibly.
