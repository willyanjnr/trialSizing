# Plot MCM fits for several trials (paginated grid)

Plot MCM fits for several trials (paginated grid)

## Usage

``` r
# S3 method for class 'mcm_multi'
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

  an `"mcm_multi"` object.

- save, file, format, dpi, width, height, units, compression:

  saving controls.

- ...:

  styling arguments forwarded to
  [`plot.mcm_fit()`](https://willyanjnr.github.io/trialSizing/reference/plot.mcm_fit.md).

## Value

A list of page objects, invisibly.
