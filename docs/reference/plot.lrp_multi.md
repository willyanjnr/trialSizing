# Plot LRP fits for several trials (paginated grid)

Arranges the per-trial plots in a grid of at most 6 panels (3 rows by 2
columns). With more than 6 trials the panels are paginated: each page is
a separate 3x2 figure, and with `save = TRUE` each page is written to
its own file. Requires the patchwork package.

## Usage

``` r
# S3 method for class 'lrp_multi'
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

  an object of class `"lrp_multi"`.

- save:

  logical; write the page(s) to disk (default FALSE).

- file:

  base file name; page number and extension are appended when there is
  more than one page.

- format, dpi, width, height, units, compression:

  passed to the saver.

- ...:

  styling arguments forwarded to
  [`plot.lrp_fit()`](https://willyanjnr.github.io/trialSizing/reference/plot.lrp_fit.md)
  (for example `decimal_mark`, `cond_word`, `label_size`).

## Value

A list of page objects, invisibly.
