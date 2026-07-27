# cran-comments

## Submission

This is a new submission.

## Test environments

* local Windows 11 x64, R 4.6.1, `R CMD check --as-cran`

## R CMD check results

Status: 1 NOTE.

```
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Willyan Bandeira <bandeira.wjab@gmail.com>'
  New submission
```

This is the expected note for a first submission. There were no ERRORs or
WARNINGs.

The local run was made with `--no-manual`: the LaTeX installation on the
checking machine has no `makeindex`, so the PDF manual step cannot complete
there. The Rd files were verified separately with
`R CMD Rd2pdf --no-index`, which builds the manual without error.

## Downstream dependencies

There are none, this being a new submission.
