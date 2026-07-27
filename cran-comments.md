# cran-comments

## Resubmission

This is a resubmission. The pre-test of 2026-07-27 reported two items; both are
addressed below.

### Overall checktime 26 min > 10 min

Fixed. The cost was concentrated in `fit_lrp()` and `fit_qrp()`, whose
breakpoint is found by a grid search whose resolution is set by `step`
(default `0.001`). Tests, vignettes and examples were calling them without
`step`, so every illustrative fit paid for full precision.

`step` is now explicit everywhere. The default `0.001` is kept only where the
third decimal is actually asserted, namely the two "stable to three decimals"
regression tests and the test that checks `search_range` leaves the optimum
untouched. Everywhere else the fits use `step = 0.01`, which resolves the
breakpoint to two decimals and is what the package's own documentation
recommends for routine use. The bootstrap examples were reduced from 500 to
200 replicates, matching `compare_methods()`.

No user-facing default changed: `fit_lrp()` and `fit_qrp()` still default to
`step = 0.001`.

Measured locally, on the timings the check itself reports:

| Stage                     | Before | After |
|---------------------------|-------:|------:|
| tests                     |  390 s |  88 s |
| re-building vignettes     |  341 s |  94 s |
| examples (--run-donttest) |  110 s |  57 s |
| whole check, wall clock   |  ~15 m | 5 m 21 s |

### Possibly misspelled words in DESCRIPTION

These are spelled as intended and are not errors:

* `Cargnelutti`, `Ferreira`, `Filho`, `Lessman`, `Morais`, `Paranaiba` --
  surnames of the authors of the cited methods. `Paranaiba` is the ASCII
  transliteration of "Paranaíba", used so the DESCRIPTION stays ASCII.
* `LRP`, `QRP` -- acronyms for "linear response plateau" and "quadratic
  response plateau", both expanded in the DESCRIPTION where they first appear.
* `semivariogram`, `kriged` -- standard geostatistical terms.

`Language: en` was added to the DESCRIPTION.

## Test environments

* local Windows 11 x64, R 4.6.1, `R CMD check --as-cran`

## R CMD check results

No ERRORs or WARNINGs. One NOTE, the expected `New submission`.

The local run uses `--no-manual`: the LaTeX installation on this machine has no
`makeindex`, so the PDF manual step cannot complete here. The Rd files were
verified with `R CMD Rd2pdf --no-index`, which builds the manual without error.

## Downstream dependencies

There are none, this being a new submission.
