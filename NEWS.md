# cofechar 1.0.1 (development)

## Bug fixes

* `dpl_barplot(layout = "page")` advanced pages by 100 years while each page
  spans 400 years (8 columns x 50 rows), so every year was printed up to four
  times and the columns were not anchored on COFECHA's 400-year grid. Pages
  now follow the DPL BARPL layout exactly.
* COFECHA Part 4 and `dpl_barplot` now reproduce the DPL bar plot line for
  line: decade separator rows after every 10 years, no blank line after the
  column header, page header repeated per 400-year page, no trailing footer
  in Part 4.
* Decile cut-points used Rs `round()` (half-to-even) where Holmes' Fortran
  uses `NINT` (half-away-from-zero); bars at exact-.5 index boundaries were
  one dash short. New internal `.cof_nint()` fixes this.
* `dpl_barplot()` now returns the documented list of formatted line vectors
  (one element per series for `layout = "page"`, `all` for `"column"`).

## Internal

* Part 4 and `dpl_barplot` share one BARPL engine (`.cof_barpl_cuts()`,
  `.cof_barpl_car()`, `.cof_barpl_pages()`) so the two cannot drift apart.
* New regression test against the DPL COFECHA 4.04P reference run on the
  PEL dataset (`inst/extdata/benchmark/`): Part 7 statistics, Part 4 line
  for line, and `dpl_barplot` bar cells.

# cofechar 1.0.0

* Initial CRAN-ready package release.
* Roxygen2 documentation for all 12 exported functions.
* Full test suite (testthat) covering EDT operations and COF internals.
* Vignette: "Getting started with cofechar".
* All internal statistical routines verified against PUE (Puerco) DPL
  benchmark: Part 4 bar plot matches 67/67 bars; segment correlations
  match reference output.