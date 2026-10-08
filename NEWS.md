
# cofechar 1.0.0.9000 (development)

## New features

* `dpl_short()`: anchored dating check for series too short for COFECHA's
  sliding segments (microcores, xylogenesis samples, short anatomical
  sequences). Tests lag 0 against +/- `max_lag` on first differences of the
  log series, reporting r, an exact p-value (n >= 5), Gleichlaeufigkeit with
  its binomial p-value, and a verdict (ok / shifted / weak).
* `dpl_cof()` gains `min_length` (default 10): series with fewer measured
  years are excluded from the master and from segment testing, listed in
  Part 1 and returned in `$short`. COFECHA's CRIT99 table starts at 10 years;
  shorter series produced meaningless unflagged correlations (r > 0.98 on
  3 or 4 values). Legacy behaviour: `min_length = 1`.
* `dpl_read_dec()` gains `na_val` (default `c(-999, -9999)`): within-series
  missing markers become `NA` instead of -9.99 mm.

## Bug fixes

* `.cof_spline()` crashed with "subscript out of bounds" on series of
  exactly 4 years (R's `a:b` counts backwards where a Fortran DO loop would
  not execute).

# cofechar 1.0.0

First release: a faithful R port of Holmes' COFECHA and EDT modules from the
Dendrochronology Program Library (DPL), validated line for line against the
DPL COFECHA 4.04P reference run on the PEL dataset (`inst/extdata/benchmark/`)
and against the PUE benchmark.

* 12 exported functions (`dpl_read`, `dpl_read_dec`, `dpl_write`, `dpl_merge`,
  `dpl_edt`, `dpl_display`, `dpl_edit_file`, `dpl_cof`, `dpl_cof_diag`,
  `dpl_dateme`, `dpl_barplot`, `plot_cof_barplot`), all taking and returning
  dplR `rwl` data.frames.
* Internal routines (cubic smoothing spline, Burg AR modelling, log-transform,
  normalisation, BARPL decile bars using Fortran `NINT` rounding) translated
  directly from the Fortran source.
* Part 4 and `dpl_barplot` share one BARPL page engine (400-year pages,
  8 columns x 50 rows). One intentional departure from DPL: the decade
  separator is printed as a blank line, as in Part 3, instead of a row of
  `----`.
* Example datasets: `CL-MIR.rwl` (36 *Fitzroya cupressoides* series,
  Cerro Mirador, 1406--2002; Barichivich 2005) and the PEL benchmark.
* Test suite (testthat): EDT operations, COF internals, and the PEL
  benchmark regression (Part 7 statistics, Part 4 line for line).
