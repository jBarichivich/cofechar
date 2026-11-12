# cofechar 1.0.0.9018 (development)

## Mixed-unit files

* `dpl_read_dec(unit = "auto")` now scales each series by its own
  terminator instead of choosing one unit for the whole file: a series ending
  in 999 is read in 1/100 mm and one ending in -9999 in 1/1000 mm, as in the
  Tucson convention. Before, a 1/1000 mm series in a file of 1/100 mm series
  came back 10 times too large, and its -9999 terminator was read as a
  missing ring, adding an empty final year. The package's own `CL-MIR.rwl`
  read with default arguments had both problems; it is now correct without
  naming the unit.
* `stop_val` defaults to `c(999L, -9999L)`. A stop code is the terminator
  only when it is the last value of its series. Elsewhere it is data, so
  999 inside a 1/1000 mm series is a 0.999 mm ring and no longer silently
  drops the rest of its line. -9999 inside a series is still a missing ring.
* The unit of each series is returned in `attr(x, "units")` and printed when
  a file mixes units. Forcing `unit = "0.01mm"` or `"0.001mm"` warns about
  series whose terminator implies the other unit.

## Multi-panel figures

* `plot_cof_barplot()` and `plot_cof_coverage()` take one panel of the
  caller's `par(mfrow)` or `layout()` figure. They used to save and restore
  the whole `par()` state, which rewound the device to the panel they started
  in, so the next plot overwrote them; `plot_cof_coverage()` also called
  `layout()`, replacing the caller's grid. The coverage plot now stacks its
  series and depth panels inside its own panel.
* `plot_cof_barplot()` draws the sample-depth line in the bar plot's own
  coordinates. It was drawn on a separate year axis that did not line up with
  the bars.
* `plot_cof_barplot()`'s default title uses a hyphen between the years: the
  en dash it used printed as dots on `pdf()` devices.


## Years before the common era

* `dpl_read_dec()` reads negative (BCE) years in every layout: spaced
  (`AMI02   -446 ...`), compact with the year glued to the label
  (`AMI02XYZ-446 ...`), and five-character years that eat into the label
  field (`AMI02-1000 ...`). Series may run across year 0. Years are read as
  written (astronomical numbering, as in the Tucson format). A line that could
  be either a dashed ID or a negative year (`BAR-01  1930`) is resolved per
  file by which reading keeps consecutive lines in the same series.
* `label_length` now does what the error message promised: the first token is
  split after `label_length` characters and the rest is taken as the year,
  signed or not.
* **Behaviour change:** `base_century` now defaults to `NULL` (years read
  literally). The old default of 1900 silently moved every year below 1000:
  446 CE became 2346 and -446 became 1454. Pass `base_century = 1900L` for
  files that store only the year within the century.
* `dpl_cof()`, `dpl_barplot()`, `plot_cof_rings()` and `plot_cof_coverage()`
  are tested end to end on a BCE collection.

## Accented text in graphics

* The plotting functions (`plot_cof_rings()`, `plot_cof_coverage()`,
  `plot_cof_barplot()`, `dpl_print_pdf()`) switch `LC_CTYPE` to UTF-8 for
  the duration of the call when the session runs in a C/POSIX locale
  (Rscript on a server or container, RStudio Server with `LANG` unset), where
  every accented character was drawn as dots ("Regi..n"). The session locale
  is restored on exit; Latin-1 and Windows sessions are left untouched.
* `dpl_print_pdf()` writes with the WinAnsi (Windows-1252) encoding, so en
  and em dashes, curly quotes and the euro sign print as well as Latin-1
  accents.
* `dpl_print_pdf()` no longer opens a stray `Rplots.pdf` device on exit, and
  `plot_cof_barplot()` no longer discards other exit handlers when restoring
  `par()`.

## New features

* `dpl_short()`: anchored dating check for series too short for COFECHA's
  sliding segments (microcores, xylogenesis samples, short anatomical
  sequences). Tests lag 0 against +/- `max_lag` on first differences of the
  log series, reporting r, an exact p-value (n >= 5), Gleichlaeufigkeit with
  its binomial p-value, and a verdict (ok / shifted / weak / untestable).
  `pool = TRUE` also tests the mean log-difference series of all short
  series from the stand, which can be tested even when the individual
  series are too short. The help page documents the workflow and the
  number of years needed for each level of evidence, and what
  Gleichlaeufigkeit measures, its exact binomial significance and how to
  weigh it against r on short series. `glk_p0` is reported in the summary.
* `dpl_short_barplot()`: bar plot of the master and short series side by side,
  one row per year, in COFECHA notation. Samples are not ranked over their
  own few years but expressed on the master's scale (log first differences
  standardised by the master's window mean and SD), so bars are comparable
  across columns and pointer years line up. `pool = TRUE` adds the stand
  mean; `lag` previews a shift suggested by `dpl_short()`.
* `dpl_trim()`: drop the leading and trailing all-`NA` years of an `rwl`,
  so a column subset (`rwl[, c("A", "B")]`, which keeps the full year axis)
  spans only the selected series; `dpl_trim(rwl, series = 1:5)` subsets
  and trims in one step. `dpl_merge()` now trims by default (`trim = FALSE`
  for the old behaviour); `dpl_edt()` already did.
* `dpl_edt()` and `dpl_edit_file()` gain `keep` and `drop`: select the
  samples to copy (or to leave out) by ID or position in one argument, e.g.
  `dpl_edit_file("site.rwl", output_path = "sub.rwl", keep = c("A01", "A02"))`
  or `keep = 1:5`. Replaces the DPLEDT idiom of one `copy` edit per series
  with `default_action = "omit"`, which still works.
* `dpl_cormat()`: pairwise correlation matrix of all series printed in the
  COFECHA Part 5 style (lower triangle, F4.2 cells with the leading zero
  suppressed, `-` below `min_overlap`), with significance asterisks
  (* 0.05, ** 0.01, *** 0.001; one-sided, n - 2 df on each pair's own
  overlap) and per-series summary rows (pairs, mean r, number significant)
  plus the collection rbar. `type = "raw"` uses the measurements,
  `type = "transformed"` the COFECHA-filtered series.
* `dpl_cof()` returns `$filtered`: the fully transformed series (spline,
  log, AR, normalised) as an `rwl` on the master's year axis - the values
  Part 5 correlates - for `dpl_cormat()` and for the user's own analyses.
* `dpl_print_pdf()`: print the COFECHA text output (`cof$output`) or any
  ASCII bar plot (`dpl_barplot`, `dpl_short_barplot`) to an A4/Letter PDF
  in a monospaced font. The font size is computed so that the longest line
  fits the page width and a complete 400-year bar-plot page fits one sheet;
  page breaks fall on PART headers, series blocks and bar-plot pages, never
  inside a block. Orientation is chosen automatically (landscape for
  132-column output, portrait for narrower text). `parts = 4` prints only
  the master bar plot of a COFECHA output; any subset of 1:7 can be chosen.
* `dpl_barplot()`, `dpl_short()` and `dpl_short_barplot()` accept integer
  positions in `series` (e.g. `1:7`) as well as IDs. For `dpl_barplot()`
  positions index the columns of `rwl`; for the two short-series functions
  they index the default set (`cof$short`), so `1:7` means the first seven
  short series. Unknown IDs warn and are dropped; positions out of range
  are an error.
* `dpl_cof()` gains `min_length` (default 10): series with fewer measured
  years are excluded from the master and from segment testing, listed in
  Part 1 and returned in `$short`. COFECHA's CRIT99 table starts at 10 years;
  shorter series produced meaningless unflagged correlations (r > 0.98 on
  3 or 4 values). Legacy behaviour: `min_length = 1`.
* `dpl_read_dec()` gains `na_val` (default `c(-999, -9999)`): within-series
  missing markers become `NA` instead of -9.99 mm.

* `plot_cof_rings()`: ring-width diagram of a collection in the TSAP-Win
  layout -- one bar per core spanning its years, one cell per ring shaded by
  its within-series decile (dark = narrow) or by the COFECHA-transformed
  index, pointer years marked where most cores agree, optional vertical rules
  at known event years, sorting by first/last year or length, and grouping of
  cores by tree. `xaxis = "length"` draws every cell as wide as its ring, so
  each bar is a scale drawing of the core, aligned at the bark or the pith
  (`align`), with pointer years as coloured strips inside the cells and year
  ticks above each bar. Returns the decile matrix and the pointer-year table.

* `plot_cof_coverage()`: timelines of the series (one line per core,
  stacked, with the sample depth underneath) on the calendar axis, or aligned
  at the first ring (cambial age, the RCS view, with optional pith offsets
  and a replication threshold line) or at the last ring.

## Bug fixes

* `dpl_cof()` failed with "subscript out of bounds" in Part 5 when the
  last series of the collection were set aside (empty or shorter than
  `min_length`); the Part 5 row list was not pre-sized.
* `.cof_spline()` crashed with "subscript out of bounds" on series of
  exactly 4 years (R's `a:b` counts backwards where a Fortran DO loop would
  not execute).

## Documentation

* `dpl_edt()` examples now cover every DPLEDT operation on the CL-MIR data:
  replace, first_year / last_year redating, insert and delete with both
  `move` directions, trim_start / trim_end, rename followed by a positional
  edit, a multi-edit call, and `keep` combined with `edits`.

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
