# Regression test against the DPL COFECHA 4.04P reference run on the PEL
# dataset (7 Pilgerodendron series, 1862-2022; inst/extdata/benchmark).

bench_lines <- function() {
  f <- system.file("extdata", "benchmark", "PELCOF.OUT", package = "cofechar")
  skip_if(f == "", "benchmark output not installed")
  # DPL output contains form-feed page breaks: split them into lines
  unlist(strsplit(readLines(f, warn = FALSE), "\f", fixed = TRUE))
}

pel_cof <- function() {
  f   <- system.file("extdata", "benchmark", "PEL.rwl", package = "cofechar")
  rwl <- dpl_read_dec(f, stop_val = 999L, unit = "0.01mm")
  dpl_cof(rwl, seg_length = 20L, seg_lag = 10L, spline_period = 32L,
          verbose = FALSE)
}

extract_part <- function(lines, part) {
  a <- grep(sprintf("^\\s*PART %d", part), lines)[1L]
  b <- grep(sprintf("^\\s*PART %d", part + 1L), lines)[1L]
  sub("\\s+$", "", lines[a:(b - 1L)])
}

bar_cells <- function(lines) {
  m <- regmatches(lines, gregexpr(" (\\d{4})(-*)([A-Za-z@<>])", lines))
  m <- unlist(m)
  yr  <- as.integer(substr(m, 2L, 5L))
  bar <- substr(m, 6L, nchar(m))
  setNames(bar, yr)
}

test_that("PEL benchmark: Part 7 series statistics reproduce DPL output", {
  cof <- pel_cof()
  expect_equal(cof$crit, 0.5155)
  expect_equal(round(cof$stats$r_master, 3),
               c(.639, .563, .696, .541, .484, .556, .462))
  expect_equal(cof$stats$n_flags, c(1L, 2L, 1L, 3L, 3L, 2L, 2L))
  expect_equal(cof$stats$n_segs,  c(10L, 8L, 8L, 10L, 6L, 7L, 6L))
})

test_that("PEL benchmark: Part 4 bar plot matches DPL line for line", {
  cof   <- pel_cof()
  bench <- extract_part(bench_lines(), 4L)
  ours  <- extract_part(cof$output, 4L)

  # Page header line carries a timestamp in DPL and the run title here;
  # the horizontal rule differs by Fortran's carriage-control column.
  mask <- function(x) {
    x <- sub("^\\s*PART 4:.*$", "PART 4", x)
    sub("^ ?-{131,132}$", "RULE", x)
  }
  expect_equal(length(ours), length(bench))
  expect_equal(mask(ours), mask(bench))

  b <- bar_cells(bench); o <- bar_cells(ours)
  expect_equal(length(o), 161L)
  expect_equal(o[names(b)], b)
})

test_that("dpl_barplot(cof) reproduces Part 4 bars and prints each year once", {
  cof   <- pel_cof()
  bench <- bar_cells(extract_part(bench_lines(), 4L))
  out   <- dpl_barplot(cof, quiet = TRUE)[["master"]]
  o     <- bar_cells(out)
  expect_equal(length(o), 161L)
  expect_false(any(duplicated(names(o))))
  expect_equal(o[names(bench)], bench)
  # decade separators present, header not followed by a blank line
  expect_true(any(grepl("^   ----            ----", out)))
})

test_that(".cof_nint rounds half away from zero like Fortran NINT", {
  expect_equal(cofechar:::.cof_nint(c(80.5, -80.5, 2.4, -2.6, 0.5, -0.5)),
               c(81, -81, 2, -3, 1, -1))
  expect_equal(round(80.5), 80)   # documents why round() cannot be used
})
