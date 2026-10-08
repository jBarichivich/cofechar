# Short series: .cof_spline must not crash on 4-year input, dpl_read_dec must
# turn within-series -999 into NA, dpl_cof must set aside series shorter than
# min_length, and dpl_short must give an honest verdict on them.

mk_rwl <- function(years, ...) {
  cols <- list(...)
  m <- matrix(NA_real_, length(years), length(cols),
              dimnames = list(as.character(years), names(cols)))
  for (nm in names(cols)) {
    v <- cols[[nm]]; m[seq_along(v) + attr(v, "offset"), nm] <- v
  }
  structure(as.data.frame(m), class = c("rwl", "data.frame"))
}
ser <- function(x, offset = 0L) { attr(x, "offset") <- offset; x }

test_that(".cof_spline runs for series of 3 to 8 years", {
  sp <- cofechar:::.cof_spline
  x  <- c(1.70, 0.63, 1.89, 1.20, 0.90, 1.40, 1.10, 1.30)
  for (n in 3:8) expect_length(sp(x[1:n], 32, 0.5), n)
})

test_that("dpl_read_dec maps within-series -999 to NA and keeps the years", {
  lines <- c("S1      2000    10    20  -999    40   999",
             "S2      2000    11    22    33   999")
  rwl <- dpl_read_dec(lines, unit = "mm")
  expect_true(is.na(rwl["2002", "S1"]))
  expect_equal(rwl["2003", "S1"], 40)
  expect_equal(sum(!is.na(rwl$S2)), 3L)
  # na_val = NULL keeps the raw value
  raw <- dpl_read_dec(lines, unit = "mm", na_val = NULL)
  expect_equal(raw["2002", "S1"], -999)
})

test_that("dpl_cof sets aside series shorter than min_length", {
  set.seed(1)
  sig <- cumsum(rnorm(60)); sig <- abs(sig / max(abs(sig)) * 2) + 0.5
  rwl <- mk_rwl(1950:2009,
                A = ser(sig + rnorm(60, sd = 0.2)),
                B = ser(sig + rnorm(60, sd = 0.2)),
                C = ser(sig + rnorm(60, sd = 0.2)),
                S = ser(sig[55:60], offset = 54L))        # 6-year series
  cof <- dpl_cof(rwl, verbose = FALSE, parts = 1L)
  expect_equal(cof$short$series, "S")
  expect_equal(cof$short$n, 6L)
  expect_false("S" %in% cof$stats$series)
  expect_equal(nrow(cof$stats), 3L)
  expect_true(any(grepl("SERIES NOT USED", cof$output)))
  expect_equal(cof$options$min_length, 10L)
  # min_length = 1 brings it back
  cof2 <- dpl_cof(rwl, verbose = FALSE, parts = integer(0), min_length = 1L)
  expect_equal(nrow(cof2$short), 0L)
  expect_true("S" %in% cof2$stats$series)
})

test_that("dpl_short finds lag 0 for a correctly dated short series and the shift otherwise", {
  set.seed(2)
  sig <- cumsum(rnorm(80)); sig <- abs(sig / max(abs(sig)) * 2) + 0.5
  good <- sig[71:80] * exp(rnorm(10, sd = 0.05))      # 10 yr, dated right
  shft <- sig[70:79] * exp(rnorm(10, sd = 0.05))      # same rings labelled 1 yr too old
  rwl <- mk_rwl(1931:2010,
                A = ser(sig * exp(rnorm(80, sd = 0.15))),
                B = ser(sig * exp(rnorm(80, sd = 0.15))),
                C = ser(sig * exp(rnorm(80, sd = 0.15))),
                D = ser(sig * exp(rnorm(80, sd = 0.15))),
                G = ser(good, offset = 70L),
                H = ser(shft, offset = 70L))
  cof <- dpl_cof(rwl, verbose = FALSE, parts = integer(0), min_length = 11L)
  expect_setequal(cof$short$series, c("G", "H"))
  chk <- dpl_short(rwl, cof, max_lag = 2L)
  s <- chk$summary
  expect_equal(s$verdict[s$series == "G"], "ok")
  expect_equal(s$best_lag[s$series == "G"], 0L)
  expect_equal(s$verdict[s$series == "H"], "shifted")
  expect_equal(s$best_lag[s$series == "H"], -1L)   # true year = label + lag
  # p-values only from min_n differences
  expect_true(all(is.na(chk$lags$p[chk$lags$n < 5])))
})

test_that("dpl_short pool = TRUE tests the stand mean and labels 3-year series untestable", {
  set.seed(3)
  sig <- cumsum(rnorm(80)); sig <- abs(sig / max(abs(sig)) * 2) + 0.5
  mk <- function(from, to) ser(sig[from:to] * exp(rnorm(to - from + 1, sd = 0.10)),
                                offset = from - 1L)
  rwl <- mk_rwl(1931:2010,
                A = ser(sig * exp(rnorm(80, sd = 0.15))),
                B = ser(sig * exp(rnorm(80, sd = 0.15))),
                C = ser(sig * exp(rnorm(80, sd = 0.15))),
                D = ser(sig * exp(rnorm(80, sd = 0.15))),
                s1 = mk(76, 80), s2 = mk(75, 80), s3 = mk(74, 80),   # 5-7 yr
                s4 = mk(77, 80), s5 = mk(78, 80))                     # 4, 3 yr
  cof <- dpl_cof(rwl, verbose = FALSE, parts = integer(0))
  expect_setequal(cof$short$series, paste0("s", 1:5))
  chk <- dpl_short(rwl, cof, pool = TRUE, pool_min = 2L)
  s <- chk$summary
  expect_equal(s$series[1L], "POOL")
  expect_equal(s$best_lag[s$series == "POOL"], 0L)
  expect_equal(s$verdict[s$series == "POOL"], "ok")
  # the pool is testable although its members are not
  expect_false(is.na(s$p0[s$series == "POOL"]))
  expect_equal(s$verdict[s$series == "s5"], "untestable")
  expect_true(is.na(s$p0[s$series == "s4"]))          # 4 yr: no p-value
  expect_true(!is.null(chk$pool))
  expect_true(all(chk$pool$depth >= 2L))
  # pooled years cover 1976..1980 differences (>= 2 series)
  expect_equal(range(as.integer(names(chk$pool$series))), c(2006L, 2010L))
})

test_that("dpl_cof runs when the last series in the collection are set aside", {
  set.seed(4)
  sig <- cumsum(rnorm(60)); sig <- abs(sig / max(abs(sig)) * 2) + 0.5
  rwl <- mk_rwl(1950:2009,
                A = ser(sig * exp(rnorm(60, sd = 0.15))),
                B = ser(sig * exp(rnorm(60, sd = 0.15))),
                C = ser(sig * exp(rnorm(60, sd = 0.15))),
                S1 = ser(sig[53:60], offset = 52L),     # 8 yr, last columns
                S2 = ser(sig[55:60], offset = 54L))     # 6 yr
  # min_length larger than the short series, Part 5 and 6 requested
  cof <- dpl_cof(rwl, verbose = FALSE, min_length = 15L)
  expect_setequal(cof$short$series, c("S1", "S2"))
  expect_equal(nrow(cof$stats), 3L)
  expect_true(any(grepl("^PART 5", cof$output)))
  expect_true(any(grepl("^PART 7", cof$output)))
})

test_that("dpl_short_barplot prints master and samples on the master's scale", {
  set.seed(5)
  sig <- cumsum(rnorm(60)); sig <- abs(sig / max(abs(sig)) * 2) + 0.5
  rwl <- mk_rwl(1950:2009,
                A = ser(sig * exp(rnorm(60, sd = 0.15))),
                B = ser(sig * exp(rnorm(60, sd = 0.15))),
                C = ser(sig * exp(rnorm(60, sd = 0.15))),
                S1 = ser(sig[54:60], offset = 53L),
                S2 = ser(sig[56:60], offset = 55L))
  cof <- dpl_cof(rwl, verbose = FALSE, parts = integer(0))
  out <- dpl_short_barplot(rwl, cof, quiet = TRUE)
  expect_true(any(grepl("MASTER", out)))
  expect_true(any(grepl("POOL", out)))
  expect_true(any(grepl("^2009  ", out)))
  # first year of a sample is blank (differences): the 2003 row has no S1 bar
  r2003 <- out[grepl("^2003  ", out)]
  expect_true(nchar(r2003) < 6 + 2 * 12 + 1)
  # lag shifts the column
  out2 <- dpl_short_barplot(rwl, cof, series = "S2", pool = FALSE,
                         lag = c(S2 = -1), quiet = TRUE)
  expect_true(any(grepl("2004-2008", out2)))
  # explicit window
  out3 <- dpl_short_barplot(rwl, cof, years = c(2000, 2009), quiet = TRUE)
  expect_false(any(grepl("^1999  ", out3)))
  # file output
  f <- tempfile(); dpl_short_barplot(rwl, cof, quiet = TRUE, output_file = f)
  expect_equal(readLines(f), out)
})
