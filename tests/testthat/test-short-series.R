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
