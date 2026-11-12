dec <- function(id, yr, vals) paste0(sprintf("%-8s%5d", id, yr), paste(sprintf("%6d", vals), collapse = ""))
read_q <- function(...) { r <- NULL; invisible(capture.output(r <- dpl_read_dec(...))); r }
w <- c(0.38, 0.41, 0.35, 0.29, 0.44, 0.37, 0.40, 0.33, 0.36, 0.39, 0.42, 0.31)
mixed <- c(dec("A01", -446, round(w[1:6] * 100)),  dec("A01", -440, c(round(w[7:12] * 100), 999)),
           dec("A02", -446, round(w[1:6] * 1000)), dec("A02", -440, c(round(w[7:12] * 1000), -9999)),
           dec("A03", -446, round(w[1:6] * 100)),  dec("A03", -440, c(round(w[7:12] * 100), 999)))

test_that("unit = 'auto' scales each series by its own terminator", {
  r <- read_q(mixed)
  expect_equal(r$A01, w[c(1:12)], tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(r$A02, r$A01)                          # 0.001 mm series not 10x
  expect_equal(r$A03, r$A01)
  expect_equal(unname(attr(r, "units")), c("0.01mm", "0.001mm", "0.01mm"))
  expect_equal(range(as.integer(rownames(r))), c(-446L, -435L))  # terminators dropped
  out <- capture.output(dpl_read_dec(mixed))
  expect_true(any(grepl("0.001mm 1 series: A02", out)))
})

test_that("a forced unit warns about series whose terminator disagrees", {
  expect_warning(r <- read_q(mixed, unit = "0.01mm"), "A02")
  expect_equal(r$A02, 10 * r$A01)
  expect_silent(read_q(mixed, unit = "mm"))            # raw values: nothing to contradict
})

test_that("999 inside a 1/1000 mm series is a ring, not the end of the line", {
  x <- c(dec("B1", 1990, c(1200, 999, 850, 1010, 940, 999, 870, 905, 960, 880)),
         dec("B1", 2000, c(910, -9999)))
  r <- read_q(x)
  expect_equal(nrow(r), 11L)
  expect_equal(r["1991", "B1"], 0.999)
  expect_equal(r["1995", "B1"], 0.999)
  expect_equal(r["2000", "B1"], 0.910)
})

test_that("-9999 inside a series is still a missing ring (microcore exports)", {
  x <- c(dec("M1", 2001, c(120, -9999, 95, 999)))
  r <- read_q(x)
  expect_true(is.na(r["2002", "M1"]))
  expect_equal(r["2003", "M1"], 0.95)
  expect_equal(unname(attr(r, "units")), "0.01mm")
})

test_that("unterminated series follow the terminated ones of the file", {
  x <- c(dec("T1", 2000, c(1200, 1100, -9999)), dec("U1", 2000, c(1300, 1250)))
  r <- read_q(x)
  expect_equal(unname(attr(r, "units")), c("0.001mm", "0.001mm"))
  expect_equal(r["2000", "U1"], 1.3)
})

test_that("CL-MIR still reads at 1/1000 mm without naming the unit", {
  f <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
  a <- read_q(f, stop_val = -9999L, unit = "0.001mm")
  b <- read_q(f)
  expect_equal(as.data.frame(b), as.data.frame(a), ignore_attr = TRUE)
})
