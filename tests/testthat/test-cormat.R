pel_rwl <- function() {
  f <- system.file("extdata", "benchmark", "PEL.rwl", package = "cofechar")
  dpl_read_dec(f, stop_val = 999L, unit = "0.01mm")
}

test_that("dpl_cof returns the transformed series as $filtered on the master axis", {
  rwl <- pel_rwl()
  cof <- dpl_cof(rwl, seg_length = 20L, seg_lag = 10L, verbose = FALSE, parts = integer(0))
  f <- cof$filtered
  expect_s3_class(f, "rwl")
  expect_equal(rownames(f), names(cof$master))
  expect_equal(colnames(f), cof$stats$series)
  # normalised: mean ~0, population sd ~1 per series
  z <- f$PEL11A[!is.na(f$PEL11A)]
  expect_equal(mean(z), 0, tolerance = 1e-5)
  expect_equal(sqrt(mean((z - mean(z))^2)), 1, tolerance = 1e-5)
})

test_that("dpl_cormat raw and transformed agree with cor() and are symmetric", {
  rwl <- pel_rwl()
  cm  <- dpl_cormat(rwl, quiet = TRUE)
  expect_equal(cm$type, "raw")
  expect_true(isSymmetric(cm$r))
  expect_equal(unname(diag(cm$r)), rep(1, ncol(rwl)))
  expect_equal(cm$r["PEL06B", "PEL01A"],
               cor(rwl$PEL06B, rwl$PEL01A, use = "pairwise.complete.obs"), tolerance = 1e-6)
  cof <- dpl_cof(rwl, seg_length = 20L, seg_lag = 10L, verbose = FALSE, parts = integer(0))
  ct  <- dpl_cormat(cof, quiet = TRUE)
  expect_equal(ct$type, "transformed")
  expect_equal(ct$r["PEL06B", "PEL01A"],
               cor(cof$filtered$PEL06B, cof$filtered$PEL01A, use = "pairwise.complete.obs"),
               tolerance = 1e-6)
  # all PEL pairs are positive; the strongest are significant, and the two
  # weakest involve PEL14A / PEL07A, the series DPL also rates lowest
  expect_true(all(ct$r[upper.tri(ct$r)] > 0))
  expect_lt(ct$p["PEL06B", "PEL02A"], 0.001)
  expect_gt(ct$p["PEL01A", "PEL14A"], 0.05)
  expect_equal(ct$n["PEL11A", "PEL11A"], 161L)
  expect_equal(ct$n["PEL14A", "PEL07A"], 55L)
  expect_equal(nrow(ct$summary), 7L)
  expect_true(ct$rbar > 0.3 && ct$rbar < 0.9)
})

test_that("dpl_cormat prints Part-5 style cells and flags a duplicated sample", {
  rwl <- pel_rwl()
  dup <- rwl; dup$PEL06B_copy <- rwl$PEL06B * 1.02      # same core, re-measured
  cm  <- dpl_cormat(dup, quiet = TRUE)
  expect_gt(cm$r["PEL06B", "PEL06B_copy"], 0.999)
  expect_true(any(grepl("1.00\\*\\*\\*", cm$lines)))
  # leading zero suppressed in the matrix cells, Fortran style
  cells <- cm$lines[grepl("^PEL", cm$lines)]
  expect_true(any(grepl(" \\.[0-9]{2}\\*", cells)))
  expect_false(any(grepl(" 0\\.[0-9]{2}", cells)))
  # min_overlap: no PEL pair has < 50 common years except with the copy
  cm2 <- dpl_cormat(rwl, min_overlap = 60L, quiet = TRUE)
  expect_true(any(is.na(cm2$r[upper.tri(cm2$r)])))
  expect_true(any(grepl("  -$| - ", cm2$lines)))
  # series by position, digits = 3, file output
  f <- tempfile()
  cm3 <- dpl_cormat(rwl, series = 1:3, digits = 3, quiet = TRUE, output_file = f)
  expect_equal(dim(cm3$r), c(3L, 3L))
  expect_true(any(grepl("\\.[0-9]{3}", cm3$lines)))
  expect_equal(readLines(f), cm3$lines)
})
