test_that(".cof_normts returns mean=0 sd=1 (population)", {
  x   <- c(1.1, 0.9, 1.3, 0.8, 1.2)
  res <- cofechar:::.cof_normts(x, k = 0L)
  expect_equal(mean(res$z), 0, tolerance = 1e-10)
  expect_equal(sd(res$z) * sqrt((length(x) - 1) / length(x)),
               1, tolerance = 1e-10)
})

test_that(".cof_spline returns a vector of the same length", {
  x  <- c(1.2, 0.9, 1.4, 1.1, 0.8, 1.3, 1.0, 1.2, 0.9, 1.1)
  sp <- cofechar:::.cof_spline(x, zsp = 5, pvp = 0.5)
  expect_equal(length(sp), length(x))
  expect_true(is.numeric(sp))
})

test_that(".cof_correl matches stats::cor on a clean pair", {
  set.seed(42)
  x <- rnorm(30); y <- x * 0.7 + rnorm(30, sd = 0.3)
  r_cof   <- cofechar:::.cof_correl(x, y)
  r_stats <- cor(x, y)
  expect_equal(r_cof, r_stats, tolerance = 1e-6)
})

test_that(".cof_correl returns -9 for length-1 input", {
  expect_equal(cofechar:::.cof_correl(1, 1), -9.0)
})

test_that(".cof_crit99 returns decreasing values for increasing segment length", {
  crits <- sapply(c(15, 25, 50, 100), cofechar:::.cof_crit99)
  expect_true(all(diff(crits) < 0))
})

test_that(".cof_mempr returns residuals of same length as input", {
  x   <- cumsum(rnorm(50)) # random walk for AR structure
  res <- cofechar:::.cof_mempr(x, lg = 5L)
  expect_equal(length(res$residuals), length(x))
  expect_true(res$order >= 0L)
})

test_that("dpl_cof runs without error on minimal synthetic rwl", {
  # Two 50-year series with slight correlation
  set.seed(7)
  signal <- cumsum(rnorm(50))
  m <- matrix(c(signal + rnorm(50, sd = 0.2),
                signal + rnorm(50, sd = 0.2)),
               nrow = 50, ncol = 2,
               dimnames = list(as.character(1950:1999), c("SYN1","SYN2")))
  # Scale to mm-like values
  m <- abs(m / max(abs(m)) * 2) + 0.5
  rwl <- structure(as.data.frame(m), class = c("rwl","data.frame"))

  cof <- dpl_cof(rwl, verbose = FALSE, parts = c(1L, 7L))

  # Return structure
  expect_true(is.list(cof))
  expect_true(!is.null(cof$master))
  expect_true(!is.null(cof$stats))
  expect_true(!is.null(cof$segments))
  expect_true(!is.null(cof$problems))

  # Stats table has one row per series
  expect_equal(nrow(cof$stats), 2L)
  expect_true(all(c("series","r_master","n_flags") %in% names(cof$stats)))

  # master is named and same length as year span
  expect_named(cof$master)
  expect_equal(length(cof$master), 50L)

  # crit is a reasonable value
  expect_gt(cof$crit, 0)
  expect_lt(cof$crit, 1)
})

test_that("dpl_cof rejects fewer than two series", {
  m <- matrix(runif(10), nrow = 10, ncol = 1,
               dimnames = list(as.character(1900:1909), "ONLY"))
  rwl <- structure(as.data.frame(m), class = c("rwl","data.frame"))
  expect_error(dpl_cof(rwl, verbose = FALSE), "at least two")
})
