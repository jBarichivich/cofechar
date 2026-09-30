test_that("dpl_edt copy/omit work on a minimal rwl", {
  # Build a tiny rwl in-memory
  m <- matrix(c(1.2, 1.5, 1.1, 0.9,
                0.8, 1.3, 1.6, 1.0),
               nrow = 4, ncol = 2,
               dimnames = list(as.character(1900:1903), c("S1", "S2")))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))

  # copy both → same dimensions
  out <- dpl_edt(rwl, verbose = FALSE)
  expect_equal(dim(out), dim(rwl))
  expect_equal(colnames(out), c("S1", "S2"))

  # omit S2
  out2 <- dpl_edt(rwl,
                    edits = list(list(series = "S2", op = "omit")),
                    verbose = FALSE)
  expect_equal(colnames(out2), "S1")
})

test_that("dpl_edt replace changes exactly one value", {
  m <- matrix(c(1.2, 1.5, 1.1,
                0.8, 1.3, 1.6),
               nrow = 3, ncol = 2,
               dimnames = list(as.character(1900:1902), c("A", "B")))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))

  out <- dpl_edt(rwl,
                   edits = list(list(series = "A", op = "replace",
                                     year = 1901, value = 9.99)),
                   verbose = FALSE)
  # 9.99 → 10.0 (DPL absent-ring sentinel replacement)
  expect_equal(out["1901", "A"], 10.0)
  # Other values unchanged
  expect_equal(out["1900", "A"], 1.2)
  expect_equal(out["1902", "A"], 1.1)
})

test_that("dpl_edt insert increases series length by one", {
  m <- matrix(c(1.2, 1.5, 1.1), nrow = 3, ncol = 1,
               dimnames = list(as.character(1900:1902), "X"))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))

  out <- dpl_edt(rwl,
                   edits = list(list(series = "X", op = "insert",
                                     year = 1901, value = 0.5, move = "back")),
                   verbose = FALSE)
  ok <- !is.na(out[, "X"])
  expect_equal(sum(ok), 4L)          # one more ring
  expect_equal(as.integer(rownames(out)[which(ok)[1]]), 1899L)  # shifted back
})

test_that("dpl_edt delete decreases series length by one", {
  m <- matrix(c(1.2, 1.5, 1.1, 0.9), nrow = 4, ncol = 1,
               dimnames = list(as.character(1900:1903), "Y"))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))

  out <- dpl_edt(rwl,
                   edits = list(list(series = "Y", op = "delete",
                                     year = 1901, move = "forward")),
                   verbose = FALSE)
  ok <- !is.na(out[, "Y"])
  expect_equal(sum(ok), 3L)
})

test_that("dpl_edt trim_start and trim_end work", {
  m <- matrix(1:6 / 10, nrow = 6, ncol = 1,
               dimnames = list(as.character(1900:1905), "Z"))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))

  out <- dpl_edt(rwl,
                   edits = list(
                     list(series = "Z", op = "trim_start", first_year = 1902),
                     list(series = "Z", op = "trim_end",   last_year  = 1904)
                   ),
                   verbose = FALSE)
  ok <- !is.na(out[, "Z"])
  yrs <- as.integer(rownames(out)[ok])
  expect_equal(yrs, 1902:1904)
})

test_that("dpl_edt default_action = 'omit' keeps only explicitly copied series", {
  m <- matrix(runif(12), nrow = 4, ncol = 3,
               dimnames = list(as.character(1900:1903), c("A", "B", "C")))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))

  out <- dpl_edt(rwl,
                   edits = list(list(series = "B", op = "copy")),
                   default_action = "omit",
                   verbose = FALSE)
  expect_equal(colnames(out), "B")
})

test_that("dpl_merge joins two rwl on union year axis", {
  m1 <- matrix(c(1.1, 1.2), nrow = 2, ncol = 1,
                dimnames = list(c("1900","1901"), "S1"))
  m2 <- matrix(c(0.9, 1.3), nrow = 2, ncol = 1,
                dimnames = list(c("1901","1902"), "S2"))
  r1 <- structure(as.data.frame(m1), class = c("rwl","data.frame"))
  r2 <- structure(as.data.frame(m2), class = c("rwl","data.frame"))

  out <- dpl_merge(list(r1, r2))
  expect_equal(rownames(out), c("1900","1901","1902"))
  expect_equal(colnames(out), c("S1","S2"))
  expect_true(is.na(out["1902","S1"]))
  expect_true(is.na(out["1900","S2"]))
  expect_equal(out["1901","S1"], 1.2)
})

test_that("dpl_merge dup_action = 'error' fires on duplicate IDs", {
  m <- matrix(1, nrow = 2, ncol = 1,
               dimnames = list(c("1900","1901"), "DUP"))
  r <- structure(as.data.frame(m), class = c("rwl","data.frame"))
  expect_error(dpl_merge(list(r, r)), "Duplicate")
})
