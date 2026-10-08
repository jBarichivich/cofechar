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

test_that("dpl_edt keep/drop select series by name or position before editing", {
  m <- matrix(runif(20), nrow = 4, ncol = 5,
               dimnames = list(as.character(1900:1903), c("A", "B", "C", "D", "E")))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))
  expect_equal(colnames(dpl_edt(rwl, keep = c("B", "D"), verbose = FALSE)), c("B", "D"))
  expect_equal(colnames(dpl_edt(rwl, keep = 1:3, verbose = FALSE)), c("A", "B", "C"))
  expect_equal(colnames(dpl_edt(rwl, keep = c(4, 2), verbose = FALSE)), c("B", "D"))  # input order
  expect_equal(colnames(dpl_edt(rwl, drop = "C", verbose = FALSE)), c("A", "B", "D", "E"))
  expect_equal(colnames(dpl_edt(rwl, drop = c(1, 5), verbose = FALSE)), c("B", "C", "D"))
  expect_error(dpl_edt(rwl, keep = "A", drop = "B", verbose = FALSE), "either")
  expect_error(dpl_edt(rwl, keep = 9, verbose = FALSE), "out of range")
  expect_warning(dpl_edt(rwl, keep = c("A", "ZZ"), verbose = FALSE), "not found")
  # edits apply to the kept series; positions in edits refer to the kept set
  out <- dpl_edt(rwl, keep = c("B", "D"),
                 edits = list(list(series = "D", op = "replace", year = 1901, value = 9)),
                 verbose = FALSE)
  expect_equal(out["1901", "D"], 9)
  expect_equal(colnames(out), c("B", "D"))
})

test_that("dpl_edit_file copies a subset of samples to a new file", {
  m <- matrix(round(runif(30, 0.5, 3), 2), nrow = 6, ncol = 5,
               dimnames = list(as.character(1990:1995), c("A", "B", "C", "D", "E")))
  rwl <- structure(as.data.frame(m), class = c("rwl", "data.frame"))
  src <- tempfile(fileext = ".rwl"); dst <- tempfile(fileext = ".rwl")
  dpl_write(rwl, src, format = "tucson")
  dpl_edit_file(src, output_path = dst, keep = c("A", "E"), format = "tucson",
                verbose = FALSE)
  back <- dpl_read(dst, format = "tucson")
  expect_equal(colnames(back), c("A", "E"))
  expect_equal(unname(unlist(back["1992", ])), unname(unlist(rwl["1992", c("A", "E")])))
  dpl_edit_file(src, output_path = dst, drop = 2:4, format = "tucson", verbose = FALSE)
  expect_equal(colnames(dpl_read(dst, format = "tucson")), c("A", "E"))
})

test_that("dpl_trim and dpl_merge(trim = TRUE) span only the selected series", {
  f   <- system.file("extdata", "benchmark", "PEL.rwl", package = "cofechar")
  rwl <- dpl_read_dec(f, stop_val = 999L, unit = "0.01mm")        # 1862-2022
  yrs <- function(x) range(as.integer(rownames(x)))
  one <- rwl[, "PEL14A", drop = FALSE]                            # [ ] keeps 1862
  expect_equal(yrs(one), c(1862L, 2022L))
  expect_equal(yrs(dpl_trim(one)), c(1968L, 2022L))
  expect_equal(dim(dpl_trim(one)), c(55L, 1L))
  expect_equal(yrs(dpl_trim(rwl, series = c("PEL14A", "PEL07A"))), c(1962L, 2022L))
  expect_equal(yrs(dpl_trim(rwl, series = 5)), c(1968L, 2022L))
  # keep/drop in dpl_edt already trim
  expect_equal(yrs(dpl_edt(rwl, keep = "PEL14A", verbose = FALSE)), c(1968L, 2022L))
  expect_equal(yrs(dpl_edt(rwl, drop = "PEL11A", verbose = FALSE)), c(1929L, 2022L))
  # dpl_merge of [ ]-subsets now trims by default; trim = FALSE keeps the union
  m <- dpl_merge(list(rwl[, "PEL14A", drop = FALSE], rwl[, "PEL07A", drop = FALSE]))
  expect_equal(yrs(m), c(1962L, 2022L))
  m2 <- dpl_merge(list(rwl[, "PEL14A", drop = FALSE], rwl[, "PEL07A", drop = FALSE]), trim = FALSE)
  expect_equal(yrs(m2), c(1862L, 2022L))
  # interior gaps are never removed; all-NA gives 0 rows
  g <- rwl[, "PEL06B", drop = FALSE]; g["1950", 1] <- NA
  expect_equal(nrow(dpl_trim(g)), 94L)
  expect_equal(nrow(dpl_trim(rwl[, "PEL06B", drop = FALSE][0, , drop = FALSE])), 0L)
  expect_s3_class(dpl_trim(one), "rwl")
})
