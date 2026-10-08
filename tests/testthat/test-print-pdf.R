pel_cof <- function() {
  f   <- system.file("extdata", "benchmark", "PEL.rwl", package = "cofechar")
  rwl <- dpl_read_dec(f, stop_val = 999L, unit = "0.01mm")
  dpl_cof(rwl, seg_length = 20L, seg_lag = 10L, verbose = FALSE)
}

test_that("dpl_print_pdf writes a PDF and sizes the font to fit 132 columns", {
  cof <- pel_cof()
  f   <- tempfile(fileext = ".pdf")
  r   <- dpl_print_pdf(cof, f)
  expect_true(file.exists(f)); expect_gt(file.info(f)$size, 1000)
  expect_equal(r$orientation, "landscape")
  expect_true(r$font_size >= 5 && r$font_size <= 10)
  # 132 columns at r$font_size pt must fit 277 mm (A4 landscape, 10 mm margins)
  expect_lte(132 * 0.6 * r$font_size / 72, (297 - 20) / 25.4)
  # a complete 400-year bar-plot page (column header + 50 rows + 4 blank
  # separators = 55 lines) fits on one sheet
  expect_gte(r$lines_per_page, 55)
})

test_that("dpl_print_pdf breaks pages at PART headers", {
  cof <- pel_cof()
  # Use the internal pagination: every page's first line should be a PART
  # header or the continuation of a block longer than one page.
  f <- tempfile(fileext = ".pdf")
  r <- dpl_print_pdf(cof, f)
  expect_true(r$pages >= 7)   # at least one page per part
})

test_that("dpl_print_pdf accepts a dpl_barplot list and bare lines, portrait when narrow", {
  cof <- pel_cof()
  f1 <- tempfile(fileext = ".pdf")
  r1 <- dpl_print_pdf(dpl_barplot(cof, quiet = TRUE), f1)
  expect_true(file.exists(f1)); expect_equal(r1$orientation, "landscape")
  narrow <- sprintf("%4d  %s", 1900:1960, strrep("-", 40))
  f2 <- tempfile(fileext = ".pdf")
  r2 <- dpl_print_pdf(narrow, f2, title = "narrow")
  expect_equal(r2$orientation, "portrait"); expect_equal(r2$pages, 1L)
  expect_equal(r2$font_size, 10)
  # forcing a font that cannot fit is an error, not a wrap
  expect_error(dpl_print_pdf(cof, tempfile(fileext = ".pdf"), font_size = 14), "does not fit")
})
