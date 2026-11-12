test_that("plots render accents and restore LC_CTYPE in a C/POSIX session", {
  skip_on_os("windows")
  old <- Sys.getlocale("LC_CTYPE")
  on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE", old)), add = TRUE)
  skip_if_not(nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", "C"))), "no C locale")
  # the fix needs some UTF-8 locale on the machine to switch to
  probe <- function() { sw <- .local_utf8_ctype(); sw }
  skip_if_not(isTRUE(probe()), "no UTF-8 locale available to switch to")
  expect_equal(Sys.getlocale("LC_CTYPE"), "C")            # probe restored it

  f <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
  invisible(capture.output(rwl <- dpl_read_dec(f, stop_val = -9999L, unit = "0.001mm")))
  ttl <- "Región de Los Ríos, Laña"         # Latin-1 range only
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_warning(plot_cof_coverage(rwl[, 1:4], main = ttl, depth = FALSE))
  expect_no_warning(plot_cof_rings(rwl[, 1:4], main = ttl))
  expect_equal(Sys.getlocale("LC_CTYPE"), "C")            # session left as found
})
