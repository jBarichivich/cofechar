test_that("plot_cof_coverage stacks series and counts sample depth", {
  f <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
  rwl <- dpl_read_dec(f, label_length = NULL, stop_val = -9999L, unit = "0.001mm")
  pdf(NULL); on.exit(dev.off())
  r <- plot_cof_coverage(rwl, highlight = "ACC026B")
  expect_equal(nrow(r$series), 36L)
  expect_equal(r$series$row[r$series$series == "ACC026B"], 1L)      # earliest at bottom
  expect_equal(r$depth$n[r$depth$x == 1950], sum(!is.na(rwl["1950", ])))
  r2 <- plot_cof_coverage(rwl, align = "first", depth_min = 5, labels = TRUE,
                          pith_offset = c(ACC026B = 10))
  expect_equal(r2$series$x0[r2$series$series == "ACC026B"], 11)
  expect_equal(max(r2$depth$n), 36L)
  r3 <- plot_cof_coverage(rwl, align = "last", reverse = TRUE, depth = FALSE)
  expect_true(all(r3$series$x1 == 0))
  expect_equal(r3$depth$n[r3$depth$x == 0], 36L)
})
