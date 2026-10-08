test_that("plot_cof_rings draws and returns classes, pointer years and order", {
  f <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
  rwl <- dpl_read_dec(f, label_length = NULL, stop_val = -9999L,
                      unit = "0.001mm")
  pdf(NULL); on.exit(dev.off())
  r <- plot_cof_rings(rwl, sort = "first", tree_fun = function(id) substr(id, 1, 6))
  expect_equal(dim(r$deciles), dim(rwl))
  expect_true(all(r$deciles >= 1 & r$deciles <= 10, na.rm = TRUE))
  expect_equal(r$order[1], "ACC026B")             # oldest inner ring on top
  expect_true(1944 %in% r$pointer$year[r$pointer$type == "narrow"])
  r2 <- plot_cof_rings(rwl, series = 1:5, years = c(1900, 2000), shade = "none",
                       mark_years = c(Fire = 1960))
  expect_equal(r2$order, colnames(rwl)[1:5])
  expect_equal(nrow(r2$pointer), 0L)
  r3 <- plot_cof_rings(rwl, xaxis = "length", sort = "length", mark_years = 1944)
  expect_equal(r3$order[1], "ACC026B")
  r4 <- plot_cof_rings(rwl, xaxis = "length", align = "inner", years = c(1850, 2002))
  expect_equal(nrow(r4$pointer), sum(r$pointer$year >= 1850))
  expect_error(suppressWarnings(plot_cof_rings(rwl, series = "NOPE")), "No series")
})
