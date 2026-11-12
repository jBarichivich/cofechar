# Lines are built with sprintf so column layout is explicit.
dec <- function(id, yr, vals, w = 8, sep = TRUE) {
  lab <- if (sep) sprintf("%-*s%5d", w, id, yr) else paste0(id, yr)
  paste0(lab, paste(sprintf("%6d", vals), collapse = ""))
}
read_q <- function(...) {
  r <- NULL
  invisible(capture.output(r <- suppressMessages(dpl_read_dec(...))))
  r
}

test_that("spaced format with negative (BCE) years is read literally", {
  lines <- c(dec("AMI02", -446, 101:104),           # -446..-443 (partial decade)
             dec("AMI02", -440, 111:120),
             dec("AMI02", -430, c(121:122, 999)))
  r <- read_q(lines, stop_val = 999L, unit = "mm")
  expect_equal(range(as.integer(rownames(r))), c(-446L, -429L))
  expect_equal(r["-446", "AMI02"], 101)
  expect_equal(r["-440", "AMI02"], 111)
  expect_equal(r["-429", "AMI02"], 122)
})

test_that("compact format with the negative year glued to the ID", {
  lines <- c(paste0("AMI02XYZ-446", paste(sprintf("%6d", 1:4), collapse = "")),
             paste0("AMI02XYZ-440", paste(sprintf("%6d", 11:20), collapse = "")),
             paste0("AMI02XYZ-430", paste(sprintf("%6d", c(21, 999)), collapse = "")))
  r <- read_q(lines, stop_val = 999L, unit = "mm")
  expect_equal(colnames(r), "AMI02XYZ")
  expect_equal(range(as.integer(rownames(r))), c(-446L, -430L))
})

test_that("five-character years (< -999) eating into the label", {
  lines <- c(paste0("AMI02-1005", paste(sprintf("%6d", 1:5), collapse = "")),
             paste0("AMI02-1000", paste(sprintf("%6d", 11:20), collapse = "")),
             paste0("AMI02 -990", paste(sprintf("%6d", c(21, 999)), collapse = "")))
  r <- read_q(lines, stop_val = 999L, unit = "mm")
  expect_equal(colnames(r), "AMI02")
  expect_equal(r["-1005", "AMI02"], 1)
  expect_equal(r["-990", "AMI02"], 21)
})

test_that("a series crossing from BCE to CE keeps its years", {
  lines <- c(dec("X1", -20, 1:10), dec("X1", -10, 11:20),
             dec("X1", 0, 21:30), dec("X1", 10, c(31, 999)))
  r <- read_q(lines, stop_val = 999L, unit = "mm")
  expect_equal(as.integer(rownames(r)), -20:10)
  expect_equal(r["0", "X1"], 21)
})

test_that("CE years below 1000 are no longer shifted by base_century", {
  lines <- c(dec("MED1", 446, 1:4), dec("MED1", 450, c(5, 999)))
  r <- read_q(lines, stop_val = 999L, unit = "mm")
  expect_equal(range(as.integer(rownames(r))), c(446L, 450L))
})

test_that("century-relative years still work when base_century is given", {
  lines <- c("OLD1    0   1   2   3", "OLD1   10  11 999")
  r <- read_q(lines, stop_val = 999L, unit = "mm", base_century = 1900L)
  expect_equal(range(as.integer(rownames(r))), c(1900L, 1910L))
})

test_that("dashed IDs in spaced format are not mistaken for negative years", {
  lines <- c(dec("BAR-01", 1930, 1:10), dec("BAR-01", 1940, c(11, 999)),
             dec("BAR-02", 1930, 21:30), dec("BAR-02", 1940, c(31, 999)))
  r <- read_q(lines, stop_val = 999L, unit = "mm")
  expect_setequal(colnames(r), c("BAR-01", "BAR-02"))
  expect_equal(range(as.integer(rownames(r))), c(1930L, 1940L))
})

test_that("label_length resolves compact negative years explicitly", {
  lines <- c(paste0("AMI02", "-446", paste(sprintf("%6d", 1:4), collapse = "")),
             paste0("AMI02", "-440", paste(sprintf("%6d", c(5, 999)), collapse = "")),
             dec("B1", -446, c(7, 8, 999)))
  r <- read_q(lines, label_length = 5L, stop_val = 999L, unit = "mm")
  expect_setequal(colnames(r), c("AMI02", "B1"))
  expect_equal(r["-446", "AMI02"], 1)
  expect_equal(r["-446", "B1"], 7)
})

test_that("a BCE collection runs through dpl_cof, barplots and plots", {
  set.seed(42)
  yrs <- -446:-100
  sig <- rnorm(length(yrs))
  mk <- function(a, b) {
    k <- yrs >= a & yrs <= b
    round(100 * exp(0.3 * sig[k] + rnorm(sum(k), 0, 0.15)))
  }
  spans <- list(c(-446, -150), c(-420, -100), c(-400, -120), c(-380, -110), c(-350, -100))
  lines <- character(0)
  for (i in seq_along(spans)) {
    v <- mk(spans[[i]][1], spans[[i]][2]); y0 <- spans[[i]][1]; yy <- y0:spans[[i]][2]
    dec0 <- yy - (yy %% 10)
    for (d in unique(dec0)) {
      vv <- v[dec0 == d]
      if (d == max(dec0)) vv <- c(vv, 999)
      lines <- c(lines, dec(sprintf("AMI%02d", i), min(yy[dec0 == d]), vv))
    }
  }
  r <- read_q(lines, stop_val = 999L, unit = "0.01mm")
  expect_equal(min(as.integer(rownames(r))), -446L)
  cof <- dpl_cof(r, verbose = FALSE, min_length = 10L)
  expect_true(all(is.finite(cof$master)))
  expect_true(any(grepl("-446", cof$output)))
  expect_gt(mean(cof$stats$r_master), 0.5)
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot_cof_rings(r))
  expect_silent(plot_cof_coverage(r))
  bp <- dpl_barplot(r, series = 1, quiet = TRUE)
  expect_true(any(grepl("-44[0-9]", unlist(bp))))
})
