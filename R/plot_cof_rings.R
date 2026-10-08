
# =============================================================================
# plot_cof_rings  --  ring-width diagram of a collection (TSAP-style bars)
# =============================================================================

#' Ring-width diagram of a collection: one bar per core, one cell per ring
#'
#' @description
#' Draws the layout familiar from TSAP-Win: each series is a horizontal bar
#' spanning its calendar years, divided into one cell per ring. Cells are
#' shaded by ring width relative to the series itself (darker = narrower), so
#' that narrow rings shared by the stand appear as dark vertical stripes and a
#' misdated core shows up as a bar whose stripes are offset from its
#' neighbours. Years in which a large share of the cores formed a narrow ring
#' are marked along the top as pointer years.
#'
#' @details
#' \strong{Shading.} With \code{shade = "decile"} (default) each ring is placed
#' in a decile of its own series, exactly as COFECHA's Part 4 bar plot
#' (\code{\link{dpl_barplot}}) does for the master: the 1st decile is the
#' darkest cell, the 10th the lightest. The shading is therefore a
#' within-series rank and does not depend on absolute growth rate or age
#' trend, which makes bars of fast- and slow-growing trees comparable.
#' \code{shade = "index"} uses the COFECHA-transformed series
#' (\code{cof_result$filtered}: detrended, log, prewhitened) and shades by
#' the index value clipped to +/- 2.5 SD; \code{shade = "none"} draws only the
#' spans.
#'
#' \strong{Pointer years.} A year is marked when at least \code{pointer_min}
#' series cover it and at least \code{pointer_frac} of them have the ring in
#' their lowest \code{pointer_decile} deciles (narrow pointer years, triangle
#' pointing down) or highest deciles (wide, triangle pointing up). These are
#' event years in the sense of Schweingruber et al. (1990), computed on
#' within-series deciles rather than on a fixed percentage change.
#'
#' \strong{Marked years.} \code{mark_years} draws thin vertical rules through
#' all bars at the given years: a known fire, a frost year, a candidate
#' dating problem. A named integer vector labels them.
#'
#' @param rwl A dplR \code{rwl} data.frame.
#' @param cof_result The list returned by \code{\link{dpl_cof}} on \code{rwl},
#'   needed for \code{shade = "index"} (its \code{$filtered} is used); run
#'   automatically when \code{NULL}.
#' @param series Series to draw: character IDs or integer positions. Default
#'   all.
#' @param sort \code{"none"} (default, file order), \code{"first"} (oldest
#'   inner ring on top), \code{"last"}, \code{"length"}, or \code{"id"}.
#' @param shade \code{"decile"} (default), \code{"index"} or \code{"none"}.
#' @param palette Character vector of colours, light to dark, for the
#'   shading ramp; or a single hue name among \code{"grey"} (default),
#'   \code{"blue"}, \code{"brown"}, \code{"green"}.
#' @param pointer Logical (default \code{TRUE}). Mark pointer years.
#' @param pointer_frac Numeric (default \code{0.6}). Fraction of covering
#'   series that must agree.
#' @param pointer_min Integer (default \code{5L}). Minimum covering series.
#' @param pointer_decile Integer (default \code{2L}). Deciles counted as
#'   narrow (1..\code{pointer_decile}) or wide.
#' @param mark_years Integer vector (optionally named) of years to rule.
#' @param tree_fun Function mapping series IDs to tree IDs, or \code{NULL}
#'   (default). When given, cores of the same tree are grouped and a thin
#'   bracket drawn beside their labels. E.g. \code{function(id) substr(id, 1,
#'   6)} for \code{ACC014A}/\code{ACC014B}.
#' @param years Integer vector of length 2, year range to draw. Default: the
#'   span of the selected series.
#' @param grid_by Integer. Spacing of year gridlines; \code{NULL} (default)
#'   picks 10, 50, 100 or 200 from the span.
#' @param bar_height Numeric in (0, 1] (default \code{0.72}). Fraction of the
#'   row occupied by the bar; the remainder is the gap between cores.
#' @param cell_border Logical or \code{NULL} (default). Draw a hairline
#'   between rings; \code{NULL} draws it only when the span is under 300
#'   years so that cells are wide enough to see.
#' @param main Title; \code{NULL} (default) builds one from the span.
#' @param legend Logical (default \code{TRUE}).
#' @param cex_id Label size (default \code{0.7}).
#' @param ... Passed to \code{\link[graphics]{plot.window}}; unused otherwise.
#'
#' @return Invisibly, a list with \code{deciles} (matrix years x series, the
#'   shading classes), \code{pointer} (data.frame of marked years with
#'   \code{year}, \code{type}, \code{n}, \code{frac}) and \code{order} (the
#'   series IDs top to bottom).
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' plot_cof_rings(rwl)
#' plot_cof_rings(rwl, sort = "first", tree_fun = function(id) substr(id, 1, 6))
#' plot_cof_rings(rwl, years = c(1850, 2002), palette = "blue",
#'                mark_years = c(Fire = 1860, 1912))
#'
#' # COFECHA-transformed indices instead of raw deciles
#' cof <- dpl_cof(rwl, parts = integer(0), verbose = FALSE)
#' plot_cof_rings(rwl, cof, shade = "index", years = c(1900, 2002))
#'
#' pdf("MIR_rings.pdf", width = 11.7, height = 8.3)   # A4 landscape
#' plot_cof_rings(rwl, sort = "first")
#' dev.off()
#' }
#'
#' @seealso \code{\link{dpl_barplot}}, \code{\link{plot_cof_barplot}},
#'   \code{\link{dpl_cof}}
#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom graphics par plot.new plot.window rect axis mtext text segments
#'   points abline strwidth legend
#' @export
plot_cof_rings <- function(rwl,
                           cof_result     = NULL,
                           series         = NULL,
                           sort           = c("none", "first", "last", "length", "id"),
                           shade          = c("decile", "index", "none"),
                           palette        = "grey",
                           pointer        = TRUE,
                           pointer_frac   = 0.6,
                           pointer_min    = 5L,
                           pointer_decile = 2L,
                           mark_years     = NULL,
                           tree_fun       = NULL,
                           years          = NULL,
                           grid_by        = NULL,
                           bar_height     = 0.72,
                           cell_border    = NULL,
                           main           = NULL,
                           legend         = TRUE,
                           cex_id         = 0.7,
                           ...) {

  sort  <- match.arg(sort)
  shade <- match.arg(shade)

  # ---- data ------------------------------------------------------------------
  if (!.is_rwl(rwl)) stop("'rwl' must be a dplR rwl data.frame.")
  idx <- if (shade == "index") {
    if (is.null(cof_result)) cof_result <- dpl_cof(rwl, parts = integer(0), verbose = FALSE)
    if (is.null(cof_result$filtered)) stop("'cof_result' has no $filtered element; re-run dpl_cof().")
    cof_result$filtered
  } else NULL

  ids <- .resolve_series(series, colnames(rwl))
  if (!length(ids)) stop("No series selected.")
  yrs <- as.integer(rownames(rwl))
  dat <- rwl[, ids, drop = FALSE]
  first <- vapply(dat, function(v) yrs[which(!is.na(v))[1L]], integer(1))
  last  <- vapply(dat, function(v) max(yrs[!is.na(v)]), integer(1))
  ord <- switch(sort,
    none   = seq_along(ids),
    first  = order(first, last),
    last   = order(-last, first),
    length = order(-(last - first)),
    id     = order(ids))
  if (!is.null(tree_fun)) {                     # keep cores of a tree together
    tree <- vapply(ids, tree_fun, character(1))
    ord  <- ord[order(match(tree[ord], unique(tree[ord])))]
  }
  ids <- ids[ord]; dat <- dat[, ord, drop = FALSE]
  first <- first[ord]; last <- last[ord]
  n <- length(ids)

  if (is.null(years)) years <- c(min(first), max(last))
  years <- as.integer(years)
  rows_in <- yrs >= years[1L] & yrs <= years[2L]

  # ---- shading classes ------------------------------------------------------
  # decile: rank within series over its WHOLE length (not just the window),
  #         so the shading is a property of the core, not of the view
  nclass <- 10L
  cls <- matrix(NA_integer_, nrow(dat), n, dimnames = list(rownames(dat), ids))
  if (shade == "decile") {
    for (j in seq_len(n)) {
      v <- dat[[j]]; ok <- !is.na(v)
      if (sum(ok) >= 2L) {
        r <- rank(v[ok], ties.method = "average")
        cls[ok, j] <- pmin(nclass, as.integer(ceiling(r / sum(ok) * nclass)))
      } else cls[ok, j] <- 5L
    }
  } else if (shade == "index") {
    for (j in seq_len(n)) {
      v <- idx[rownames(dat), ids[j]]; ok <- !is.na(v)
      z <- pmax(-2.5, pmin(2.5, v[ok]))
      cls[ok, j] <- pmin(nclass, as.integer(floor((z + 2.5) / 5 * nclass)) + 1L)
    }
  }

  # ---- colours ----------------------------------------------------------------
  # sequential single hue, light = wide, dark = narrow; row surface is paper
  hues <- list(
    grey  = c("#F2F2F0", "#1B1B1B"),
    blue  = c("#E8EEF5", "#102A4C"),
    brown = c("#F4ECE2", "#4A2C12"),
    green = c("#E9F0E6", "#143D1E"))
  if (length(palette) == 1L && palette %in% names(hues)) palette <- hues[[palette]]
  ramp <- grDevices::colorRampPalette(palette)(nclass)   # [1] light .. [n] dark
  col_class <- rev(ramp)                                  # class 1 (narrow) = dark
  col_span   <- "#D9D9D6"
  col_border <- "#FFFFFF"
  col_grid   <- "#C8C8C4"
  col_text   <- "#333333"
  col_muted  <- "#777777"
  col_narrow <- "#B3261E"
  col_wide   <- "#2E6E9E"
  col_mark   <- "#B3261E"

  # ---- pointer years ---------------------------------------------------------
  ptr <- data.frame(year = integer(0), type = character(0), n = integer(0),
                    frac = numeric(0), stringsAsFactors = FALSE)
  if (isTRUE(pointer) && shade != "none") {
    cover  <- rowSums(!is.na(cls))
    narrow <- rowSums(cls <= pointer_decile, na.rm = TRUE)
    wide   <- rowSums(cls >  nclass - pointer_decile, na.rm = TRUE)
    okc <- cover >= pointer_min
    fn <- ifelse(okc, narrow / pmax(1L, cover), 0)
    fw <- ifelse(okc, wide   / pmax(1L, cover), 0)
    yn <- yrs[fn >= pointer_frac]; yw <- yrs[fw >= pointer_frac]
    ptr <- rbind(
      data.frame(year = yn, type = rep("narrow", length(yn)), n = cover[fn >= pointer_frac],
                 frac = fn[fn >= pointer_frac], stringsAsFactors = FALSE),
      data.frame(year = yw, type = rep("wide", length(yw)), n = cover[fw >= pointer_frac],
                 frac = fw[fw >= pointer_frac], stringsAsFactors = FALSE))
    ptr <- ptr[ptr$year >= years[1L] & ptr$year <= years[2L], , drop = FALSE]
    ptr <- ptr[order(ptr$year), , drop = FALSE]; rownames(ptr) <- NULL
  }

  # ---- layout -----------------------------------------------------------------
  span <- diff(years) + 1L
  if (is.null(grid_by)) grid_by <- if (span <= 120) 10L else if (span <= 400) 50L else if (span <= 1200) 100L else 200L
  if (is.null(cell_border)) cell_border <- span < 300L
  lab_w <- max(graphics::strwidth(ids, units = "inches", cex = cex_id)) + 0.25
  tree_w <- if (is.null(tree_fun)) 0 else 0.12
  op <- graphics::par(mai = c(0.45, lab_w + tree_w + 0.1, 0.5, 0.25),
                      xaxs = "i", yaxs = "i", family = graphics::par("family"))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(years[1L] - 0.5, years[2L] + 0.5),
                        ylim = c(n + 0.5, 0.5), ...)          # row 1 on top

  # grid
  gy <- seq(ceiling(years[1L] / grid_by) * grid_by, years[2L], by = grid_by)
  graphics::abline(v = gy - 0.5, col = col_grid, lwd = 0.5, lty = "13")

  # bars
  h <- bar_height / 2
  for (j in seq_len(n)) {
    # span background (also what shade = "none" shows)
    graphics::rect(max(first[j], years[1L]) - 0.5, j - h,
                   min(last[j],  years[2L]) + 0.5, j + h,
                   col = col_span, border = NA)
    if (shade != "none") {
      k  <- which(rows_in & !is.na(cls[, j]))
      if (length(k)) {
        graphics::rect(yrs[k] - 0.5, j - h, yrs[k] + 0.5, j + h,
                       col = col_class[cls[k, j]],
                       border = if (cell_border) col_border else NA,
                       lwd = 0.3)
      }
    }
  }

  # marked years
  if (length(mark_years)) {
    my <- as.integer(mark_years); my <- my[my >= years[1L] & my <= years[2L]]
    graphics::abline(v = my, col = grDevices::adjustcolor(col_mark, 0.7), lwd = 0.8)
    if (!is.null(names(mark_years))) {
      lab <- names(mark_years)[as.integer(mark_years) %in% my]
      lab[is.na(lab)] <- ""
      graphics::mtext(lab, side = 3, at = my, line = 0.1, cex = 0.6, col = col_mark, adj = 0)
    }
  }

  # pointer triangles above the top row
  if (nrow(ptr)) {
    nn <- ptr$type == "narrow"
    graphics::points(ptr$year[nn], rep(0.5, sum(nn)) - 0.03 * n, pch = 25, cex = 0.55,
                     col = col_narrow, bg = col_narrow, xpd = NA)
    graphics::points(ptr$year[!nn], rep(0.5, sum(!nn)) - 0.03 * n, pch = 24, cex = 0.55,
                     col = col_wide, bg = col_wide, xpd = NA)
  }

  # axes
  graphics::axis(1, at = gy, labels = gy, col = NA, col.ticks = col_grid, tck = -0.012,
                 cex.axis = 0.7, col.axis = col_text, mgp = c(2, 0.3, 0), lwd.ticks = 0.5)
  graphics::axis(2, at = seq_len(n), labels = ids, las = 1, col = NA, tick = FALSE,
                 cex.axis = cex_id, col.axis = col_text, mgp = c(2, 0.2, 0), hadj = 1)
  if (!is.null(tree_fun)) {
    tree <- vapply(ids, tree_fun, character(1))
    grp  <- rle(tree)
    ends <- cumsum(grp$lengths); starts <- ends - grp$lengths + 1L
    xl <- graphics::grconvertX(0.08, from = "inches", to = "user")
    for (g in seq_along(starts)) if (grp$lengths[g] > 1L)
      graphics::segments(xl, starts[g] - 0.3, xl, ends[g] + 0.3, col = col_muted,
                         lwd = 1.2, xpd = NA, lend = 1)
  }

  # title and sample depth line
  if (is.null(main)) main <- sprintf("%d series, %d-%d", n, years[1L], years[2L])
  graphics::mtext(main, side = 3, line = 1.7, adj = 0, cex = 0.9, col = col_text, font = 2)

  # legend: narrow -> wide ramp plus pointer symbols, top right beside the title
  if (isTRUE(legend) && shade != "none") {
    usr <- graphics::par("usr")
    top_in <- graphics::grconvertY(usr[4L], from = "user", to = "inches")
    y0 <- graphics::grconvertY(top_in + 0.30, from = "inches", to = "user")
    y1 <- graphics::grconvertY(top_in + 0.40, from = "inches", to = "user")
    ym <- (y0 + y1) / 2
    w  <- 0.009 * span
    ptr_lab <- if (nrow(ptr)) "pointer year" else ""
    lx <- usr[2L] - nclass * w - (if (nrow(ptr)) 0.09 * span else 0.035 * span)
    for (i in seq_len(nclass))
      graphics::rect(lx + (i - 1L) * w, y0, lx + i * w, y1,
                     col = col_class[i], border = NA, xpd = NA)
    graphics::text(lx - 0.004 * span, ym, "narrow", adj = c(1, 0.5), cex = 0.6,
                   col = col_muted, xpd = NA)
    graphics::text(lx + nclass * w + 0.004 * span, ym, "wide", adj = c(0, 0.5),
                   cex = 0.6, col = col_muted, xpd = NA)
    if (nrow(ptr)) {
      px <- lx + nclass * w + 0.04 * span
      graphics::points(px, ym, pch = 25, cex = 0.55, col = col_narrow, bg = col_narrow, xpd = NA)
      graphics::points(px + 0.008 * span, ym, pch = 24, cex = 0.55, col = col_wide, bg = col_wide, xpd = NA)
      graphics::text(px + 0.015 * span, ym, ptr_lab, adj = c(0, 0.5), cex = 0.6,
                     col = col_muted, xpd = NA)
    }
  }

  invisible(list(deciles = cls, pointer = ptr, order = ids))
}
