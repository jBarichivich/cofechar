
# =============================================================================
# plot_cof_coverage  --  timelines of the series (sample coverage)
# =============================================================================

#' Timelines of the series: sample coverage by calendar year or by ring age
#'
#' @description
#' Draws one thin horizontal line per series spanning its years, stacked in
#' the order chosen, with the sample depth drawn underneath. On the calendar
#' axis this is the classical "timelines of the included series" figure of a
#' chronology; aligned at the first ring it becomes the regional-curve
#' standardisation (RCS) view of the collection, where the x axis is ring
#' number from the pith (cambial age) and the stack shows how many series
#' support each age; aligned at the last ring it shows how many series reach
#' back a given number of years from the bark.
#'
#' @details
#' With \code{align = "first"} the x value of a ring is its ring number
#' counted from the innermost ring present, plus \code{pith_offset} for that
#' series when given, so series with a known distance to pith can be placed
#' at their true cambial age. With \code{align = "last"} the x value is the
#' number of years before the outermost ring (0 at the bark, growing to the
#' right by default, or leftwards with \code{reverse = TRUE}).
#'
#' The sample-depth panel is the count of series present at each x; on the
#' cambial-age axis this is the replication of the regional curve and the
#' usual argument for truncating it where fewer than, say, 5 or 10 series
#' remain (\code{depth_min} draws that threshold).
#'
#' @param rwl A dplR \code{rwl} data.frame.
#' @param series Series to draw: character IDs or integer positions. Default
#'   all.
#' @param align \code{"year"} (default), \code{"first"} or \code{"last"}.
#' @param sort Order of the stack from bottom to top: \code{"first"}
#'   (default; earliest start at the bottom, as in the classical figure),
#'   \code{"last"}, \code{"length"} (shortest at the bottom), \code{"none"}
#'   (file order) or \code{"id"}.
#' @param pith_offset Named numeric vector (series ID = rings missing to the
#'   pith) or a vector in the order of \code{series}; only used with
#'   \code{align = "first"}. Missing entries count as 0.
#' @param reverse Logical (default \code{FALSE}). With \code{align = "last"},
#'   put the bark at the right and count leftwards.
#' @param depth Logical (default \code{TRUE}). Draw the sample-depth panel.
#' @param depth_min Integer or \code{NULL}. Horizontal reference line in the
#'   depth panel (e.g. \code{5}); \code{NULL} draws none.
#' @param highlight Series IDs or positions drawn in the accent colour, on
#'   top of the others.
#' @param col Line colour for the series (default \code{"#3A6EA5"}).
#' @param col_highlight Colour for \code{highlight} (default \code{"#B3261E"}).
#' @param lwd Line width (default \code{NULL}: chosen from the number of
#'   series so that lines neither overlap nor leave gaps).
#' @param labels Logical (default \code{FALSE}). Print the series IDs at the
#'   left of the lines; useful up to a few dozen series.
#' @param xlab,main Axis label and title; \code{NULL} builds them from
#'   \code{align}.
#' @param cex_id Label size when \code{labels = TRUE}.
#'
#' @return Invisibly, a list with \code{series} (data.frame: \code{series},
#'   \code{first}, \code{last}, \code{n}, \code{x0}, \code{x1}, \code{row})
#'   and \code{depth} (data.frame: \code{x}, \code{n}).
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' plot_cof_coverage(rwl)                                   # calendar timelines
#' plot_cof_coverage(rwl, align = "first", depth_min = 5)   # RCS view: cambial age
#' plot_cof_coverage(rwl, align = "last", reverse = TRUE)   # years before bark
#' plot_cof_coverage(rwl, labels = TRUE, highlight = c("ACC026B", "ACC026"))
#'
#' # with known distances to pith for some cores
#' po <- c(ACC026B = 12, ACC026 = 40)
#' plot_cof_coverage(rwl, align = "first", pith_offset = po)
#' }
#'
#' @seealso \code{\link{plot_cof_rings}}, \code{\link{plot_cof_barplot}},
#'   \code{\link{dpl_cof}}
#' @importFrom graphics layout polygon lines
#' @export
plot_cof_coverage <- function(rwl,
                              series        = NULL,
                              align         = c("year", "first", "last"),
                              sort          = c("first", "last", "length", "none", "id"),
                              pith_offset   = NULL,
                              reverse       = FALSE,
                              depth         = TRUE,
                              depth_min     = NULL,
                              highlight     = NULL,
                              col           = "#3A6EA5",
                              col_highlight = "#B3261E",
                              lwd           = NULL,
                              labels        = FALSE,
                              xlab          = NULL,
                              main          = NULL,
                              cex_id        = 0.6) {

  align <- match.arg(align)
  sort  <- match.arg(sort)
  if (!.is_rwl(rwl)) stop("'rwl' must be a dplR rwl data.frame.")
  ids <- .resolve_series(series, colnames(rwl))
  if (!length(ids)) stop("No series selected.")
  yrs <- as.integer(rownames(rwl))
  dat <- rwl[, ids, drop = FALSE]
  n   <- length(ids)

  first <- vapply(dat, function(v) yrs[which(!is.na(v))[1L]], integer(1))
  last  <- vapply(dat, function(v) max(yrs[!is.na(v)]), integer(1))
  nyr   <- last - first + 1L

  # ---- x position of each series ---------------------------------------------
  po <- rep(0, n); names(po) <- ids
  if (!is.null(pith_offset) && align == "first") {
    if (!is.null(names(pith_offset))) {
      m <- match(names(pith_offset), ids)
      po[m[!is.na(m)]] <- pith_offset[!is.na(m)]
    } else if (length(pith_offset) == n) po[] <- pith_offset
    else stop("'pith_offset' must be named or have one value per series.")
  }
  x0 <- switch(align,
    year  = first,
    first = 1 + po,
    last  = if (reverse) -(nyr - 1L) else 0)
  x1 <- switch(align,
    year  = last,
    first = nyr + po,
    last  = if (reverse) 0 else nyr - 1L)
  if (align == "last" && !reverse) { x0 <- rep(0, n); x1 <- nyr - 1 }
  if (align == "last" &&  reverse) { x0 <- -(nyr - 1); x1 <- rep(0, n) }

  # ---- stacking order (row 1 at the bottom) ------------------------------------
  ord <- switch(sort,
    first  = if (align == "last") order(-nyr, first) else order(x0, x1),
    last   = order(x1, x0),
    length = order(nyr, first),
    none   = seq_len(n),
    id     = order(ids))
  row <- integer(n); row[ord] <- seq_len(n)

  # ---- sample depth -------------------------------------------------------------
  xr <- range(c(x0, x1))
  xs <- seq(floor(xr[1L]), ceiling(xr[2L]))
  dep <- integer(length(xs))
  for (j in seq_len(n)) {
    k <- xs >= x0[j] & xs <= x1[j]
    # calendar axis: count only years actually measured (gaps inside a series)
    if (align == "year") k <- k & !is.na(dat[[j]][match(xs, yrs)])
    dep[k] <- dep[k] + 1L
  }

  # ---- colours and widths --------------------------------------------------------
  hl <- if (is.null(highlight)) character(0) else .resolve_series(highlight, ids)
  col_grid <- "#C8C8C4"; col_text <- "#333333"; col_muted <- "#777777"
  col_depth <- "#D5DCE6"

  # ---- layout ------------------------------------------------------------------
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (isTRUE(depth)) graphics::layout(matrix(1:2, 2L), heights = c(3.2, 1))
  lab_in <- if (labels) max(graphics::strwidth(ids, units = "inches", cex = cex_id)) + 0.15 else 0
  mai_l <- 0.55 + lab_in
  graphics::par(mai = c(if (depth) 0.1 else 0.5, mai_l, 0.5, 0.25), xaxs = "i", yaxs = "i")

  if (is.null(xlab)) xlab <- switch(align,
    year  = "Year",
    first = if (any(po > 0)) "Cambial age (rings from pith)" else "Ring number from first ring",
    last  = "Years before last ring")
  if (is.null(main)) main <- switch(align,
    year  = sprintf("%d series, %d-%d", n, min(first), max(last)),
    first = sprintf("%d series aligned at the first ring", n),
    last  = sprintf("%d series aligned at the last ring", n))

  pad <- 0.01 * diff(xr)
  xlim <- c(xr[1L] - pad, xr[2L] + pad)
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = c(0, n + 1))
  gx <- pretty(xlim, n = 8)
  graphics::abline(v = gx, col = col_grid, lwd = 0.5, lty = "13")
  if (!labels) graphics::abline(h = pretty(c(0, n), n = 5), col = col_grid, lwd = 0.5, lty = "13")

  # line width: fill the row height, leaving ~25% gap, capped for small n
  if (is.null(lwd)) {
    row_in <- (graphics::par("pin")[2L]) / (n + 1)
    lwd <- max(0.4, min(3, row_in * 72 * 0.75))
  }
  main_col <- ifelse(ids %in% hl, NA, col)
  k <- which(!is.na(main_col))
  graphics::segments(x0[k], row[k], x1[k], row[k], col = col, lwd = lwd, lend = 1)
  if (length(hl)) {
    k <- which(ids %in% hl)
    graphics::segments(x0[k], row[k], x1[k], row[k], col = col_highlight, lwd = lwd, lend = 1)
  }
  if (labels)
    graphics::text(graphics::par("usr")[1L], row, ids, adj = c(1.05, 0.5), cex = cex_id,
                   col = ifelse(ids %in% hl, col_highlight, col_text), xpd = NA)

  if (!labels) {
    graphics::axis(2, at = pretty(c(0, n), n = 5), las = 1, col = NA, col.ticks = col_grid,
                   tck = -0.012, cex.axis = 0.7, col.axis = col_text, mgp = c(2, 0.4, 0))
    graphics::mtext("Series", side = 2, line = 2.0, cex = 0.75, col = col_muted)
  }
  if (!depth) {
    graphics::axis(1, at = gx, col = NA, col.ticks = col_grid, tck = -0.012, cex.axis = 0.7,
                   col.axis = col_text, mgp = c(2, 0.3, 0))
    graphics::mtext(xlab, side = 1, line = 1.6, cex = 0.75, col = col_muted)
  } else {
    graphics::axis(1, at = gx, labels = FALSE, col = NA, col.ticks = col_grid, tck = -0.012)
  }
  graphics::mtext(main, side = 3, line = 0.6, adj = 0, cex = 0.9, col = col_text, font = 2)

  # ---- depth panel -------------------------------------------------------------
  if (isTRUE(depth)) {
    graphics::par(mai = c(0.5, mai_l, 0.05, 0.25))
    graphics::plot.new()
    ymax <- max(dep) * 1.08
    graphics::plot.window(xlim = xlim, ylim = c(0, ymax))
    graphics::abline(v = gx, col = col_grid, lwd = 0.5, lty = "13")
    graphics::polygon(c(xs[1L] - 0.5, rep(xs, each = 2) + c(-0.5, 0.5), xs[length(xs)] + 0.5),
                      c(0, rep(dep, each = 2), 0), col = col_depth, border = NA)
    graphics::lines(rep(xs, each = 2) + c(-0.5, 0.5), rep(dep, each = 2), col = col, lwd = 1)
    if (!is.null(depth_min)) {
      graphics::abline(h = depth_min, col = col_highlight, lwd = 0.8, lty = "22")
      graphics::text(xlim[2L], depth_min, sprintf("n = %d", as.integer(depth_min)),
                     adj = c(1.05, -0.3), cex = 0.6, col = col_highlight, xpd = NA)
    }
    yt <- pretty(c(0, max(dep)), n = 3)
    graphics::axis(2, at = yt, las = 1, col = NA, col.ticks = col_grid, tck = -0.03,
                   cex.axis = 0.7, col.axis = col_text, mgp = c(2, 0.4, 0))
    graphics::axis(1, at = gx, col = NA, col.ticks = col_grid, tck = -0.03, cex.axis = 0.7,
                   col.axis = col_text, mgp = c(2, 0.3, 0))
    graphics::mtext("Sample depth", side = 2, line = 2.0, cex = 0.7, col = col_muted)
    graphics::mtext(xlab, side = 1, line = 1.6, cex = 0.75, col = col_muted)
  }

  invisible(list(
    series = data.frame(series = ids, first = first, last = last, n = nyr,
                        x0 = x0, x1 = x1, row = row, stringsAsFactors = FALSE),
    depth  = data.frame(x = xs, n = dep)))
}
