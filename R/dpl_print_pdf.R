
# =============================================================================
# dpl_print_pdf  --  print COFECHA text output and ASCII bar plots to PDF
# =============================================================================

#' Print COFECHA output or ASCII bar plots to a PDF for the printer
#'
#' @description
#' Lays out the fixed-width text produced by \code{\link{dpl_cof}}
#' (\code{$output}), \code{\link{dpl_barplot}}, \code{\link{dpl_short_barplot}}
#' or any character vector of lines on A4 or Letter pages in a monospaced
#' font, with the font size chosen so that the longest line fits the page
#' width without wrapping or truncation, and page breaks placed at the
#' natural block boundaries of the output (part headers, series blocks,
#' 400-year bar-plot pages) rather than in the middle of a block.
#'
#' @details
#' DPL output was written for a 132-column line printer. On A4 landscape with
#' 10 mm margins, 132 columns fit at 9 pt Courier; on A4 portrait at 7 pt.
#' \code{orientation = "auto"} (the default) therefore chooses landscape when
#' the longest line exceeds \code{portrait_max} characters and portrait
#' otherwise, so \code{\link{dpl_short_barplot}} output and the column
#' layout of \code{\link{dpl_barplot}} come out portrait and the Part 4 page
#' layout landscape. \code{font_size} overrides the computed size; if the
#' text then does not fit, the function stops rather than wrap.
#'
#' Page breaks are taken, in order of preference, at lines beginning
#' \code{"PART "}, at the \code{=====} rules that open a series block in
#' \code{dpl_barplot} and Part 6, and at the blank line that closes a
#' 400-year bar-plot page; a block longer than a page is split at a blank
#' line, or as a last resort at any line. Every page carries a running
#' header with the title and the page number.
#'
#' @param x Character vector of lines, or the list returned by
#'   \code{\link{dpl_cof}} (its \code{$output} is used) or by
#'   \code{\link{dpl_barplot}} (elements are concatenated).
#' @param file Output PDF path.
#' @param title Character. Running header on every page. Default: the file
#'   name without extension.
#' @param paper \code{"a4"} (default) or \code{"letter"}.
#' @param orientation \code{"auto"} (default), \code{"portrait"} or
#'   \code{"landscape"}.
#' @param portrait_max Integer (default \code{100L}). Under \code{"auto"},
#'   the widest text printed portrait.
#' @param font_size Numeric or \code{NULL} (default). Font size in points;
#'   computed from the longest line when \code{NULL}.
#' @param margin_mm Numeric (default \code{10}). Page margin.
#' @param family Font family (default \code{"Courier"}; any monospaced
#'   family known to \code{\link[grDevices]{pdf}}).
#' @param open Logical (default \code{FALSE}). Open the PDF after writing.
#'
#' @return Invisibly, a list with \code{file}, \code{pages}, \code{font_size},
#'   \code{orientation} and \code{lines_per_page}.
#'
#' @examples
#' \dontrun{
#' cof <- dpl_cof(rwl_all, output_file = "MAI.COF")
#' dpl_print_pdf(cof, "MAI_COF.pdf")                    # full 7-part output
#' dpl_print_pdf(cof$output[grep("^PART 4", cof$output)[1]:length(cof$output)],
#'               "MAI_part4.pdf")                        # from Part 4 on
#'
#' bp <- dpl_barplot(rwl_all, series = 1:8, quiet = TRUE)
#' dpl_print_pdf(bp, "MAI_barplots.pdf")                # one series block per page
#'
#' sb <- dpl_short_barplot(rwl_all, cof, series = 1:8, quiet = TRUE)
#' dpl_print_pdf(sb, "MAI_short.pdf", orientation = "portrait")
#' }
#'
#' @seealso \code{\link{dpl_barplot}}, \code{\link{dpl_short_barplot}},
#'   \code{\link{dpl_cof}}
#' @export
dpl_print_pdf <- function(x, file,
                          title        = NULL,
                          paper        = c("a4", "letter"),
                          orientation  = c("auto", "portrait", "landscape"),
                          portrait_max = 100L,
                          font_size    = NULL,
                          margin_mm    = 10,
                          family       = "Courier",
                          open         = FALSE) {

  paper       <- match.arg(paper)
  orientation <- match.arg(orientation)

  # ---- lines -----------------------------------------------------------------
  lines <- if (is.list(x) && !is.null(x$output)) x$output
           else if (is.list(x)) unlist(x, use.names = FALSE)
           else as.character(x)
  lines <- sub("\\s+$", "", lines)
  lines <- unlist(strsplit(lines, "\n", fixed = TRUE))   # embedded newlines
  if (length(lines) == 0L) stop("Nothing to print.")
  if (is.null(title)) title <- sub("\\.[Pp][Dd][Ff]$", "", basename(file))

  # ---- page geometry (inches) ---------------------------------------------
  dims <- switch(paper, a4 = c(8.27, 11.69), letter = c(8.5, 11))
  maxw <- max(nchar(lines))
  if (orientation == "auto")
    orientation <- if (maxw > portrait_max) "landscape" else "portrait"
  if (orientation == "landscape") dims <- rev(dims)
  mar  <- margin_mm / 25.4
  usable_w <- dims[1L] - 2 * mar
  usable_h <- dims[2L] - 2 * mar

  # ---- block structure ------------------------------------------------------
  # Hard boundaries: PART headers and the ===== rule opening a series block.
  # Bar-plot pages (8 cols x 50 rows + separators) must not be split, so the
  # longest run between hard boundaries that contains a bar-plot column
  # header also constrains the font height.
  n        <- length(lines)
  is_part  <- grepl("^\\s*PART \\d", lines)
  is_rule  <- grepl("^={20,}$", lines)
  is_blank <- !nzchar(lines)
  is_colh  <- grepl("^\\s+Year Rel value", lines)
  hard     <- which(is_part | is_rule)
  starts   <- unique(c(1L, hard)); ends <- c(starts[-1L] - 1L, n)
  block_len <- ends - starts + 1L
  has_bar   <- vapply(seq_along(starts), function(i) any(is_colh[starts[i]:ends[i]]), logical(1))
  # a bar-plot block may hold several 400-year pages: the unbreakable unit is
  # one page = from one column header to the blank line before the next
  unit <- 0L
  for (i in which(has_bar)) {
    ch <- which(is_colh[starts[i]:ends[i]]) + starts[i] - 1L
    ch_end <- c(ch[-1L] - 1L, ends[i])
    unit <- max(unit, ch_end - ch + 1L)
  }
  # lines of header above the first column header of each block count too
  if (unit > 0L) unit <- unit + 4L

  # ---- font size: width first, then height of the unbreakable unit --------
  # Courier advance is 0.6 em: width of N chars at s points = N * 0.6 * s / 72 in
  fit_w <- 72 * usable_w / (0.6 * maxw)
  fit_h <- if (unit > 0L) 72 * usable_h / (1.15 * (unit + 2L)) else Inf
  fit_size <- min(fit_w, fit_h)
  if (is.null(font_size)) {
    font_size <- min(10, floor(fit_size * 2) / 2)       # at most 10 pt, in 0.5 steps
    if (font_size < 5) stop(sprintf(
      "Output needs %.1f pt to fit (width %d chars, longest block %d lines); split it or use landscape.",
      fit_size, maxw, unit))
  } else if (font_size > fit_w + 1e-6) {
    stop(sprintf("font_size = %g does not fit %d columns on this page (max %.1f pt).",
                 font_size, maxw, fit_w))
  }
  line_h  <- 1.15 * font_size / 72                       # leading
  hdr_h   <- 2 * line_h
  lpp     <- floor((usable_h - hdr_h) / line_h)          # lines per page

  # ---- page breaks -----------------------------------------------------------
  # Inside a bar-plot block, each 400-year page (column header .. next header)
  # is also a hard boundary.
  pages <- list(); start <- 1L
  while (start <= n) {
    end <- min(n, start + lpp - 1L)
    if (end < n) {
      cand <- seq(start + 1L, end + 1L)
      # prefer a hard boundary (part / series rule / next bar-plot page)
      # anywhere after the first line of the page
      hb <- cand[is_part[cand] | is_rule[cand] | is_colh[cand]]
      hb <- hb[hb > start]
      if (length(hb)) {
        # take the LAST boundary that still leaves the unit intact
        end <- max(hb) - 1L
        # a column header is preceded by its own part/rule/title lines: back up
        # to the start of that block so the header travels with its page
        if (is_colh[end + 1L]) {
          k <- end
          while (k > start && !is_blank[k] && !is_part[k] && !is_rule[k]) k <- k - 1L
          if (k > start) end <- k - 1L
        }
      }
      else {
        soft <- cand[is_blank[cand]]
        soft <- soft[soft > start + lpp %/% 2]
        if (length(soft)) end <- max(soft) - 1L
      }
    }
    pages[[length(pages) + 1L]] <- lines[start:end]
    start <- end + 1L
    while (start <= n && is_blank[start]) start <- start + 1L   # no leading blanks
  }

  # ---- draw --------------------------------------------------------------
  grDevices::pdf(file, width = dims[1L], height = dims[2L], family = family,
                 title = title, paper = "special")
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", family = family)
  on.exit(graphics::par(op), add = TRUE)
  cex <- font_size / 12          # pdf() base pointsize is 12
  for (p in seq_along(pages)) {
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, dims[1L]), ylim = c(0, dims[2L]))
    y0 <- dims[2L] - mar
    graphics::text(mar, y0 - 0.6 * line_h, title, adj = c(0, 0.5), cex = cex, font = 2)
    graphics::text(dims[1L] - mar, y0 - 0.6 * line_h,
                   sprintf("page %d / %d", p, length(pages)), adj = c(1, 0.5), cex = cex)
    graphics::segments(mar, y0 - hdr_h + 0.3 * line_h, dims[1L] - mar,
                       y0 - hdr_h + 0.3 * line_h, lwd = 0.5, col = "grey50")
    y <- y0 - hdr_h - seq_along(pages[[p]]) * line_h + 0.3 * line_h
    graphics::text(rep(mar, length(y)), y, pages[[p]], adj = c(0, 0), cex = cex)
  }
  grDevices::dev.off(); on.exit(graphics::par(op))   # close before opening the file

  if (isTRUE(open)) {
    if (.Platform$OS.type == "windows") shell.exec(file)
    else system2(if (Sys.info()[["sysname"]] == "Darwin") "open" else "xdg-open", shQuote(file))
  }
  invisible(list(file = file, pages = length(pages), font_size = font_size,
                 orientation = orientation, lines_per_page = lpp))
}
