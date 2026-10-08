
# =============================================================================
# dpl_short_barplot  --  bar plot of short series side by side with the master
# =============================================================================

#' Bar plot of short series beside the COFECHA master
#'
#' @description
#' Prints the master dating series and one or more short series side by side,
#' one row per year, in the COFECHA bar-plot notation (\code{\link{dpl_barplot}}),
#' so that the year-to-year pattern of a microcore or other short sample can be
#' compared by eye with the site's pattern -- the visual counterpart of
#' \code{\link{dpl_short}}. Pointer years stand out as rows where every column
#' carries a short bar (narrow) or a long one (wide).
#'
#' @details
#' A series of a few years has no distribution of its own to rank against, so
#' unlike \code{\link{dpl_barplot}} the samples are \emph{not} normalised over
#' their own span. Instead every value is expressed on the master's scale:
#' the sample is log-transformed and first-differenced like the master, and
#' both are standardised with the master's mean and standard deviation over
#' the displayed window. A bar then reads \dQuote{this ring grew more (less)
#' than the previous one by so many site-SDs}, directly comparable across
#' columns. Because differences are used, the first year of each sample is
#' blank.
#'
#' Bar length encodes the decile of the master's own changes over the
#' window and the end symbol the standardised value (\code{A}--\code{Z}
#' positive, \code{a}--\code{z} negative, \code{@} near zero), exactly as in
#' Part 4 of the COFECHA output, so a reader of DPL bar plots needs nothing
#' new. A column \code{POOL} (the stand mean of the samples, as in
#' \code{dpl_short(pool = TRUE)}) is added when \code{pool = TRUE}.
#'
#' @param rwl A dplR \code{rwl} data.frame holding the short series.
#' @param cof_result The list returned by \code{\link{dpl_cof}}; its
#'   \code{$master_raw} is displayed and sets the scale.
#' @param series Character vector of series to show. Default: those in
#'   \code{cof_result$short}, else all columns of \code{rwl}.
#' @param years Integer vector \code{c(first, last)}. Display window. Default:
#'   the years covered by at least two of the selected series, extended by
#'   \code{margin} years on each side for context; one unusually long sample
#'   therefore does not stretch the plot.
#' @param margin Integer (default \code{3L}). Context years added on each side
#'   of the series span when \code{years} is \code{NULL}.
#' @param pool Logical (default \code{TRUE}). Add the stand-mean column.
#' @param lag Named integer vector, e.g. \code{c("6143" = -1)}: shift the named
#'   series by that many years before plotting, to see what a correction
#'   suggested by \code{\link{dpl_short}} would look like. Not applied to the
#'   data.
#' @param output_file Character or \code{NULL}. Also write the plot to a file.
#' @param quiet Logical (default \code{FALSE}). Suppress console output.
#'
#' @return Invisibly, the character vector of printed lines.
#'
#' @examples
#' \dontrun{
#' cof <- dpl_cof(rwl_all)                      # short series set aside
#' dpl_short(rwl_all, cof, pool = TRUE)         # the numbers
#' dpl_short_barplot(rwl_all, cof)              # the picture
#'
#' # what would 6143 look like shifted back one year?
#' dpl_short_barplot(rwl_all, cof, series = c("6143", "6144"), lag = c("6143" = -1))
#' }
#'
#' @seealso \code{\link{dpl_short}}, \code{\link{dpl_barplot}}
#' @export
dpl_short_barplot <- function(rwl, cof_result,
                           series      = NULL,
                           years       = NULL,
                           margin      = 3L,
                           pool        = TRUE,
                           lag         = NULL,
                           output_file = NULL,
                           quiet       = FALSE) {

  if (!.is_rwl(rwl)) stop("'rwl' must be a dplR rwl data.frame.")
  if (is.null(cof_result$master_raw))
    stop("'cof_result' must be the output of dpl_cof().")
  if (is.null(series)) {
    series <- if (!is.null(cof_result$short) && nrow(cof_result$short) > 0L)
      cof_result$short$series else colnames(rwl)
  }
  series <- intersect(series, colnames(rwl))
  if (length(series) == 0L) stop("No series to plot.")

  yrs_rwl <- as.integer(rownames(rwl))
  m_yrs   <- as.integer(names(cof_result$master_raw))
  dm      <- setNames(diff(as.numeric(cof_result$master_raw)), m_yrs[-1L])

  # ---- per-series log differences, optionally shifted ------------------------
  diffs <- list()
  for (sid in series) {
    ok <- which(!is.na(rwl[[sid]]))
    if (length(ok) < 2L) next
    y  <- yrs_rwl[ok]
    d  <- setNames(diff(log(pmax(rwl[[sid]][ok], 0.001))), y[-1L])
    d  <- d[diff(y) == 1L]
    if (!is.null(lag) && sid %in% names(lag))
      names(d) <- as.integer(names(d)) + as.integer(lag[[sid]])
    diffs[[sid]] <- d
  }
  if (length(diffs) == 0L) stop("No series with two or more consecutive years.")

  if (isTRUE(pool) && length(diffs) > 1L) {
    all_y <- sort(unique(unlist(lapply(diffs, function(d) as.integer(names(d))))))
    mat   <- sapply(diffs, function(d) d[as.character(all_y)])
    mat   <- matrix(mat, nrow = length(all_y))
    keep  <- rowSums(!is.na(mat)) >= 2L
    if (any(keep))
      diffs <- c(list(POOL = setNames(rowMeans(mat[keep, , drop = FALSE], na.rm = TRUE),
                                      all_y[keep])), diffs)
  }

  # ---- window --------------------------------------------------------------
  # Default: the years covered by at least two series (or by the single
  # series), so one unusually long sample does not stretch the plot.
  if (is.null(years)) {
    all_y <- unlist(lapply(diffs[names(diffs) != "POOL"], function(d) as.integer(names(d))))
    cnt   <- table(all_y)
    core  <- as.integer(names(cnt)[cnt >= min(2L, length(diffs[names(diffs) != "POOL"]))])
    if (length(core) == 0L) core <- as.integer(names(cnt))
    years <- c(min(core) - margin - 1L, max(core) + margin)
  }
  years <- as.integer(years)
  win   <- seq(years[1L], years[2L])
  win   <- win[win %in% as.integer(names(dm))]
  if (length(win) < 3L) stop("Display window overlaps the master in fewer than 3 years.")

  # ---- master scale over the window ----------------------------------------
  mw  <- dm[as.character(win)]
  nr  <- .cof_normts(as.numeric(mw), k = 0L)
  std <- function(v) if (nr$sd > 0) (v - nr$mean) / nr$sd else v - nr$mean
  Z   <- .cof_barpl_cuts(nr$z)

  cell <- function(yr, v) {
    if (is.na(v)) return(strrep(" ", 11L))
    substr(.cof_barpl_car(yr, std(v), Z), 6L, 16L)   # bar only, 11 chars
  }

  cols <- c("MASTER", names(diffs))
  getv <- function(cn, yr) {
    if (cn == "MASTER") unname(dm[as.character(yr)]) else {
      d <- diffs[[cn]]; v <- d[as.character(yr)]; if (length(v)) unname(v) else NA_real_
    }
  }

  # ---- assemble -----------------------------------------------------------
  BAR_W <- 12L
  hdr   <- paste0("      ", paste(formatC(cols, width = BAR_W - 1L, flag = "-"), collapse = " "))
  span  <- vapply(cols, function(cn) {
    if (cn == "MASTER") sprintf("%d-%d", min(m_yrs), max(m_yrs))
    else { y <- as.integer(names(diffs[[cn]])); sprintf("%d-%d", min(y) - 1L, max(y)) }
  }, character(1))
  span_line <- paste0("      ", paste(formatC(span, width = BAR_W - 1L, flag = "-"), collapse = " "))
  rule <- strrep("-", min(132L, 6L + length(cols) * BAR_W))

  out <- c(
    strrep("=", min(132L, 6L + length(cols) * BAR_W)),
    sprintf(" Short series vs master, %d-%d: log first differences on the master's scale (window mean/SD)",
            min(win), max(win)),
    rule, hdr, span_line, rule)
  for (i in seq_along(win)) {
    yr <- win[i]
    cells <- vapply(cols, function(cn) {
      formatC(sub("\\s+$", "", cell(yr, getv(cn, yr))), width = BAR_W - 1L, flag = "-")
    }, character(1))
    out <- c(out, paste0(sprintf("%4d  ", yr), paste(cells, collapse = " ")))
    if (i %% 10L == 0L && i < length(win)) out <- c(out, "")
  }
  out <- c(out, rule,
           " Bar: decile of the master's year-to-year changes over the window; symbol: change in",
           " master SDs (A-Z up, a-z down, @ ~0). First year of each sample is blank (differences).")
  out <- sub("\\s+$", "", out)

  if (!quiet) cat(out, sep = "\n")
  if (!is.null(output_file)) writeLines(out, output_file)
  invisible(out)
}
