
# =============================================================================
# dpl_cormat  --  correlation matrix of series in COFECHA style
# =============================================================================

#' Correlation matrix of series, printed in COFECHA style
#'
#' @description
#' Computes the pairwise Pearson correlations between all series of a
#' collection over their common years and prints them as a lower-triangular
#' matrix in the fixed-width style of the COFECHA output, with significance
#' marked by asterisks. Correlations can be computed on the raw measurements
#' or on the COFECHA-transformed series, in which case the year-to-year
#' variability drives the result exactly as in Part 5 of
#' \code{\link{dpl_cof}}.
#'
#' @details
#' With \code{type = "transformed"} (the default when \code{x} is a
#' \code{dpl_cof} result) the series are taken from
#' \code{cof_result$filtered}: spline-detrended with the run's
#' \code{spline_period}, log-transformed, autoregressively prewhitened and
#' normalised -- the values COFECHA correlates with the master. Raw
#' correlations (\code{type = "raw"}) are dominated by the shared growth
#' trend and age curve and are usually much higher; they say little about
#' dating but are useful to find duplicated samples (see
#' \code{Note}).
#'
#' Each pair uses its own overlap, so \code{n} differs from cell to cell.
#' Significance is the one-sided test of \eqn{r > 0} on \eqn{n - 2} degrees
#' of freedom: \code{*} \eqn{p < 0.05}, \code{**} \eqn{p < 0.01},
#' \code{***} \eqn{p < 0.001}; a pair with fewer than \code{min_overlap}
#' common years is printed as \code{-}. For the transformed series a
#' positive correlation is the expectation under correct dating; a negative
#' or near-zero value between two long series from the same site is a
#' dating warning.
#'
#' The summary rows give each series' mean correlation with all others
#' (\code{rbar} of the collection is their average) and the number of
#' significant pairs.
#'
#' @note Two columns that are the \emph{same} sample measured twice, or the
#' same core entered under two labels, show \eqn{r} close to 1 on the raw
#' series (typically \eqn{> 0.97}) as well as on the transformed ones; two
#' radii of the same tree show high but clearly lower values. A column of
#' raw \eqn{r} above 0.95 is worth checking before building a chronology, as
#' duplicates inflate every collection statistic.
#'
#' @param x A dplR \code{rwl} data.frame, or the list returned by
#'   \code{\link{dpl_cof}}.
#' @param type \code{"transformed"} or \code{"raw"}. Default: transformed
#'   when \code{x} is a \code{dpl_cof} result, raw when it is an \code{rwl}.
#'   Transformed from an \code{rwl} runs \code{\link{dpl_cof}} first with its
#'   defaults (pass \code{...} to change them).
#' @param series Series to include: character IDs or integer positions
#'   (columns of the \code{rwl}, or of \code{cof_result$filtered}). Default
#'   all.
#' @param min_overlap Integer (default \code{10L}). Pairs with fewer common
#'   years are shown as \code{-}.
#' @param digits Integer, 2 (default) or 3 decimals.
#' @param output_file Character or \code{NULL}. Also write the printed text to
#'   a file.
#' @param quiet Logical (default \code{FALSE}). Suppress console output.
#' @param ... Passed to \code{\link{dpl_cof}} when \code{x} is an \code{rwl}
#'   and \code{type = "transformed"}.
#'
#' @return Invisibly, a list with:
#' \describe{
#'   \item{\code{r}}{Symmetric correlation matrix (diagonal 1, \code{NA}
#'     below \code{min_overlap}).}
#'   \item{\code{n}}{Matrix of overlap lengths.}
#'   \item{\code{p}}{Matrix of one-sided p-values.}
#'   \item{\code{summary}}{\code{data.frame}: \code{series}, \code{n_pairs},
#'     \code{mean_r}, \code{n_sig} (pairs with \eqn{p < 0.05}).}
#'   \item{\code{rbar}}{Mean inter-series correlation over all valid pairs.}
#'   \item{\code{type}}{\code{"raw"} or \code{"transformed"}.}
#'   \item{\code{lines}}{The printed text.}
#' }
#'
#' @examples
#' \dontrun{
#' cof <- dpl_cof(rwl)
#' dpl_cormat(cof)                       # transformed, as in Part 5
#' dpl_cormat(rwl)                       # raw: look for duplicates near 1
#' cm <- dpl_cormat(cof, series = 1:8, digits = 3, quiet = TRUE)
#' cm$summary
#' dpl_print_pdf(cm$lines, "cormat.pdf")
#' }
#'
#' @seealso \code{\link{dpl_cof}}, \code{\link{dpl_print_pdf}}
#' @export
dpl_cormat <- function(x,
                       type        = NULL,
                       series      = NULL,
                       min_overlap = 10L,
                       digits      = 2L,
                       output_file = NULL,
                       quiet       = FALSE,
                       ...) {

  is_cof <- is.list(x) && !is.data.frame(x) && !is.null(x$master)
  if (is.null(type)) type <- if (is_cof) "transformed" else "raw"
  type <- match.arg(type, c("transformed", "raw"))

  if (type == "transformed") {
    if (!is_cof) x <- dpl_cof(x, verbose = FALSE, parts = integer(0), ...)
    if (is.null(x$filtered))
      stop("'x' has no $filtered element: re-run dpl_cof() with this version of cofechar.")
    dat <- x$filtered
    what <- sprintf("COFECHA-transformed series (%d-yr spline, log, AR, normalised)",
                    x$options$spline_period)
  } else {
    if (is_cof) stop("type = \"raw\" needs the rwl data.frame, not a dpl_cof() result.")
    if (!.is_rwl(x)) stop("'x' must be a dplR rwl data.frame or a dpl_cof() result.")
    dat  <- x
    what <- "raw measurements"
  }
  ids <- .resolve_series(series, colnames(dat))
  if (length(ids) < 2L) stop("At least two series are needed.")
  dat <- dat[, ids, drop = FALSE]
  k   <- length(ids)
  digits <- if (digits >= 3L) 3L else 2L

  # ---- correlations over pairwise overlap ---------------------------------
  R <- matrix(NA_real_, k, k, dimnames = list(ids, ids))
  N <- matrix(0L, k, k, dimnames = list(ids, ids))
  P <- matrix(NA_real_, k, k, dimnames = list(ids, ids))
  diag(R) <- 1; diag(N) <- colSums(!is.na(dat))
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k) {
    ok <- !is.na(dat[[i]]) & !is.na(dat[[j]])
    n  <- sum(ok); N[i, j] <- N[j, i] <- n
    if (n >= max(3L, min_overlap)) {
      r <- .cof_correl(dat[[i]][ok], dat[[j]][ok])
      R[i, j] <- R[j, i] <- r
      tt <- r * sqrt((n - 2L) / max(1e-12, 1 - r^2))
      P[i, j] <- P[j, i] <- stats::pt(tt, n - 2L, lower.tail = FALSE)
    }
  }
  stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""

  # ---- summary ---------------------------------------------------------------
  off <- R; diag(off) <- NA
  summ <- data.frame(
    series  = ids,
    n_pairs = rowSums(!is.na(off)),
    mean_r  = rowMeans(off, na.rm = TRUE),
    n_sig   = rowSums(!is.na(off) & P < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE)
  summ$mean_r[summ$n_pairs == 0L] <- NA_real_
  rbar <- mean(off[upper.tri(off)], na.rm = TRUE)

  # ---- ASCII -------------------------------------------------------------
  # Cell: F4.2 (' .51') or F5.3 (' .512') with leading zero suppressed, as in
  # Part 5, followed by up to 3 stars; cells are right-aligned in CW chars.
  CW  <- if (digits == 2L) 9L else 10L
  fmt <- function(r, p) {
    if (is.na(r)) return(formatC("-", width = CW))
    v <- sprintf(if (digits == 2L) "%.2f" else "%.3f", r)
    v <- sub("^0\\.", ".", v); v <- sub("^-0\\.", "-.", v)   # Fortran-style
    formatC(paste0(v, stars(p)), width = CW)
  }
  lab_w  <- max(8L, max(nchar(ids)))
  per_blk <- max(1L, (132L - lab_w - 2L) %/% CW)
  blocks  <- split(seq_len(k), ceiling(seq_len(k) / per_blk))

  out <- c(
    sprintf("CORRELATION MATRIX OF SERIES:  %s", what),
    strrep("-", 132L),
    sprintf(" Pearson r over the common years of each pair; '-' = fewer than %d common years.",
            as.integer(min_overlap)),
    " Significance (one-sided, r > 0):  * p < 0.05   ** p < 0.01   *** p < 0.001",
    "")
  for (b in seq_along(blocks)) {
    cols <- blocks[[b]]
    out <- c(out,
      paste0(formatC("", width = lab_w + 2L),
             paste(formatC(substr(ids[cols], 1L, CW - 1L), width = CW), collapse = "")),
      paste0(formatC("", width = lab_w + 2L),
             paste(rep(formatC(strrep("-", CW - 1L), width = CW), length(cols)), collapse = "")))
    for (i in seq_len(k)) {
      cells <- vapply(cols, function(j) {
        if (j > i) formatC("", width = CW)            # upper triangle blank
        else if (j == i) formatC("1", width = CW)
        else fmt(R[i, j], P[i, j])
      }, character(1))
      if (all(trimws(cells) == "")) next              # row entirely above diagonal
      out <- c(out, paste0(formatC(ids[i], width = lab_w, flag = "-"), "  ",
                           paste(cells, collapse = "")))
    }
    out <- c(out, "")
  }
  out <- c(out,
    sprintf(" %-*s  %7s  %7s  %7s", lab_w, "Series", "pairs", "mean r", "sig"),
    sprintf(" %-*s  %7s  %7s  %7s", lab_w, strrep("-", lab_w), "-----", "------", "---"))
  for (i in seq_len(k))
    out <- c(out, sprintf(" %-*s  %7d  %7s  %7d", lab_w, ids[i], summ$n_pairs[i],
                          if (is.na(summ$mean_r[i])) "-" else .cof_suppress_zero(sprintf("%6.3f", summ$mean_r[i]), 6L),
                          summ$n_sig[i]))
  out <- c(out, "",
    sprintf(" Mean inter-series correlation (rbar) over %d pairs: %s",
            sum(!is.na(off[upper.tri(off)])),
            .cof_suppress_zero(sprintf("%6.3f", rbar), 6L)))
  out <- sub("\\s+$", "", out)

  if (!quiet) cat(out, sep = "\n")
  if (!is.null(output_file)) writeLines(out, output_file)
  invisible(list(r = R, n = N, p = P, summary = summ, rbar = rbar,
                 type = type, lines = out))
}
