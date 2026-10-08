
# =============================================================================
# dpl_short  --  anchored dating check for very short series
# =============================================================================

#' Check very short series against a COFECHA master (anchored dating)
#'
#' @description
#' Tests series that are too short for COFECHA's sliding-segment procedure
#' (typically 5--20 years: microcores, xylogenesis samples, short anatomical
#' sequences) against the master dating series of a completed
#' \code{\link{dpl_cof}} run. Instead of sliding a segment along the whole
#' master, the series is assumed to be \emph{anchored} -- its last year is
#' known, as for living trees sampled on a known date -- and the question is
#' only whether the sequence is consistent with the master at lag 0, or
#' whether a missing or double ring near the bark has shifted it by one or
#' two years.
#'
#' @details
#' For each series the raw values are log-transformed and converted to
#' first differences (year-to-year changes), which removes the level and
#' the trend that a handful of years cannot estimate, and the same is done
#' to the master. At each lag in \code{-max_lag:max_lag} three quantities are
#' computed on the overlapping years:
#' \describe{
#'   \item{\code{r}}{Pearson correlation of the first differences.}
#'   \item{\code{p}}{One-sided exact p-value of \code{r} (t distribution with
#'     \code{n - 3} degrees of freedom, one fewer than usual because first
#'     differencing uses one year). \code{NA} when fewer than 5 differences.}
#'   \item{\code{glk}}{Gleichläufigkeit: the proportion of year-to-year
#'     changes with the same sign in series and master (Eckstein & Bauch
#'     1969), with its binomial p-value \code{glk_p}.}
#' }
#' \strong{Sign of the lag.} \code{best_lag} is the correction to apply to
#' the series' years: true year = labelled year + lag. A negative lag means
#' the rings are labelled too old (one ring too many was counted, e.g. a
#' false ring); a positive lag means they are labelled too young (a ring is
#' missing from the count).
#'
#' The verdict is \code{"ok"} when lag 0 has the highest \code{r} and
#' \code{p < alpha}; \code{"shifted"} when another lag is both higher and
#' significant (the lag is reported); \code{"weak"} when nothing is
#' significant, which for 5--8 years is the usual honest answer and should be
#' read together with \code{glk}.
#'
#' The master used is \code{cof_result$master_raw}, the mean of the normalised,
#' spline-filtered series. When the short series was itself part of the run,
#' pass the \code{cof_result} of a run with it excluded (the default
#' \code{min_length = 10} in \code{\link{dpl_cof}} does this for series under
#' 10 years), otherwise the test is circular.
#'
#' @param rwl A dplR \code{rwl} data.frame holding the short series (other
#'   series present are ignored unless named in \code{series}).
#' @param cof_result The list returned by \code{\link{dpl_cof}}.
#' @param series Character vector of series to test. Default: those listed in
#'   \code{cof_result$short}, or all columns of \code{rwl} if that is empty.
#' @param max_lag Integer (default \code{2L}). Shifts tested on each side of
#'   lag 0.
#' @param min_n Integer (default \code{5L}). Minimum overlapping years for a
#'   p-value to be reported.
#' @param alpha Numeric (default \code{0.05}). Significance level for the
#'   verdict.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{summary}}{One row per series: \code{series}, \code{jyr},
#'     \code{lyr}, \code{n}, \code{r0}, \code{p0}, \code{glk0}, \code{best_lag},
#'     \code{r_best}, \code{verdict}.}
#'   \item{\code{lags}}{One row per series and lag: \code{series}, \code{lag},
#'     \code{n}, \code{r}, \code{p}, \code{glk}, \code{glk_p}.}
#' }
#' Printed to the console: a compact table per series.
#'
#' @references
#' Eckstein, D. & Bauch, J. (1969). Beitrag zur Rationalisierung eines
#' dendrochronologischen Verfahrens und zur Analyse seiner Aussagesicherheit.
#' \emph{Forstwissenschaftliches Centralblatt} 88:230--250.
#'
#' @examples
#' \dontrun{
#' cores <- dpl_read_dec("MAI.rwl")
#' micro <- dpl_read_dec("MAI_xylo.rwl")          # 3- to 13-year series
#' all   <- dpl_merge(list(cores, micro))
#' cof   <- dpl_cof(all)                          # min_length = 10: short
#' cof$short                                      #   series set aside here
#' chk   <- dpl_short(all, cof)
#' chk$summary
#' subset(chk$lags, series == "6162")
#' }
#'
#' @seealso \code{\link{dpl_cof}}, \code{\link{dpl_dateme}}
#' @export
dpl_short <- function(rwl, cof_result,
                      series  = NULL,
                      max_lag = 2L,
                      min_n   = 5L,
                      alpha   = 0.05) {

  if (!.is_rwl(rwl)) stop("'rwl' must be a dplR rwl data.frame.")
  if (is.null(cof_result$master_raw))
    stop("'cof_result' must be the output of dpl_cof().")

  if (is.null(series)) {
    series <- if (!is.null(cof_result$short) && nrow(cof_result$short) > 0L)
      cof_result$short$series else colnames(rwl)
  }
  series <- intersect(series, colnames(rwl))
  if (length(series) == 0L) stop("No series to test.")

  max_lag <- as.integer(max_lag)
  years   <- as.integer(rownames(rwl))

  # Master: first differences of the log master, indexed by the year of the
  # change (value at year t is x[t] - x[t-1]).
  m_yrs <- as.integer(names(cof_result$master_raw))
  m_val <- as.numeric(cof_result$master_raw)
  dm    <- setNames(diff(m_val), m_yrs[-1L])

  lag_rows <- list(); sum_rows <- list()

  for (sid in series) {
    col <- rwl[[sid]]
    ok  <- which(!is.na(col))
    if (length(ok) < 3L) next
    s_yrs <- years[ok]
    s_val <- col[ok]
    # log-transform as COFECHA does, then first differences; a gap in the
    # series (non-consecutive years) breaks the difference at that point
    s_log <- log(pmax(s_val, 0.001))
    ds    <- setNames(diff(s_log), s_yrs[-1L])
    ds    <- ds[diff(s_yrs) == 1L]

    best <- NULL
    for (lag in seq(-max_lag, max_lag)) {
      # the series' change labelled year t is compared with the master's
      # change at year t + lag, i.e. true year = labelled year + lag
      yr_s <- as.integer(names(ds))
      yr_m <- yr_s + lag
      keep <- yr_m %in% as.integer(names(dm))
      n    <- sum(keep)
      r <- p <- glk <- glk_p <- NA_real_
      if (n >= 3L) {
        a <- ds[keep]; b <- dm[as.character(yr_m[keep])]
        r <- .cof_correl(a, b)
        if (n >= min_n) {
          df <- n - 2L
          tt <- r * sqrt(df / max(1e-12, 1 - r^2))
          p  <- stats::pt(tt, df, lower.tail = FALSE)
        }
        agree <- sum(sign(a) == sign(b))
        glk   <- agree / n
        glk_p <- stats::pbinom(agree - 1L, n, 0.5, lower.tail = FALSE)
      }
      lag_rows[[length(lag_rows) + 1L]] <- data.frame(
        series = sid, lag = lag, n = n, r = r, p = p,
        glk = glk, glk_p = glk_p, stringsAsFactors = FALSE)
    }
    tab  <- do.call(rbind, lag_rows[(length(lag_rows) - 2L * max_lag):length(lag_rows)])
    r0   <- tab$r[tab$lag == 0L]; p0 <- tab$p[tab$lag == 0L]; g0 <- tab$glk[tab$lag == 0L]
    ib   <- if (all(is.na(tab$r))) NA_integer_ else which.max(tab$r)
    blag <- if (is.na(ib)) NA_integer_ else tab$lag[ib]
    rb   <- if (is.na(ib)) NA_real_ else tab$r[ib]
    pb   <- if (is.na(ib)) NA_real_ else tab$p[ib]
    verdict <- if (is.na(ib) || is.na(pb)) "weak"
      else if (blag == 0L && pb < alpha) "ok"
      else if (blag != 0L && pb < alpha) "shifted"
      else "weak"
    sum_rows[[length(sum_rows) + 1L]] <- data.frame(
      series = sid, jyr = s_yrs[1L], lyr = s_yrs[length(s_yrs)], n = length(ok),
      r0 = r0, p0 = p0, glk0 = g0, best_lag = blag, r_best = rb,
      verdict = verdict, stringsAsFactors = FALSE)
  }

  summary <- do.call(rbind, sum_rows); lags <- do.call(rbind, lag_rows)
  rownames(summary) <- NULL; rownames(lags) <- NULL

  # ---- console report -------------------------------------------------------
  cat(sprintf("\n Anchored check of %d short series against the master (lags %+d..%+d)\n",
              nrow(summary), -max_lag, max_lag))
  cat(sprintf(" %-8s %4s-%-4s %3s  %6s %6s %5s   %4s %6s   %s\n",
              "Series", "from", "to", "n", "r(0)", "p(0)", "GLK", "best", "r", "verdict"))
  for (k in seq_len(nrow(summary))) {
    z <- summary[k, ]
    cat(sprintf(" %-8s %4d-%-4d %3d  %6.2f %6s %5s   %+4d %6.2f   %s\n",
                z$series, z$jyr, z$lyr, z$n,
                ifelse(is.na(z$r0), NA, z$r0),
                ifelse(is.na(z$p0), "   -", sprintf("%.3f", z$p0)),
                ifelse(is.na(z$glk0), "  -", sprintf("%.2f", z$glk0)),
                ifelse(is.na(z$best_lag), 0L, z$best_lag),
                ifelse(is.na(z$r_best), NA, z$r_best), z$verdict))
  }
  cat("\n")
  invisible(list(summary = summary, lags = lags))
}
