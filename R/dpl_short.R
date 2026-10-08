
# =============================================================================
# dpl_short  --  anchored dating check for very short series
# =============================================================================

#' Check very short series against a COFECHA master (anchored dating)
#'
#' @description
#' Tests series that are too short for COFECHA's sliding-segment procedure
#' -- microcores, xylogenesis samples, short wood-anatomy or isotope
#' sequences of roughly 4 to 20 years -- against the master dating series of
#' a completed \code{\link{dpl_cof}} run. Instead of sliding a segment along
#' the whole master, the series is assumed to be \emph{anchored}: its last
#' year is known, as for living trees sampled on a known date, so the only
#' question is whether the sequence agrees with the master at lag 0 or
#' whether a missing or false ring near the bark has shifted it by one or two
#' years. With \code{pool = TRUE} the series are also combined into one
#' stand series, which can be tested even when the individual ones are too
#' short.
#'
#' @details
#' \subsection{Why COFECHA cannot test short series}{
#' COFECHA correlates segments of 50 years (20 at the shortest in common
#' use) with the master at 21 lags and compares the result with a critical
#' value tabulated from 10 years upwards. For a series of 3--10 years the
#' segment collapses to the series itself, the spline, autoregressive model
#' and log-transform have nothing to estimate, and the correlation at any
#' lag is dominated by chance: on 3 or 4 values \eqn{r} is routinely above
#' 0.95 whether or not the dating is right. For this reason
#' \code{\link{dpl_cof}} sets such series aside (\code{min_length}) and
#' \code{dpl_short} tests them differently.
#' }
#'
#' \subsection{What is tested}{
#' Each series is log-transformed, as in COFECHA, and converted to first
#' differences (year-to-year changes). First differences remove the level
#' and the growth trend, which a handful of years cannot estimate, and keep
#' what dating is about: whether the series went up or down in the same
#' years as the master. The master (\code{cof_result$master_raw}) is
#' differenced in the same way. At each lag in \code{-max_lag:max_lag} the
#' overlapping differences give:
#' \describe{
#'   \item{\code{r}}{Pearson correlation of the differences.}
#'   \item{\code{p}}{One-sided exact p-value of \code{r} (Student t with
#'     \code{n - 2} degrees of freedom, \code{n} = number of differences).
#'     Reported only when \code{n >= min_n}; below that no p-value is
#'     defensible.}
#'   \item{\code{glk}}{Gleichl\enc{ä}{ae}ufigkeit (Eckstein & Bauch 1969): the
#'     proportion of year-to-year changes with the same sign in series and
#'     master, with its exact binomial p-value \code{glk_p} (see the GLK
#'     section below).}
#' }
#' }
#'
#' \subsection{Gleichl\enc{ä}{ae}ufigkeit (GLK)}{
#' GLK asks one question of every pair of consecutive years: did the series
#' and the master change in the same direction? It is the fraction of years
#' where they did, ignoring by how much. This is what a dendrochronologist
#' does by eye on a skeleton plot, and it is the dating statistic of the
#' German school (TSAP, CDendro).
#'
#' Its significance is exact at any length: if the series were unrelated to
#' the master, each year's agreement would be a coin toss, so with \eqn{n}
#' differences and \eqn{k} agreements the one-sided p-value (\code{glk_p})
#' is the binomial tail \eqn{P(X \ge k \mid n, 1/2)}. Perfect agreement
#' gives \eqn{p = 0.5^n}: 0.125 for 3 differences, 0.06 for 4, 0.03 for 5,
#' 0.008 for 7, 0.001 for 10. One disagreement roughly quadruples these.
#'
#' For 4--6 years GLK is the statistic to trust: its p-value needs no
#' distributional assumption, whereas the t-test behind \code{p} for
#' \code{r} does, and it cannot be dragged by one extreme ring (a partly
#' sampled first or last ring is common in microcores). Its weakness is the
#' reverse: it discards the magnitudes, so when the signal is real it has
#' less power than \code{r}, small changes count as much as pointer years,
#' and on long series it rarely exceeds 0.75 even when dating is certain
#' (0.60--0.65 is the usual working threshold at \eqn{n \ge 50}). Read
#' both: when \code{r} and GLK agree the dating is as secure as a short
#' series allows; when they diverge, the better-looking one is probably
#' being carried by a single large ring.
#' }
#'
#' \subsection{Reading the verdict}{
#' \describe{
#'   \item{\code{"ok"}}{Lag 0 has the highest \code{r} and \code{p < alpha}:
#'     the series is consistent with its labelled dates.}
#'   \item{\code{"shifted"}}{Another lag has the highest \code{r} with
#'     \code{p < alpha}: the series is probably misdated by \code{best_lag}
#'     years. Check the sample before correcting.}
#'   \item{\code{"weak"}}{No lag is significant. For 4--6 years this is the
#'     usual and honest answer; read \code{r0} and \code{glk0} as
#'     \emph{indications} and look at the sample.}
#'   \item{\code{"untestable"}}{Fewer than 3 differences (series of 3 years
#'     or less, or broken by gaps).}
#' }
#' \strong{Sign of the lag.} \code{best_lag} is the correction to apply to
#' the series' years: true year = labelled year + lag. A negative lag means
#' the rings are labelled too old (a ring too many was counted, e.g. a false
#' ring); a positive lag means they are labelled too young (a ring is
#' missing from the count).
#' }
#'
#' \subsection{How many years are needed}{
#' With \eqn{n} measured years there are \eqn{n - 1} differences and, at a
#' non-zero lag, one fewer still. The best any series can achieve at lag 0 is:
#' \tabular{lll}{
#'   years \tab differences \tab what can be said \cr
#'   3     \tab 2 \tab nothing (\code{"untestable"}) \cr
#'   4     \tab 3 \tab \code{r} and GLK only, no p-value \cr
#'   5     \tab 4 \tab \code{r} and GLK only (default \code{min_n = 5}) \cr
#'   6     \tab 5 \tab p down to 0.03 if agreement is perfect; GLK 5/5, p = 0.03 \cr
#'   8     \tab 7 \tab p < 0.01 achievable; GLK 7/7, p = 0.008 \cr
#'   10+   \tab 9+ \tab a normal, if short, test \cr
#' }
#' In practice a single series needs about 6 years for an \code{"ok"} and 8
#' or more for a confident one. Below that, pooling is the way forward.
#' }
#'
#' \subsection{Pooling (\code{pool = TRUE})}{
#' Several short series from the same stand and the same years share the
#' pointer years, so their mean carries far more dating evidence than any
#' one of them. With \code{pool = TRUE} the log-differences of all tested
#' series are averaged year by year (years covered by at least
#' \code{pool_min} series) into one series named \code{"POOL"}, which is
#' tested like the others. Eleven 3- to 8-year microcores from 1992--2004,
#' for example, give one pooled series of up to 13 differences. The pooled
#' verdict answers \dQuote{is this \emph{set} of samples correctly dated};
#' it cannot say which individual is off, but if the pool is \code{"ok"}
#' while one series is \code{"shifted"} or has a low GLK, that series is
#' the one to examine. Pooling assumes the series were dated independently;
#' if they were all counted from the same photograph by the same person, a
#' shared error would be pooled too.
#' }
#'
#' \subsection{Workflow}{
#' \enumerate{
#'   \item Read cores and short series, merge them, run
#'     \code{dpl_cof(all)}. With the default \code{min_length = 10} the short
#'     series are set aside (\code{cof$short}) and the master is built from
#'     the cores alone -- essential, otherwise the test is circular.
#'   \item \code{dpl_short(all, cof, pool = TRUE)}: read the pooled verdict
#'     first, then the individuals.
#'   \item Examine any \code{"shifted"} series and any with low GLK at the
#'     microscope; correct with \code{\link{dpl_edt}} (insert/delete a ring)
#'     and re-run.
#'   \item Series that remain \code{"weak"} or \code{"untestable"} are not
#'     wrong -- they are too short to be checked numerically. Keep them if
#'     their dates rest on other evidence (bark year, pointer years seen on
#'     the slide, the pooled result); say so in the metadata.
#' }
#' }
#'
#' \subsection{Limits}{
#' The test assumes the last year is anchored; it is not a dating search.
#' For a floating short series use \code{\link{dpl_dateme}}, bearing in mind
#' that a 10-year fragment slid along a long master will produce spurious
#' matches. It also assumes that the short series responds to the same
#' signal as the master: a microcore variable other than ring width (cell
#' number, lumen area) may follow the ring-width master only loosely, in
#' which case a low GLK indicates a different signal rather than a dating
#' error.
#' }
#'
#' @param rwl A dplR \code{rwl} data.frame holding the short series (other
#'   series present are ignored unless named in \code{series}).
#' @param cof_result The list returned by \code{\link{dpl_cof}}. Its master
#'   must have been built without the series being tested (see Workflow).
#' @param series Character vector of series to test. Default: those listed in
#'   \code{cof_result$short}, or all columns of \code{rwl} if that is empty.
#' @param max_lag Integer (default \code{2L}). Shifts tested on each side of
#'   lag 0.
#' @param min_n Integer (default \code{5L}). Minimum number of overlapping
#'   differences for a p-value to be reported.
#' @param alpha Numeric (default \code{0.05}). Significance level for the
#'   verdict.
#' @param pool Logical (default \code{FALSE}). Also test the mean
#'   log-difference series of all tested series, reported as \code{"POOL"}.
#' @param pool_min Integer (default \code{2L}). Minimum number of series
#'   contributing to a year for it to enter the pooled series.
#'
#' @return Invisibly, a list with:
#' \describe{
#'   \item{\code{summary}}{One row per series (plus \code{"POOL"} if
#'     requested): \code{series}, \code{jyr}, \code{lyr}, \code{n} (measured
#'     years; for the pool, years covered), \code{r0}, \code{p0},
#'     \code{glk0}, \code{glk_p0}, \code{best_lag}, \code{r_best},
#'     \code{verdict}.}
#'   \item{\code{lags}}{One row per series and lag: \code{series}, \code{lag},
#'     \code{n} (differences), \code{r}, \code{p}, \code{glk}, \code{glk_p}.}
#'   \item{\code{pool}}{When \code{pool = TRUE}: the pooled log-difference
#'     series (named numeric vector, names = years) and the number of series
#'     contributing to each year (\code{depth}).}
#' }
#' A compact table is printed to the console.
#'
#' @references
#' Eckstein, D. & Bauch, J. (1969). Beitrag zur Rationalisierung eines
#' dendrochronologischen Verfahrens und zur Analyse seiner Aussagesicherheit.
#' \emph{Forstwissenschaftliches Centralblatt} 88:230--250.
#'
#' @examples
#' \dontrun{
#' cores <- dpl_read_dec("MAI.rwl")
#' micro <- dpl_read_dec("MAI_xylo.rwl")          # 3- to 13-year microcores
#' all   <- dpl_merge(list(cores, micro))
#' cof   <- dpl_cof(all)                          # min_length = 10 sets the
#' cof$short                                      #   short ones aside
#'
#' chk <- dpl_short(all, cof, pool = TRUE)
#' chk$summary                                    # POOL row first, then each
#' subset(chk$lags, series == "POOL")             # pooled lag table
#'
#' # a series flagged "shifted" by -1: one ring too many near the bark;
#' # inspect the sample, then e.g. delete the suspect ring and re-check
#' micro2 <- dpl_edt(micro, edits = list(
#'   list(series = "6143", op = "delete", year = 2003, move = "back")))
#' dpl_short(dpl_merge(list(cores, micro2)), cof, series = "6143")
#' }
#'
#' @seealso \code{\link{dpl_cof}}, \code{\link{dpl_dateme}},
#'   \code{\link{dpl_edt}}
#' @export
dpl_short <- function(rwl, cof_result,
                      series   = NULL,
                      max_lag  = 2L,
                      min_n    = 5L,
                      alpha    = 0.05,
                      pool     = FALSE,
                      pool_min = 2L) {

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
  dm    <- setNames(diff(as.numeric(cof_result$master_raw)), m_yrs[-1L])

  # ---- per-series log first differences -------------------------------------
  # A gap (non-consecutive years) breaks the difference at that point.
  diffs <- list(); spans <- list()
  for (sid in series) {
    col <- rwl[[sid]]
    ok  <- which(!is.na(col))
    if (length(ok) == 0L) next
    s_yrs <- years[ok]
    s_log <- log(pmax(col[ok], 0.001))
    ds    <- setNames(diff(s_log), s_yrs[-1L])
    diffs[[sid]] <- ds[diff(s_yrs) == 1L]
    spans[[sid]] <- c(s_yrs[1L], s_yrs[length(s_yrs)], length(ok))
  }

  # ---- pooled stand series ---------------------------------------------------
  pool_out <- NULL
  if (isTRUE(pool) && length(diffs) >= 1L) {
    all_yrs <- sort(unique(unlist(lapply(diffs, function(d) as.integer(names(d))))))
    mat <- sapply(diffs, function(d) d[as.character(all_yrs)])
    mat <- matrix(mat, nrow = length(all_yrs))
    depth <- rowSums(!is.na(mat))
    keep  <- depth >= pool_min
    if (any(keep)) {
      pooled <- setNames(rowMeans(mat[keep, , drop = FALSE], na.rm = TRUE),
                         all_yrs[keep])
      pool_out <- list(series = pooled, depth = setNames(depth[keep], all_yrs[keep]))
      diffs <- c(list(POOL = pooled), diffs)
      spans <- c(list(POOL = c(min(all_yrs[keep]) - 1L, max(all_yrs[keep]),
                               sum(keep) + 1L)), spans)
    }
  }

  # ---- test each series at each lag -----------------------------------------
  lag_rows <- list(); sum_rows <- list()
  for (sid in names(diffs)) {
    ds <- diffs[[sid]]
    tab <- vector("list", 2L * max_lag + 1L); k <- 0L
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
      k <- k + 1L
      tab[[k]] <- data.frame(series = sid, lag = lag, n = n, r = r, p = p,
                             glk = glk, glk_p = glk_p, stringsAsFactors = FALSE)
    }
    tab <- do.call(rbind, tab)
    lag_rows[[sid]] <- tab

    r0 <- tab$r[tab$lag == 0L]; p0 <- tab$p[tab$lag == 0L]
    g0 <- tab$glk[tab$lag == 0L]; gp0 <- tab$glk_p[tab$lag == 0L]
    ib <- if (all(is.na(tab$r))) NA_integer_ else which.max(tab$r)
    blag <- if (is.na(ib)) NA_integer_ else tab$lag[ib]
    rb   <- if (is.na(ib)) NA_real_ else tab$r[ib]
    pb   <- if (is.na(ib)) NA_real_ else tab$p[ib]
    verdict <- if (is.na(ib)) "untestable"
      else if (is.na(pb)) "weak"
      else if (blag == 0L && pb < alpha) "ok"
      else if (blag != 0L && pb < alpha) "shifted"
      else "weak"
    sp <- spans[[sid]]
    sum_rows[[sid]] <- data.frame(
      series = sid, jyr = sp[1L], lyr = sp[2L], n = sp[3L],
      r0 = r0, p0 = p0, glk0 = g0, glk_p0 = gp0, best_lag = blag, r_best = rb,
      verdict = verdict, stringsAsFactors = FALSE)
  }

  summary <- do.call(rbind, sum_rows); lags <- do.call(rbind, lag_rows)
  rownames(summary) <- NULL; rownames(lags) <- NULL

  # ---- console report -------------------------------------------------------
  cat(sprintf("\n Anchored check of %d short series against the master (lags %+d..%+d)\n",
              sum(summary$series != "POOL"), -max_lag, max_lag))
  cat(sprintf(" %-8s %4s-%-4s %3s  %6s %6s %5s %6s   %4s %6s   %s\n",
              "Series", "from", "to", "n", "r(0)", "p(0)", "GLK", "p", "best", "r", "verdict"))
  fmt <- function(v, f, na = "     -") if (is.na(v)) na else sprintf(f, v)
  for (k in seq_len(nrow(summary))) {
    z <- summary[k, ]
    cat(sprintf(" %-8s %4d-%-4d %3d  %6s %6s %5s %6s   %4s %6s   %s\n",
                z$series, z$jyr, z$lyr, z$n,
                fmt(z$r0, "%6.2f"), fmt(z$p0, "%6.3f"), fmt(z$glk0, "%5.2f", "    -"),
                fmt(z$glk_p0, "%6.3f"),
                fmt(z$best_lag, "%+4d", "   -"), fmt(z$r_best, "%6.2f"), z$verdict))
    if (z$series == "POOL" && nrow(summary) > 1L)
      cat(sprintf(" %s\n", strrep("-", 79L)))
  }
  if (!is.null(pool_out))
    cat(sprintf("\n POOL: mean log-difference of %d series, %d years with >= %d series\n",
                length(diffs) - 1L, length(pool_out$series), pool_min))
  cat("\n")
  invisible(list(summary = summary, lags = lags, pool = pool_out))
}
