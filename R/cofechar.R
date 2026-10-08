# =============================================================================
# xDPL-CATES — R port of the DPL (Dendrochronology Program Library)
# Original Fortran code (C) 1986–1994 Richard L. Holmes
# R translation: CATES ERC / CNRS-LSCE, 2026
#
# Modules: EDT — Edit ring-measurement series   (DPLEDT, APR 1988 / FEB 1993)
#          COF — COFECHA quality control         (DPLCOF, JAN 1982 / MAY 1994)
#
# See ?cofechar for full package documentation and workflow examples.
# =============================================================================

# =============================================================================
# Internal helpers
# =============================================================================
#
# These functions are not exported and should not be called directly.
# They translate the core Fortran utility routines (LARGO, DATR, DATW, TRRW,
# TRDISP) and provide the dplR conversion layer.


# -----------------------------------------------------------------------------
# .largo  —  trailing-whitespace trimmer  (LARGO equivalent)
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Returns the effective length of a string up to the last non-space
#   character, mirroring the Fortran LARGO subroutine used throughout DPL
#   for right-trimming fixed-length character fields before output.
#
# ARGUMENTS
#   x   Character scalar.
#
# RETURNS
#   Integer: number of characters up to and including the last non-space.
#
.largo <- function(x) {
  nchar(trimws(x, which = "right"))
}


# -----------------------------------------------------------------------------
# .rwl_to_series  —  dplR rwl to xDPL internal series list
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Converts a dplR-format data.frame into the xDPL internal list of series
#   so that the editing engine (dpl_edt) and file writers can operate on a
#   single canonical structure regardless of where the data came from.
#
#   Leading and trailing NAs in each column are stripped: yr_start is set
#   to the first non-NA year, and values spans only up to the last non-NA
#   year.  Internal NAs (i.e. absent rings within a series' span) are
#   preserved exactly as-is.
#
# ARGUMENTS
#   rwl   A dplR-style data.frame.  Row names must be integer-coercible
#         calendar years forming a consecutive sequence.  Each column is one
#         series; cells outside the series' span are NA.
#
# RETURNS
#   A list of length ncol(rwl).  Each element is a named list:
#     id        Character. Series identifier (taken from the column name).
#     yr_start  Integer. First non-NA calendar year.
#     values    Numeric vector spanning yr_start to the last non-NA year
#               (internal NAs retained).
#   Empty columns (all NA) produce yr_start = NA_integer_ and values = numeric(0).
#
.rwl_to_series <- function(rwl) {
  years       <- as.integer(rownames(rwl))
  series_list <- vector("list", ncol(rwl))

  for (j in seq_len(ncol(rwl))) {
    vals <- rwl[[j]]
    ok   <- which(!is.na(vals))
    if (length(ok) == 0L) {
      series_list[[j]] <- list(id       = colnames(rwl)[j],
                               yr_start = NA_integer_,
                               values   = numeric(0))
      next
    }
    series_list[[j]] <- list(
      id       = colnames(rwl)[j],
      yr_start = years[ok[1L]],
      values   = vals[ok[1L]:ok[length(ok)]]
    )
  }
  series_list
}


# -----------------------------------------------------------------------------
# .series_to_rwl  —  xDPL internal series list to dplR rwl data.frame
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Converts the xDPL internal series list back into a dplR-format data.frame
#   after editing.  The output spans the full calendar-year range of all
#   surviving series; cells outside each series' span are filled with NA.
#   The "rwl" class is explicitly assigned so that all dplR S3 generics
#   (time(), rwl.stats(), treering.plot(), detrend(), chron(), etc.) dispatch
#   correctly without any further conversion by the user.
#
# ARGUMENTS
#   series_list   A list of series in xDPL internal format (list with id,
#                 yr_start, values).  Empty series (length(values) == 0) are
#                 silently dropped before building the output frame.
#
# RETURNS
#   An object of class c("rwl", "data.frame"):
#     row names   Character strings of consecutive integer years.
#     columns     One column per series, named by series$id.
#     values      NA outside each series' span; ring widths within.
#   Returns an empty rwl object if all series are empty.
#
# NOTE
#   yr_start is coerced to integer explicitly to guard against numeric/integer
#   type ambiguity that can arise when edits are applied inside dpl_edt().
#   Character row names are required because all dplR generics recover the
#   year axis via as.numeric(rownames(x)).
#
.series_to_rwl <- function(series_list) {
  series_list <- Filter(function(s) length(s$values) > 0L, series_list)
  if (length(series_list) == 0L)
    return(structure(data.frame(), class = c("rwl", "data.frame")))

  yr_starts <- vapply(series_list, function(s) as.integer(s$yr_start), integer(1))
  yr_ends   <- vapply(series_list, function(s)
                 as.integer(s$yr_start) + length(s$values) - 1L, integer(1))

  yr_min    <- min(yr_starts)
  yr_max    <- max(yr_ends)
  all_years <- seq(yr_min, yr_max)
  n_yr      <- length(all_years)

  ids <- vapply(series_list, `[[`, character(1), "id")

  out <- as.data.frame(
    matrix(NA_real_, nrow = n_yr, ncol = length(series_list),
           dimnames = list(as.character(all_years), ids))
  )

  for (s in series_list) {
    if (length(s$values) == 0L) next
    iyr  <- as.integer(s$yr_start)
    lyr  <- iyr + length(s$values) - 1L
    rows <- match(seq(iyr, lyr), all_years)
    out[rows, s$id] <- s$values
  }

  class(out) <- c("rwl", "data.frame")
  out
}


# -----------------------------------------------------------------------------
# .is_rwl  —  detect whether an object is a dplR rwl data.frame
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Lightweight type check used throughout the public API to dispatch between
#   dplR data.frame input and xDPL internal series list input.
#
# ARGUMENTS
#   x   Any R object.
#
# RETURNS
#   TRUE if x is a data.frame whose row names are all coercible to integer
#   (i.e. look like calendar years); FALSE otherwise.
#
.is_rwl <- function(x) {
  is.data.frame(x) &&
    !is.null(rownames(x)) &&
    !anyNA(suppressWarnings(as.integer(rownames(x))))
}


# =============================================================================
# File I/O internals  (translations of DPL Fortran I/O routines)
# =============================================================================


# -----------------------------------------------------------------------------
# .read_compact  —  read a DPL compact ('~') file  [DATR equivalent]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Parses the DPL compact binary-integer format in which each series occupies
#   two logical records: a fixed-column header line ending in '~' and one or
#   more data lines of packed integers.  This is Holmes' native internal
#   format, produced by the DPL DATW routine and consumed by most DPL modules.
#
# COMPACT FILE LAYOUT
#   Header line (80 columns):
#     cols  1- 8   N         number of values (integer)
#     cols 11-18   JYR       first calendar year (integer)
#     cols 22-61   series ID (up to 40 characters)
#     cols 69-71   NC        decimal exponent  (stored_int x 10^NC = value)
#     cols 72-79   format    Fortran format string, e.g. "(16F5.0)"
#     col  80      '~'       header sentinel
#   Data lines: packed fixed-width integers, width derived from format string.
#
# ARGUMENTS
#   con   Character file path or open text connection.
#
# RETURNS
#   A list of series; each element is list(id, yr_start, values).
#
.read_compact <- function(con) {
  lines <- readLines(con, warn = FALSE)

  series_list <- list()
  i <- 1L

  while (i <= length(lines)) {
    ln <- lines[[i]]

    if (nchar(ln) == 0 || substr(ln, nchar(ln), nchar(ln)) != "~") {
      i <- i + 1L
      next
    }

    n_vals   <- as.integer(trimws(substr(ln,  1,  8)))
    yr_start <- as.integer(trimws(substr(ln, 11, 18)))
    id       <- trimws(substr(ln, 22, 61))
    nc       <- as.integer(trimws(substr(ln, 69, 71)))
    fmt_str  <- trimws(substr(ln, 72, 79))

    fld_w <- as.integer(regmatches(fmt_str,
               regexpr("(?<=[FIfI])\\d+", fmt_str, perl = TRUE)))
    if (length(fld_w) == 0 || is.na(fld_w)) fld_w <- 5L

    raw_vals     <- integer(0)
    i            <- i + 1L
    chars_needed <- n_vals * fld_w
    chars_read   <- 0L

    while (i <= length(lines) && chars_read < chars_needed) {
      dl  <- lines[[i]]
      pos <- 1L
      while (pos + fld_w - 1L <= nchar(dl)) {
        raw_vals <- c(raw_vals,
                      as.integer(trimws(substr(dl, pos, pos + fld_w - 1L))))
        pos <- pos + fld_w
      }
      chars_read <- chars_read + nchar(dl)
      i <- i + 1L
    }

    fac  <- 10^nc
    vals <- raw_vals[seq_len(n_vals)] * fac

    series_list[[length(series_list) + 1L]] <- list(
      id       = id,
      yr_start = yr_start,
      values   = vals
    )
  }

  series_list
}


# -----------------------------------------------------------------------------
# .read_tucson  —  read a Tucson .rwl file  [TRRW 'MR' branch equivalent]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Parses the ITRDB Tucson measurement format, the de-facto standard for
#   ring-width exchange.  Each decade is printed on one line; values are
#   stored as F6.2 (0.01 mm precision).  Series end with a sentinel value
#   of 9.99 (or -9999 for 0.001 mm precision data, which is auto-divided
#   by 10 on read).
#
# TUCSON LINE LAYOUT
#   cols  1-8    series ID (8 chars; or 7 chars if col 8 is '-' for long IDs)
#   cols  9-12   decade year (or cols 9-13 with the extended-ID variant)
#   remaining    up to 10 values of 6 characters each (F6.2)
#   sentinel     9.99 (0.01 mm) or -9999 / -99.99 (0.001 mm, divided by 10)
#
# ARGUMENTS
#   path   Character file path.
#
# RETURNS
#   A list of series; each element is list(id, yr_start, values).
#
.read_tucson <- function(path) {
  lines       <- readLines(path, warn = FALSE)
  series_list <- list()
  cur_id      <- NULL
  cur_yr      <- NULL
  cur_vals    <- numeric(0)

  for (ln in lines) {
    if (nchar(ln) < 12) next

    if (nchar(ln) >= 8 && substr(ln, 8, 8) == "-") {
      id        <- trimws(substr(ln, 1, 7))
      yr_start  <- as.integer(trimws(substr(ln, 9, 13)))
      val_start <- 14L
    } else {
      id        <- trimws(substr(ln, 1, 8))
      yr_start  <- as.integer(trimws(substr(ln, 9, 12)))
      val_start <- 13L
    }

    vals <- numeric(0)
    for (k in seq_len(10)) {
      s <- val_start + (k - 1L) * 6L
      e <- s + 5L
      if (e > nchar(ln)) break
      v <- suppressWarnings(as.numeric(trimws(substr(ln, s, e))))
      if (is.na(v)) break
      if (abs(v - 9.99) < 1e-4 || abs(v + 99.99) < 1e-4) {
        if (v < -99.98 && length(cur_vals) > 0) cur_vals <- cur_vals / 10
        if (!is.null(cur_id) && length(cur_vals) > 0) {
          series_list[[length(series_list) + 1L]] <- list(
            id       = cur_id,
            yr_start = cur_yr,
            values   = cur_vals
          )
        }
        cur_id   <- NULL
        cur_yr   <- NULL
        cur_vals <- numeric(0)
        break
      }
      vals <- c(vals, v)
    }

    if (length(vals) > 0) {
      if (is.null(cur_id)) {
        cur_id <- id
        cur_yr <- yr_start
      }
      cur_vals <- c(cur_vals, vals)
    }
  }

  if (!is.null(cur_id) && length(cur_vals) > 0) {
    series_list[[length(series_list) + 1L]] <- list(
      id       = cur_id,
      yr_start = cur_yr,
      values   = cur_vals
    )
  }

  series_list
}


# -----------------------------------------------------------------------------
# .write_compact  —  write a DPL compact ('~') file  [DATW equivalent]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Writes the xDPL internal series list to a DPL compact file, replicating
#   Holmes' DATW algorithm for automatic field-width and scaling-exponent
#   selection.  The exponent NC and field width NE are chosen to pack values
#   as compactly as possible without loss of precision.
#
# ARGUMENTS
#   series_list   List of series in xDPL internal format.
#   path          Output file path (overwritten if it exists).
#
.write_compact <- function(series_list, path) {
  con <- file(path, "wt")
  on.exit(close(con))

  for (s in series_list) {
    id  <- s$id
    n   <- length(s$values)
    jyr <- s$yr_start
    y   <- s$values
    if (n == 0) next

    fac <- max(1, max(abs(y), na.rm = TRUE), max(-y * 10, na.rm = TRUE))

    fac1 <- 1e16
    for (k in seq_len(32)) {
      fac1 <- fac1 * 0.1
      if (fac * fac1 < 999999.5) break
    }

    izz <- round(y * fac1)

    m  <- 1L
    ne <- 1L
    for (i in seq_len(6)) {
      m <- m * 10L
      if (any(izz %% m != 0)) {
        m  <- m %/% 10L
        ne <- 7L - i
        break
      }
      ne <- 1L
    }
    if (m >= 10) izz <- izz %/% m

    vmax <- max(abs(y), na.rm = TRUE)
    lmax <- which.max(abs(y))
    v    <- if (izz[lmax] != 0) abs(vmax / izz[lmax]) else 0
    nc   <- if (v > 0) round(log10(v)) else 0L

    nl      <- 80L %/% ne
    fmt_hdr <- sprintf("(%dF%d.0)", nl, ne)

    id_padded <- formatC(trimws(id), width = 40, flag = "-")
    hdr <- sprintf("%8d=N%8d=I%21s%-40s%3d%s~",
                   n, jyr, "", id_padded, nc, fmt_hdr)
    writeLines(hdr, con)

    lines_out <- character(0)
    for (start in seq(1, n, by = nl)) {
      chunk     <- izz[start:min(start + nl - 1L, n)]
      lines_out <- c(lines_out,
                     paste(formatC(chunk, width = ne, format = "d"),
                           collapse = ""))
    }
    writeLines(lines_out, con)
  }
}


# -----------------------------------------------------------------------------
# .write_tucson  —  write a Tucson .rwl file  [TRRW 'MW' branch equivalent]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Writes the xDPL internal series list to a Tucson measurement file.
#   Values are stored as F6.2 (0.01 mm precision).  Each series is terminated
#   with the standard sentinel value 9.99.  The first decade line of each
#   series is padded with blanks when the series does not begin on a decade
#   boundary, replicating the alignment behaviour of the original TRRW code.
#
# ARGUMENTS
#   series_list   List of series in xDPL internal format.
#   path          Output file path (overwritten if it exists).
#
.write_tucson <- function(series_list, path) {
  con <- file(path, "wt")
  on.exit(close(con))

  for (s in series_list) {
    id  <- formatC(trimws(s$id), width = 8, flag = "-")
    jyr <- s$yr_start
    y   <- s$values
    n   <- length(y)
    if (n == 0) next

    offset <- jyr %% 10
    pad    <- if (offset == 0) 0L else 10L - offset
    idx    <- 1L
    kyr    <- jyr - pad

    while (idx <= n) {
      if (idx == 1L && pad > 0) {
        n_line        <- min(10L - pad, n)
        prefix_blanks <- strrep("      ", pad)
      } else {
        n_line        <- min(10L, n - idx + 1L)
        prefix_blanks <- ""
      }

      vals_line <- y[idx:min(idx + n_line - 1L, n)]
      val_str   <- paste(formatC(vals_line, format = "f", digits = 2,
                                 width = 6), collapse = "")
      writeLines(sprintf("%s%4d%s%s", id, kyr, prefix_blanks, val_str), con)

      idx  <- idx + n_line
      kyr  <- kyr + 10L
      pad  <- 0L
      prefix_blanks <- ""
    }

    writeLines(sprintf("%s%4d  9.99", id, kyr), con)
  }
}


# -----------------------------------------------------------------------------
# .trdisp  —  print a ~50-year inspection window  [TRDISP equivalent]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Prints up to 50 consecutive years of ring-width values centred as closely
#   as possible on around_yr, respecting series boundaries.  Output is
#   formatted in rows of 10, matching the DPL console display, to facilitate
#   visual detection of outliers or crossdating errors before editing.
#
# ARGUMENTS
#   values      Numeric vector of ring-width values for one series.
#   yr_start    Integer. First calendar year of the series.
#   around_yr   Integer. Target year to centre the display on.
#
# RETURNS
#   Invisibly: named numeric vector of the displayed values (names = years).
#
.trdisp <- function(values, yr_start, around_yr) {
  n   <- length(values)
  lyr <- yr_start + n - 1L

  kyr <- max(yr_start + 18L, around_yr)
  kyr <- min(lyr - 18L, kyr)
  ka  <- max(yr_start, kyr - 25L)
  kz  <- min(lyr, ka + 49L)
  ka  <- max(kz - 49L, yr_start)

  years       <- seq(ka, kz)
  vals        <- values[years - yr_start + 1L]
  names(vals) <- years

  for (i in seq(1, length(years), by = 10)) {
    j <- min(i + 9L, length(years))
    cat(sprintf("%8d", years[i:j]), "\n")
    cat(sprintf("%8.2f", vals[i:j]), "\n")
  }
  invisible(vals)
}


#' Read a DPL ring-measurement file
#'
#' @description
#' Reads a ring-width measurement file in either DPL compact (`~`) or Tucson
#' (`.rwl`) format. Format is auto-detected by scanning the first 20 lines:
#' any line ending in `~` signals compact format; otherwise Tucson is assumed.
#'
#' @param path Character. Path to the input file.
#' @param format Character. One of \code{"auto"} (default), \code{"compact"} (Holmes
#'   internal `~` format), or \code{"tucson"} (ITRDB Tucson decade-per-line, F6.2).
#' @param as_rwl Logical. `TRUE` (default) returns a dplR `rwl` data.frame;
#'   `FALSE` returns the xDPL internal series list (`list(id, yr_start,
#'   values)`), useful for chaining into \code{\link{dpl_edt}} or \code{\link{dpl_write}}
#'   without a `rwl` round-trip.
#'
#' @return
#' `as_rwl = TRUE`: `c("rwl","data.frame")` with consecutive character years
#' as row names, one numeric column per series (`NA` outside each series'
#' span). Compatible with \code{\link[dplR]{time}}, \code{\link[dplR]{rwl.stats}},
#' \code{\link[dplR]{detrend}}, \code{\link[dplR]{chron}}, etc.
#'
#' `as_rwl = FALSE`: a plain list; each element: `$id` (character),
#' `$yr_start` (integer), `$values` (numeric, no leading/trailing `NA`).
#'
#' @examples
#' \dontrun{
#' # dpl_read handles Tucson .rwl files; CL-MIR.rwl uses Holmes decadal format,
#' # so use dpl_read_dec for this dataset.
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # Round-trip: write as Tucson then read back
#' tmp <- tempfile(fileext = ".rwl")
#' dpl_write(rwl, tmp, format = "tucson")
#' rwl2 <- dpl_read(tmp, format = "tucson")
#' all.equal(rwl, rwl2)
#' }
#'
#' @seealso \code{\link{dpl_write}}, \code{\link{dpl_edt}}, \code{\link{dpl_edit_file}}
#' @export
dpl_read <- function(path, format = "auto", as_rwl = TRUE) {
  if (!file.exists(path))
    stop("File not found: ", path)

  if (format == "auto") {
    lns       <- readLines(path, n = 20, warn = FALSE)
    has_tilde <- any(grepl("~$", lns))
    format    <- if (has_tilde) "compact" else "tucson"
  }

  sl <- switch(format,
    compact = .read_compact(path),
    tucson  = .read_tucson(path),
    stop("Unknown format: '", format, "'. Use 'compact', 'tucson', or 'auto'.")
  )

  if (as_rwl) .series_to_rwl(sl) else sl
}


#' Write a DPL ring-measurement file
#'
#' @description
#' Writes ring-width data to a DPL compact or Tucson file. Accepts either a
#' dplR `rwl` data.frame or the xDPL internal series list; input type is
#' detected automatically. This is the inverse of \code{\link{dpl_read}} and completes
#' the read -> edit -> write workflow.
#'
#' @param x A dplR `rwl` data.frame **or** a list of series in xDPL internal
#'   format (`list(id, yr_start, values)`).
#' @param path Character. Output file path. Existing files are overwritten.
#' @param format Character. \code{"compact"} (default, DPL `~` format) or
#'   \code{"tucson"} (ITRDB Tucson format, compatible with COFECHA, ARSTAN,
#'   CooRecorder, CDendro, and OpenDendro tools).
#'
#' @return Invisibly returns `path`, for use in pipelines.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # Write in Tucson format (compatible with COFECHA, ARSTAN, OpenDendro)
#' tmp_rw <- tempfile(fileext = ".rwl")
#' dpl_write(rwl, tmp_rw, format = "tucson")
#'
#' # Write in compact DPL format
#' tmp_dat <- tempfile(fileext = ".dat")
#' dpl_write(rwl, tmp_dat, format = "compact")
#' }
#'
#' @seealso \code{\link{dpl_read}}, \code{\link{dpl_edit_file}}
#' @export
dpl_write <- function(x, path, format = "compact") {
  sl <- if (.is_rwl(x)) .rwl_to_series(x) else x
  switch(format,
    compact = .write_compact(sl, path),
    tucson  = .write_tucson(sl,  path),
    stop("Unknown format: '", format, "'. Use 'compact' or 'tucson'.")
  )
  invisible(path)
}


#' Read a decadal-format ring-width file with any label length
#'
#' @description
#' Reads the Tucson/DPL compact decadal format where each line contains a
#' series label immediately followed by the decade year (no separator), then
#' up to 10 ring-width values. Unlike \code{\link[dplR]{read.rwl}}, makes no assumption
#' about label width: any label length (6, 8, 10+ characters) is handled by
#' parsing the last four characters of the first token as the year.
#'
#' @details
#' **Line format:** `<label><yyyy>  <v0>  <v1>  ...  <v9>`
#'
#' The end-of-series marker is `999` (Tucson standard). Multiple series per
#' file are supported (rows interleaved by decade). Also accepts a character
#' vector of lines instead of a file path.
#'
#' @param file Character. Path to the decadal format file **or** a character
#'   vector of lines (for inline data or testing).
#' @param label_length Integer or `NULL` (default). Maximum label width used
#'   only to disambiguate compact format (label+year concatenated) from spaced
#'   format (label and year as separate tokens). In spaced format, IDs shorter
#'   than `label_length` are accepted, so a single value covers files with
#'   mixed-length IDs. Supply an integer only when auto-detection fails because
#'   the label itself ends in four digits; leave `NULL` otherwise.
#' @param stop_val Integer (default `999L`). End-of-series sentinel.
#'   Set `NULL` to import all values as-is.
#' @param na_val Integer vector (default `c(-999L, -9999L)`). Values that mark
#'   a missing measurement \emph{within} a series (the year is kept, the value
#'   becomes `NA`); typical of files exported from measuring software for
#'   series with gaps, e.g. microcores not sampled in some years. Values equal
#'   to `stop_val` are always treated as end of series first. Set `NULL` to
#'   disable.
#' @param unit Character: \code{"auto"} (default), \code{"0.001mm"}, \code{"0.01mm"}, or
#'   \code{"mm"}. Under \code{"auto"}, two signals are combined: a negative `stop_val`
#'   (e.g. `-9999`) with median raw value > 20 implies 1/1000 mm (scale
#'   0.001); a positive or absent sentinel with median > 20 implies 1/100 mm
#'   (scale 0.01); otherwise values are taken as already in mm. The negative
#'   sentinel is the reliable marker of Holmes 1/1000 mm encoding, where
#'   `-9999` represents -9.999 mm.
#' @param skip_lines Integer (default `0L`) or \code{"auto"}. Number of header
#'   lines to skip. \code{"auto"} skips leading lines that do not start with a
#'   label+year token.
#' @param base_century Integer (default `1900L`). Century base for resolving
#'   century-relative year fields (e.g., \code{"0"} -> 1900, \code{"10"} -> 1910).
#'   Ignored for absolute 4-digit years (>= 1000).
#'
#' @return A `c("rwl","data.frame")` with integer row names (calendar years),
#'   one column per series (values in mm, `NA` for absent/missing rings).
#'   A one-line summary is printed to the console. Compatible with all dplR
#'   and cofechar functions.
#'
#' @examples
#' \dontrun{
#' # Read the Cerro Mirador Fitzroya dataset
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#' dim(rwl)                          # 597 years x 36 series
#' range(as.integer(rownames(rwl)))  # 1406 to 2002
#' head(rwl[, 1:4])
#'
#' # Read with explicit unit and skip check
#' rwl2 <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                      unit = "0.001mm", skip_lines = 0L)
#' identical(rwl, rwl2)
#' }
#'
#' @seealso \code{\link{dpl_read}}, \code{\link{dpl_cof}}
#' @export
dpl_read_dec <- function(file,
                           label_length = NULL,
                           stop_val     = 999L,
                           na_val       = c(-999L, -9999L),
                           unit         = c("auto", "0.001mm", "0.01mm", "mm"),
                           skip_lines   = 0L,
                           base_century = 1900L) {

  unit <- match.arg(unit)

  # ---- 1. Read lines --------------------------------------------------------
  if (length(file) > 1L || (length(file) == 1L && !file.exists(file))) {
    # Treat as character vector of lines directly
    raw_lines <- as.character(file)
  } else {
    raw_lines <- readLines(file, warn = FALSE)
  }

  # ---- 2. Skip header lines -------------------------------------------------
  if (identical(skip_lines, "auto") || identical(skip_lines, "AUTO")) {
    # Skip any leading line that is not a data line.
    # A data line has either: tokens[1] ending in 4 digits (compact),
    # or tokens[2] being a 1-4 digit year (spaced, including century-relative).
    is_data <- vapply(raw_lines, function(l) {
      tok <- strsplit(trimws(l), "\\s+")[[1L]]
      if (length(tok) == 0L) return(FALSE)
      grepl("[0-9]{1,4}$", tok[1L]) ||
        (length(tok) >= 2L && grepl("^[0-9]{1,4}$", tok[2L]))
    }, logical(1))
    first_data <- which(is_data)[1L]
    if (!is.na(first_data) && first_data > 1L)
      raw_lines <- raw_lines[first_data:length(raw_lines)]
  } else if (is.numeric(skip_lines) && skip_lines > 0L) {
    raw_lines <- raw_lines[seq(as.integer(skip_lines) + 1L, length(raw_lines))]
  }

  # Drop blank lines
  raw_lines <- raw_lines[nchar(trimws(raw_lines)) > 0L]

  if (length(raw_lines) == 0L) stop("No data lines found after header skip.")

  # ---- 3. Parse lines -------------------------------------------------------
  # Handles two sub-formats transparently:
  #   Compact: "BAR575AE1909   224   ..."  → tokens[1] = label+year concatenated
  #   Spaced:  "BAR575AE 1909   224   ..."  → tokens[1] = label, tokens[2] = year
  # Detection: if tokens[1] ends in 4 digits (and is long enough), it is compact;
  # otherwise the year is in tokens[2] and values start at tokens[3].

  data_store <- list()

  for (line in raw_lines) {
    tokens <- strsplit(trimws(line), "\\s+")[[1L]]
    if (length(tokens) < 2L) next

    tag <- tokens[1L]

    # Determine whether format is compact or spaced, then extract sid + decade_yr
    if (!is.null(label_length)) {
      ll <- as.integer(label_length)
      if (nchar(tag) == ll + 4L) {
        # Compact: "LABEL1909" — year embedded at the end of first token
        sid       <- substr(tag, 1L, ll)
        yr_str    <- substr(tag, ll + 1L, ll + 4L)
        val_strs  <- tokens[-1L]
      } else if (nchar(tag) <= ll && length(tokens) >= 2L &&
                 grepl("^[0-9]{4}$", tokens[2L])) {
        # Spaced: "LABEL 1909" — year is tokens[2]; ID may be shorter than ll
        sid       <- tag
        yr_str    <- tokens[2L]
        val_strs  <- tokens[-(1L:2L)]
      } else {
        stop(sprintf(
          paste0("Token '%s' (length %d) does not match label_length=%d: ",
                 "expected a spaced-format ID (<= %d chars) with a 4-digit ",
                 "year as the next token, or a compact-format token of %d chars."),
          tag, nchar(tag), ll, ll, ll + 4L))
      }
    } else {
      # Auto-detect: compact if last 4 chars of tokens[1] are all digits,
      # otherwise try tokens[2] as the year (spaced format).
      yr4 <- substr(tag, nchar(tag) - 3L, nchar(tag))
      if (nchar(tag) >= 5L && grepl("^[0-9]{4}$", yr4)) {
        # Compact
        sid       <- substr(tag, 1L, nchar(tag) - 4L)
        yr_str    <- yr4
        val_strs  <- tokens[-1L]
      } else if (length(tokens) >= 2L && grepl("^[0-9]{1,4}$", tokens[2L])) {
        # Spaced: label is tokens[1], year is tokens[2] (4-digit absolute or 1-3 digit relative)
        sid       <- tag
        yr_str    <- tokens[2L]
        val_strs  <- tokens[-(1L:2L)]
      } else {
        stop(sprintf(
          paste0("Cannot parse year from line starting with '%s': ",
                 "tokens[1] does not end in 4 digits and tokens[2] ('%s') ",
                 "is not a 4-digit year. Supply label_length explicitly."),
          tag, if (length(tokens) >= 2L) tokens[2L] else ""))
      }
    }

    if (!grepl("^-?[0-9]+$", trimws(yr_str)))
      stop(sprintf("Parsed year string '%s' is not numeric on line: %s",
                   yr_str, trimws(line)))

    yr_int <- as.integer(yr_str)

    # Resolve century-relative years: some files store only 2-digit decade offsets
    # (e.g., "0"=1900, "10"=1910, "80"=1980). If the year integer is < 1000,
    # treat it as an offset from base_century (default 1900).
    decade_yr <- if (yr_int < 1000L) as.integer(base_century) + yr_int else yr_int
    if (!(sid %in% names(data_store))) data_store[[sid]] <- numeric(0)

    for (i in seq_along(val_strs)) {
      v_int <- suppressWarnings(as.integer(val_strs[i]))
      if (is.na(v_int)) next
      yr <- decade_yr + (i - 1L)
      if (!is.null(stop_val) && v_int == as.integer(stop_val)) break
      if (!is.null(na_val) && v_int %in% as.integer(na_val)) v_int <- NA_integer_
      data_store[[sid]][as.character(yr)] <- v_int
    }
  }

  if (length(data_store) == 0L) stop("No series data found in file.")

  # ---- 4. Unit conversion ---------------------------------------------------
  # Collect all non-NA raw values to decide conversion
  all_raw <- unlist(lapply(data_store, as.numeric), use.names = FALSE)
  all_raw <- all_raw[!is.na(all_raw)]

  # Three-way auto-detection:
  #   negative stop_val (e.g. -9999) AND median > 20  ->  0.001 mm
  #     A negative end-of-series sentinel is the signature of Holmes
  #     1/1000 mm integer encoding (-9999 = -9.999 mm in that scale).
  #     Standard Tucson encoding uses a positive sentinel (999 or 9990).
  #   median > 20 (positive or NULL stop_val)          ->  0.01 mm
  #     Standard Tucson integer encoding (9999 ~ 99.99 mm acts as flag).
  #   otherwise                                         ->  mm
  scale <- switch(unit,
    "0.001mm" = 0.001,
    "0.01mm"  = 0.01,
    "mm"      = 1.0,
    "auto"    = {
      if (length(all_raw) == 0L) {
        1.0
      } else {
        med      <- stats::median(all_raw, na.rm = TRUE)
        neg_sent <- !is.null(stop_val) && as.integer(stop_val) < 0L
        if      (neg_sent && med > 20) 0.001
        else if (med > 20)             0.01
        else                           1.0
      }
    }
  )

  # Convert each series to mm
  for (sid in names(data_store)) {
    v <- data_store[[sid]]
    data_store[[sid]] <- setNames(as.numeric(v) * scale, names(v))
  }

  # ---- 5. Build rwl data.frame ----------------------------------------------
  # Union year range across all series
  all_yr_names <- unique(unlist(lapply(data_store, names)))
  all_yrs      <- sort(as.integer(all_yr_names))
  yr_chars     <- as.character(all_yrs)

  rwl <- as.data.frame(
    lapply(data_store, function(v) {
      out <- rep(NA_real_, length(all_yrs))
      idx <- match(names(v), yr_chars)
      ok  <- !is.na(idx)
      out[idx[ok]] <- v[ok]
      out
    }),
    row.names  = as.character(all_yrs),
    check.names = FALSE
  )

  # Attach dplR-compatible class and attributes
  class(rwl) <- c("rwl", "data.frame")
  attr(rwl, "series.ids") <- names(data_store)

  # Print a compact summary (matching dplR's read.rwl() style)
  n_ser <- ncol(rwl)
  cat(sprintf(
    "Read %d series spanning %d to %d (%d years)\n",
    n_ser, min(all_yrs), max(all_yrs), length(all_yrs)
  ))

  invisible(rwl)
}


#' Trim an rwl to the years covered by its series
#'
#' @description
#' Removes the leading and trailing years in which no series has a value, so
#' that the first row is the first year of the earliest series and the last
#' row the last year of the latest. Subsetting columns with \code{[} keeps
#' the full year axis of the original collection (rows of \code{NA} for the
#' years the selected series do not cover); \code{dpl_trim} restores a
#' frame that spans only the selected series. \code{\link{dpl_edt}} and
#' \code{\link{dpl_merge}} already return trimmed frames.
#'
#' @param rwl A dplR \code{rwl} data.frame.
#' @param series Optional: series to keep before trimming, by ID or
#'   position (a shortcut for \code{dpl_trim(rwl[, series])}).
#'
#' @return The \code{rwl} with only the years between the first and last
#'   non-missing value of any column; interior years are never removed.
#'   An all-\code{NA} input returns a 0-row frame.
#'
#' @examples
#' \dontrun{
#' sub <- rwl[, c("ACC026", "ACC026B")]   # still 1406-2002 rows, mostly NA
#' sub <- dpl_trim(sub)                    # now spans the two series only
#' sub <- dpl_trim(rwl, series = 1:5)      # same in one step
#' }
#'
#' @seealso \code{\link{dpl_edt}}, \code{\link{dpl_merge}}
#' @export
dpl_trim <- function(rwl, series = NULL) {
  if (!.is_rwl(rwl)) stop("'rwl' must be a dplR rwl data.frame.")
  if (!is.null(series))
    rwl <- rwl[, .resolve_series(series, colnames(rwl)), drop = FALSE]
  has <- rowSums(!is.na(rwl)) > 0L
  if (!any(has)) {
    out <- rwl[integer(0), , drop = FALSE]
  } else {
    r   <- range(which(has))
    out <- rwl[seq(r[1L], r[2L]), , drop = FALSE]
  }
  class(out) <- c("rwl", "data.frame")
  out
}


#' Merge two or more rwl objects onto a common year axis
#'
#' @description
#' Combines all series from a list of dplR `rwl` data.frames into a single
#' `rwl` object spanning the union of all year ranges. This is the in-memory
#' equivalent of concatenating multiple `.rwl` files before COFECHA or
#' building a site chronology from a multi-file collection.
#'
#' @param rwl_list A list of at least two dplR `rwl` data.frames (or plain
#'   `data.frame`s with integer-coercible row names and at least one column).
#' @param trim Logical (default \code{TRUE}). Drop leading and trailing
#'   years in which no series has a value, so the result spans only the
#'   merged series even when the inputs carry all-\code{NA} rows from an
#'   earlier column subset (see \code{\link{dpl_trim}}). \code{FALSE}
#'   keeps the union of the inputs' row ranges.
#' @param dup_action Character. Action when the same series ID appears in
#'   more than one input:
#'   \describe{
#'     \item{\code{"error"}}{(default) Stop with a message listing duplicate IDs.}
#'     \item{\code{"suffix"}}{Append `_2`, `_3`, ... to duplicate IDs, with a warning.}
#'     \item{\code{"keep"}}{Allow duplicate column names (may confuse dplR functions).}
#'   }
#'
#' @return A `c("rwl","data.frame")` spanning the full union year range, with
#'   `NA` outside each series' original span. Passes all dplR generics
#'   (\code{\link[dplR]{time}}, \code{\link[dplR]{rwl.stats}}, \code{\link[dplR]{detrend}}, \code{\link[dplR]{chron}})
#'   without further conversion.
#'
#' @note Measurement units are not checked for consistency across inputs.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # Split into two halves and merge back
#' half <- ncol(rwl) %/% 2L
#' rwl_a <- rwl[, 1:half]
#' rwl_b <- rwl[, (half + 1L):ncol(rwl)]
#' rwl_merged <- dpl_merge(list(rwl_a, rwl_b))
#' identical(colnames(rwl_merged), colnames(rwl))  # TRUE
#' }
#'
#' @seealso \code{\link{dpl_read}}, \code{\link{dpl_edt}}, \code{\link{dpl_write}}
#' @export
dpl_merge <- function(rwl_list, dup_action = "error", trim = TRUE) {

  # --- Input validation -----------------------------------------------------
  if (!is.list(rwl_list) || length(rwl_list) < 2L)
    stop("'rwl_list' must be a list of at least two rwl data.frames.")

  for (i in seq_along(rwl_list)) {
    x <- rwl_list[[i]]
    if (!is.data.frame(x))
      stop("Element ", i, " of 'rwl_list' is not a data.frame.")
    if (ncol(x) == 0L)
      stop("Element ", i, " of 'rwl_list' has no columns.")
    if (is.null(rownames(x)) || anyNA(suppressWarnings(as.integer(rownames(x)))))
      stop("Element ", i, " of 'rwl_list' does not have integer-coercible row names (years).")
  }

  if (!dup_action %in% c("error", "suffix", "keep"))
    stop("'dup_action' must be \"error\", \"suffix\", or \"keep\".")

  # --- Collect all series IDs and detect duplicates -------------------------
  all_ids <- unlist(lapply(rwl_list, colnames))

  dup_ids <- unique(all_ids[duplicated(all_ids)])
  if (length(dup_ids) > 0L) {
    if (dup_action == "error") {
      stop("Duplicate series IDs found across inputs: ",
           paste(dup_ids, collapse = ", "),
           ".\nUse dup_action = \"suffix\" to rename or \"keep\" to allow.")
    } else if (dup_action == "suffix") {
      # Assign unique suffixes to every occurrence of a duplicated ID
      # after the first: ABC001 stays, duplicates become ABC001_2, ABC001_3, ...
      # next_sfx[id] holds the integer suffix to apply to the *next* duplicate.
      next_sfx <- integer(0)
      new_ids  <- all_ids
      for (k in seq_along(all_ids)) {
        id <- all_ids[k]
        if (id %in% dup_ids) {
          if (is.na(next_sfx[id])) {
            # First occurrence of this duplicated ID — keep original name,
            # record that the next occurrence will get suffix _2.
            next_sfx[id] <- 2L
          } else {
            new_ids[k]   <- paste0(id, "_", next_sfx[id])
            next_sfx[id] <- next_sfx[id] + 1L
          }
        }
      }
      renamed <- new_ids != all_ids
      warning("Duplicate series IDs renamed: ",
              paste(all_ids[renamed], "->", new_ids[renamed], collapse = ", "))
      # Re-assign updated column names back to each rwl in the list
      ptr <- 1L
      for (i in seq_along(rwl_list)) {
        nc <- ncol(rwl_list[[i]])
        colnames(rwl_list[[i]]) <- new_ids[ptr:(ptr + nc - 1L)]
        ptr <- ptr + nc
      }
    }
    # dup_action == "keep": do nothing, R allows non-unique column names
  }

  # --- Determine global year range ------------------------------------------
  yr_min <- min(vapply(rwl_list, function(x) min(as.integer(rownames(x))), integer(1)))
  yr_max <- max(vapply(rwl_list, function(x) max(as.integer(rownames(x))), integer(1)))
  all_years <- as.character(seq(yr_min, yr_max))
  n_yr      <- length(all_years)

  # --- Allocate output frame ------------------------------------------------
  total_cols <- sum(vapply(rwl_list, ncol, integer(1)))
  final_ids  <- unlist(lapply(rwl_list, colnames))

  out <- as.data.frame(
    matrix(NA_real_, nrow = n_yr, ncol = total_cols,
           dimnames = list(all_years, final_ids))
  )

  # --- Fill each series into the correct rows --------------------------------
  for (rwl in rwl_list) {
    rn   <- rownames(rwl)
    # Only use years that fall within the global range (they all should, but
    # guard against manually constructed frames with out-of-range row names)
    keep <- rn[rn %in% all_years]
    out[keep, colnames(rwl)] <- rwl[keep, , drop = FALSE]
  }

  class(out) <- c("rwl", "data.frame")
  if (isTRUE(trim)) dpl_trim(out) else out
}


#' Apply editing operations to ring-measurement series
#'
#' @description
#' Programmatic, fully reproducible equivalent of Holmes' interactive DPLEDT
#' routine. Every correction DPLEDT accepted from the console --- rename, shift,
#' insert, delete, replace, trim, omit --- is expressed as a named list element
#' in `edits` and applied sequentially. Returns a dplR `rwl` data.frame ready
#' for any dplR function or for \code{\link{dpl_cof}}.
#'
#' @param x A dplR `rwl` data.frame **or** a list of series in xDPL internal
#'   format (`list(id, yr_start, values)`).
#' @param edits A list of edit instructions applied in order. Each is a named
#'   list with at minimum `series` (character ID or integer 1-based index) and
#'   `op` (character; see **Operations**). Additional fields depend on `op`.
#'   An empty list passes all series through according to `default_action`.
#' @param keep Series to keep, by ID or 1-based position in `x` (e.g.
#'   `c("ACC026", "ACC026B")` or `1:5`); all others are left out before any
#'   edit is applied. Input order is preserved. Cannot be combined with
#'   `drop`.
#' @param drop Series to leave out, by ID or position; all others are kept.
#'   Cannot be combined with `keep`.
#' @param default_action Character. \code{"copy"} (default) passes unmentioned series
#'   through; \code{"omit"} drops them. This is the DPLEDT way of extracting a
#'   subset (one `copy` edit per series to keep); `keep` and `drop` do the
#'   same in one argument.
#' @param as_rwl Logical. `TRUE` (default) returns a dplR `rwl` data.frame.
#'   `FALSE` returns the xDPL internal series list for passing to \code{\link{dpl_write}}.
#' @param verbose Logical. `TRUE` (default) prints one line per series with its
#'   ID, year range, length, and disposition. Set `FALSE` in batch scripts.
#'
#' @section Operations (`op` field):
#' \describe{
#'   \item{\code{"copy"}}{Pass through unchanged. Overrides \code{"omit"} for
#'     specific series. No additional fields.}
#'   \item{\code{"omit"}}{Drop the series. No additional fields.}
#'   \item{\code{"rename"}}{Change the series ID. Field: `new_id` (character).
#'     Subsequent edits in the same call must still use the original ID; use the
#'     1-based integer index to avoid ambiguity after renaming.}
#'   \item{\code{"first_year"}}{Reassign the first calendar year; span shifts rigidly.
#'     Field: `new_first_year` (integer).}
#'   \item{\code{"last_year"}}{Reassign the last calendar year; span shifts rigidly.
#'     Field: `new_last_year` (integer).}
#'   \item{\code{"replace"}}{Overwrite the ring width at one year. A value of 9.99
#'     (DPL absent-ring sentinel) is silently replaced with 10.0, matching
#'     Holmes' original behaviour. Fields: `year` (integer), `value` (numeric).}
#'   \item{\code{"insert"}}{Insert a new ring *before* a given year, increasing series
#'     length by 1. Fields: `year` (integer), `value` (numeric), `move`
#'     (\code{"back"} (default) shifts `first_year` back 1; \code{"forward"} shifts
#'     `last_year` forward 1).}
#'   \item{\code{"delete"}}{Remove the ring at one year, decreasing length by 1.
#'     Fields: `year` (integer), `move` (\code{"forward"} (default) shifts
#'     `first_year` forward 1; \code{"back"} shifts `last_year` back 1).}
#'   \item{\code{"trim_start"}}{Remove all rings before `first_year`.
#'     Field: `first_year` (integer).}
#'   \item{\code{"trim_end"}}{Remove all rings after `last_year`.
#'     Field: `last_year` (integer).}
#' }
#'
#' @return
#' `as_rwl = TRUE` (default): a dplR `rwl` data.frame ready for
#' \code{\link[dplR]{detrend}}, \code{\link[dplR]{chron}}, \code{\link[dplR]{rwl.stats}}, or \code{\link{dpl_cof}}.
#'
#' `as_rwl = FALSE`: a plain list of xDPL internal series for \code{\link{dpl_write}}.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # Look at the rings around a year before editing
#' dpl_display(rwl, series = "ACC014A", around_yr = 1953)
#'
#' ## --- Single operations ----------------------------------------------
#'
#' # Replace one value (a transposed reading: 0.056 entered, 0.065 measured)
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC014A", op = "replace", year = 1953, value = 0.065)
#' ))
#' rwl_ed["1953", "ACC014A"]
#'
#' # Redate a whole series from the inside: the innermost ring is 1768,
#' # not 1767. Every ring shifts by +1 (1767--1994 becomes 1768--1995).
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC014B", op = "first_year", new_first_year = 1768)
#' ))
#'
#' # Redate from the outside: the outer ring was formed in 2001, not 2000
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC013B", op = "last_year", new_last_year = 2001)
#' ))
#'
#' # Insert a ring missed on this radius before 1900 (a locally absent ring
#' # seen on another core). move = "back" (default) keeps the outer date and
#' # makes the inner part one year older (1775 -> 1774); move = "forward"
#' # keeps the pith date and moves the outer ring to 2001.
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC014A", op = "insert", year = 1900, value = 0.02)
#' ))
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC014A", op = "insert", year = 1900, value = 0.02,
#'        move = "forward")
#' ))
#'
#' # Delete a false ring counted at 1930. move = "forward" (default) keeps
#' # the outer date and makes the inner part one year younger (1775 -> 1776);
#' # move = "back" keeps the pith date and pulls the outer ring in to 1999.
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC014A", op = "delete", year = 1930)
#' ))
#'
#' # Discard rotten or unreadable rings at either end
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC14TC", op = "trim_end",   last_year  = 1940),
#'   list(series = "ACC026B", op = "trim_start", first_year = 1500)
#' ))
#'
#' # Rename a core. A later edit to the same series in the same call must
#' # use its position (or the OLD id), never the new id.
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC14TC", op = "rename",   new_id    = "ACC014C"),
#'   list(series = 21,        op = "trim_end", last_year = 1940)
#' ))
#' colnames(rwl_ed)[21]
#'
#' ## --- Several corrections in one reproducible call --------------------
#'
#' # The edit list is the whole record of what was changed: keep it in the
#' # script instead of overwriting the measurement file.
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC014A", op = "replace",    year = 1953, value = 0.065),
#'   list(series = "ACC014A", op = "delete",     year = 1930),
#'   list(series = "ACC014B", op = "first_year", new_first_year = 1768),
#'   list(series = "ACC14TC", op = "omit")
#' ), verbose = FALSE)
#'
#' # Trim the three longest series to the start of the well-replicated period
#' rwl_ed <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC026B", op = "trim_start", first_year = 1856),
#'   list(series = "ACC026",  op = "trim_start", first_year = 1856),
#'   list(series = "ACC026C", op = "trim_start", first_year = 1856)
#' ))
#'
#' ## --- Subsets: by name, by position, by exclusion ---------------------
#'
#' rwl_sub  <- dpl_edt(rwl, keep = c("ACC026", "ACC026B", "ACC026C"))
#' rwl_sub  <- dpl_edt(rwl, keep = 1:5)
#' rwl_sub  <- dpl_edt(rwl, drop = "ACC23TA")
#' long_ids <- colnames(rwl)[colSums(!is.na(rwl)) >= 300]
#' rwl_long <- dpl_edt(rwl, keep = long_ids)
#'
#' # keep and edits together: extract one core and correct it in one call.
#' # The result spans only that core's years (1775--2000), not the whole
#' # collection.
#' one <- dpl_edt(rwl, keep = "ACC014A", edits = list(
#'   list(series = "ACC014A", op = "replace", year = 1953, value = 0.065)
#' ))
#' range(as.integer(rownames(one)))
#'
#' ## --- Check the result and write it out -------------------------------
#'
#' dpl_display(rwl_ed, series = "ACC014A", around_yr = 1953)
#' cof <- dpl_cof(rwl_ed, parts = 7, verbose = FALSE)
#' cof$stats
#'
#' # Internal series list for dpl_write(), e.g. to feed another program
#' ser <- dpl_edt(rwl, edits = list(
#'   list(series = "ACC014A", op = "delete", year = 1930)
#' ), as_rwl = FALSE, verbose = FALSE)
#' dpl_write(ser, "CL-MIR_edited.rwl")
#' }
#'
#' @seealso \code{\link{dpl_display}}, \code{\link{dpl_edit_file}}, \code{\link{dpl_cof}}
#' @export
dpl_edt <- function(x,
                     edits          = list(),
                     keep           = NULL,
                     drop           = NULL,
                     default_action = "copy",
                     as_rwl         = TRUE,
                     verbose        = TRUE) {

  if (!default_action %in% c("copy", "omit"))
    stop('"default_action" must be "copy" or "omit".')
  if (!is.null(keep) && !is.null(drop))
    stop("Use either 'keep' or 'drop', not both.")

  if (.is_rwl(x)) {
    series_list <- .rwl_to_series(x)
  } else if (is.list(x)) {
    series_list <- x
  } else {
    stop('"x" must be a dplR data.frame or an xDPL series list.')
  }

  # ---- keep / drop: select series before any edit is applied ---------------
  # Names or 1-based positions in the input; input order is preserved
  # (keep = c(3, 1) still writes series 1 before series 3).
  if (!is.null(keep) || !is.null(drop)) {
    ids <- trimws(vapply(series_list, `[[`, character(1), "id"))
    sel <- if (!is.null(keep)) .resolve_series(keep, ids, "keep")
           else setdiff(ids, .resolve_series(drop, ids, "drop"))
    use <- ids %in% sel
    if (verbose)
      for (i in which(!use))
        message(sprintf("No %3d  %-8s  %d - %d (%d yr)  %s",
                        i, ids[i], series_list[[i]]$yr_start,
                        series_list[[i]]$yr_start + length(series_list[[i]]$values) - 1L,
                        length(series_list[[i]]$values),
                        if (!is.null(keep)) "NOT KEPT" else "DROPPED"))
    series_list <- series_list[use]
    if (length(series_list) == 0L) stop("No series left after keep/drop.")
  }

  n_ser      <- length(series_list)
  ops_by_seq <- vector("list", n_ser)

  for (ed in edits) {
    sel <- ed$series
    if (is.numeric(sel)) {
      idx <- as.integer(sel)
    } else {
      ids <- vapply(series_list, `[[`, character(1), "id")
      idx <- which(trimws(ids) == trimws(sel))
      if (length(idx) == 0)
        stop("Series '", sel, "' not found.")
      if (length(idx) > 1)
        warning("Series ID '", sel, "' matched ", length(idx),
                " series; applying edit to all matches.")
    }
    for (i in idx)
      ops_by_seq[[i]] <- c(ops_by_seq[[i]], list(ed))
  }

  output   <- list()
  n_copied <- 0L

  for (seq_no in seq_len(n_ser)) {
    s    <- series_list[[seq_no]]
    id   <- s$id
    y    <- s$values
    iyr  <- s$yr_start
    n    <- length(y)
    lyr  <- iyr + n - 1L
    orig <- s
    ops  <- ops_by_seq[[seq_no]]

    if (length(ops) == 0) {
      if (default_action == "omit") {
        if (verbose)
          message(sprintf("No %3d  %-8s  %d - %d (%d yr) -- OMITTED (default)",
                          seq_no, trimws(id), iyr, lyr, n))
        next
      }
      n_copied <- n_copied + 1L
      output[[n_copied]] <- s
      if (verbose)
        message(sprintf("No %3d  %-8s  %d - %d (%d yr)  COPIED",
                        n_copied, trimws(id), iyr, lyr, n))
      next
    }

    omit_flag <- FALSE

    for (ed in ops) {
      op <- ed$op

      if (op == "omit") {
        omit_flag <- TRUE
        break

      } else if (op == "copy") {
        next

      } else if (op == "rename") {
        if (is.null(ed$new_id))
          stop("'rename' requires 'new_id'.")
        id <- trimws(ed$new_id)

      } else if (op == "first_year") {
        if (is.null(ed$new_first_year))
          stop("'first_year' requires 'new_first_year'.")
        iyr <- as.integer(ed$new_first_year)
        lyr <- iyr + n - 1L

      } else if (op == "last_year") {
        if (is.null(ed$new_last_year))
          stop("'last_year' requires 'new_last_year'.")
        lyrn <- as.integer(ed$new_last_year)
        iyr  <- iyr - lyr + lyrn
        lyr  <- lyrn

      } else if (op == "replace") {
        if (is.null(ed$year) || is.null(ed$value))
          stop("'replace' requires 'year' and 'value'.")
        ichg <- as.integer(ed$year)
        if (ichg < iyr || ichg > lyr) {
          warning("replace: year ", ichg, " outside [", iyr, ",", lyr, "]. Skipped.")
          next
        }
        val <- ed$value
        if (round(val * 100) == 999) val <- 10.0
        y[ichg - iyr + 1L] <- val

      } else if (op == "insert") {
        if (is.null(ed$year) || is.null(ed$value))
          stop("'insert' requires 'year' and 'value'.")
        ichg <- as.integer(ed$year)
        if (ichg < iyr || ichg > lyr + 1L) {
          warning("insert: year ", ichg, " outside valid range [",
                  iyr, ",", lyr + 1L, "]. Skipped.")
          next
        }
        val  <- ed$value
        if (round(val * 100) == 999) val <- 10.0
        move <- if (!is.null(ed$move)) ed$move else "back"
        pos  <- ichg - iyr + 1L
        y    <- c(y[seq_len(pos - 1L)], val, y[pos:n])
        n    <- n + 1L
        if (move == "back") iyr <- iyr - 1L
        lyr  <- iyr + n - 1L

      } else if (op == "delete") {
        if (is.null(ed$year))
          stop("'delete' requires 'year'.")
        ichg <- as.integer(ed$year)
        if (ichg < iyr || ichg > lyr) {
          warning("delete: year ", ichg, " outside [", iyr, ",", lyr, "]. Skipped.")
          next
        }
        move <- if (!is.null(ed$move)) ed$move else "forward"
        pos  <- ichg - iyr + 1L
        y    <- y[-pos]
        n    <- n - 1L
        if (move == "forward") iyr <- iyr + 1L
        lyr  <- iyr + n - 1L

      } else if (op == "trim_start") {
        if (is.null(ed$first_year))
          stop("'trim_start' requires 'first_year'.")
        kyr <- as.integer(ed$first_year)
        if (kyr < iyr || kyr > lyr) {
          warning("trim_start: year ", kyr, " outside range. Skipped.")
          next
        }
        k   <- kyr - iyr
        y   <- y[(k + 1L):n]
        n   <- n - k
        iyr <- kyr
        lyr <- iyr + n - 1L

      } else if (op == "trim_end") {
        if (is.null(ed$last_year))
          stop("'trim_end' requires 'last_year'.")
        kyr <- as.integer(ed$last_year)
        if (kyr < iyr || kyr > lyr) {
          warning("trim_end: year ", kyr, " outside range. Skipped.")
          next
        }
        n   <- n - lyr + kyr
        lyr <- kyr
        y   <- y[seq_len(n)]

      } else {
        warning("Unknown op '", op, "' for series '", id, "'. Skipped.")
      }
    }

    if (omit_flag) {
      if (verbose)
        message(sprintf("No %3d  %-8s  %d - %d (%d yr)  OMITTED",
                        seq_no, trimws(orig$id), orig$yr_start,
                        orig$yr_start + length(orig$values) - 1L,
                        length(orig$values)))
      next
    }

    n_copied <- n_copied + 1L
    output[[n_copied]] <- list(id = id, yr_start = iyr, values = y)

    if (verbose)
      message(sprintf("No %3d  %-8s  %d - %d (%d yr)  COPIED",
                      n_copied, trimws(id), iyr, lyr, n))
  }

  message(sprintf("\n%d edited series returned.", n_copied))

  if (as_rwl) .series_to_rwl(output) else output
}


#' Print a ~50-year inspection window for one series
#'
#' @description
#' Prints up to 50 consecutive years of ring-width values centred on
#' `around_yr`, formatted in rows of 10 matching the DPL console display
#' (TRDISP). Useful for visually inspecting a series around a suspected
#' crossdating error before deciding which \code{\link{dpl_edt}} operation to apply.
#'
#' Accepts three input forms: (a) a multi-column `rwl` with `series` selector,
#' (b) a single-column `rwl` (`series` may be omitted), or (c) a single xDPL
#' internal series element as returned by `dpl_read(as_rwl = FALSE)[[i]]`.
#'
#' @param x A dplR `rwl` data.frame **or** a single xDPL internal series list
#'   (`list(id, yr_start, values)`).
#' @param series Character (column name) or integer (1-based index). Required
#'   when `x` has more than one column.
#' @param around_yr Integer. Year to centre the window on. Defaults to the
#'   series midpoint.
#'
#' @return Invisibly: a named numeric vector of the displayed values
#'   (names = years). Side effect: formatted output printed to the console.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # Inspect the 1950s for a series --- a known climatically active decade
#' dpl_display(rwl, series = "ACC014A", around_yr = 1950)
#'
#' # Centre on the series midpoint (default)
#' dpl_display(rwl, series = "ACC026B")
#'
#' # By column index
#' dpl_display(rwl, series = 3L, around_yr = 1800)
#' }
#'
#' @seealso \code{\link{dpl_edt}}, \code{\link{dpl_cof_diag}}
#' @export
dpl_display <- function(x, series = NULL, around_yr = NULL) {

  if (.is_rwl(x)) {
    if (is.null(series)) {
      if (ncol(x) > 1)
        stop("Supply 'series' to select one column from a multi-column rwl.")
      series <- 1L
    }
    col <- if (is.numeric(series)) as.integer(series) else
             match(series, colnames(x))
    if (is.na(col) || col < 1L || col > ncol(x))
      stop("'series' not found in data.frame.")
    sl_s <- .rwl_to_series(x[, col, drop = FALSE])[[1L]]
  } else if (is.list(x) && !is.null(x$yr_start)) {
    sl_s <- x
  } else {
    stop('"x" must be a dplR rwl data.frame or a single xDPL series list.')
  }

  y    <- sl_s$values
  iyr  <- sl_s$yr_start
  lyr  <- iyr + length(y) - 1L

  if (is.null(around_yr))
    around_yr <- round((iyr + lyr) / 2)

  cat(sprintf("\nSeries: %-8s   %d - %d  (%d yr)\n",
              trimws(sl_s$id), iyr, lyr, length(y)))
  invisible(.trdisp(y, iyr, around_yr))
}


#' Read, edit, and optionally write a ring-measurement file in one call
#'
#' @description
#' Convenience wrapper combining \code{\link{dpl_read}}, \code{\link{dpl_edt}}, and
#' \code{\link{dpl_write}}. Three common patterns:
#'
#' **(A) In-memory only** (`output_path = NULL`, default): read, edit, return
#' an `rwl` object. No file is written.
#'
#' **(B) Edit and archive** (`output_path` supplied): read, edit, write, and
#' return the edited `rwl` invisibly.
#'
#' **(C) Format conversion** (`edits = list()`, `output_path` supplied): read
#' in one format, write in another without any edits.
#'
#' @param path Character. Path to the input file.
#' @param output_path Character or `NULL`. Output file path. If `NULL`
#'   (default), no file is written and the result is returned visibly.
#' @param edits List of edit instructions; see \code{\link{dpl_edt}} for all operations.
#' @param keep,drop Series to keep or to leave out, by ID or position; see
#'   \code{\link{dpl_edt}}. The direct way to copy a subset of samples to a
#'   new file.
#' @param format Character. Input format: \code{"auto"} (default), \code{"compact"},
#'   or \code{"tucson"}.
#' @param output_format Character or `NULL`. Output format: \code{"compact"} or
#'   \code{"tucson"}. If `NULL` (default), mirrors the detected input format.
#' @param default_action Character. \code{"copy"} (default) or \code{"omit"}.
#' @param as_rwl Logical. `TRUE` (default) returns a dplR `rwl` data.frame.
#' @param verbose Logical. Print per-series edit log. Default `TRUE`.
#'
#' @return The edited data as a dplR `rwl` or xDPL internal series list.
#'   Returned *visibly* when `output_path = NULL`; *invisibly* otherwise.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#'
#' # (A) In-memory only: read, apply one correction, keep result in R
#' rwl_ed <- dpl_edit_file(mir_file,
#'   edits = list(
#'     list(series = "ACC014A", op = "replace", year = 1953, value = 0.26)
#'   ),
#'   format = "tucson"   # CL-MIR.rwl is Tucson-like decadal
#' )
#'
#' # (B) Edit and write corrected file
#' tmp <- tempfile(fileext = ".rwl")
#' dpl_edit_file(mir_file, output_path = tmp,
#'   edits = list(list(series = "ACC014A", op = "trim_start", first_year = 1856)),
#'   format = "tucson", output_format = "tucson"
#' )
#'
#' # (C) Copy a subset of samples to a new file
#' dpl_edit_file(mir_file, output_path = tempfile(fileext = ".rwl"),
#'               keep = c("ACC026", "ACC026B", "ACC026C"), format = "tucson")
#' dpl_edit_file(mir_file, output_path = tempfile(fileext = ".rwl"),
#'               keep = 1:5, format = "tucson")
#' dpl_edit_file(mir_file, output_path = tempfile(fileext = ".rwl"),
#'               drop = "ACC23TA", format = "tucson")
#' }
#'
#' @seealso \code{\link{dpl_read}}, \code{\link{dpl_edt}}, \code{\link{dpl_write}}
#' @export
dpl_edit_file <- function(path,
                            output_path    = NULL,
                            edits          = list(),
                            keep           = NULL,
                            drop           = NULL,
                            format         = "auto",
                            output_format  = NULL,
                            default_action = "copy",
                            as_rwl         = TRUE,
                            verbose        = TRUE) {

  if (format == "auto") {
    lns <- readLines(path, n = 20, warn = FALSE)
    fmt <- if (any(grepl("~$", lns))) "compact" else "tucson"
  } else {
    fmt <- format
  }

  if (is.null(output_format)) output_format <- fmt

  sl     <- dpl_read(path, format = fmt, as_rwl = FALSE)
  edited <- dpl_edt(sl,
                     edits          = edits,
                     keep           = keep,
                     drop           = drop,
                     default_action = default_action,
                     as_rwl         = FALSE,
                     verbose        = verbose)

  if (!is.null(output_path))
    dpl_write(edited, output_path, format = output_format)

  result <- if (as_rwl) .series_to_rwl(edited) else edited

  if (is.null(output_path)) result else invisible(result)
}


# =============================================================================
# =============================================================================
# COF MODULE
# COFECHA — Quality control and crossdating of tree-ring measurement series
# =============================================================================
# =============================================================================


# =============================================================================
# COF internal helpers  (direct translations of Holmes' Fortran routines)
# =============================================================================
#
# All helpers below are prefixed .cof_ and are not exported.  They mirror
# the original Fortran subroutines as closely as possible so that numerical
# output is equivalent to the DPL reference implementation.  Replacements
# with dplR or base-R equivalents can be made by swapping individual helpers
# without touching the main engine.


# -----------------------------------------------------------------------------
# .cof_suppress_zero  —  suppress leading zero in Fortran Fw.d output
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Many Fortran compilers suppress the leading zero for values |v| < 1 in
#   F-format output, giving '  .515' instead of '  0.515' and ' -.707' instead
#   of ' -0.707'.  This function replicates that behaviour in R without using
#   the perl=TRUE sub() function-replacement form (which requires R ≥ 4.0 and
#   causes 'cannot coerce type closure' errors on older builds).
#
# ARGUMENTS
#   s   Character scalar: a sprintf Fw.d formatted number string.
#   w   Integer. Total field width (used to re-justify after edit).
#
# RETURNS
#   Character scalar with leading zero removed and re-padded to width w.
#
.cof_suppress_zero <- function(s, w) {
  s <- trimws(s, "left")
  # Positive: '0.xxx' → '.xxx'
  if (grepl("^0\\.", s))
    s <- sub("^0\\.", ".", s)
  # Negative: '-0.xxx' → '-.xxx'
  else if (grepl("^-0\\.", s))
    s <- sub("^-0\\.", "-.", s)
  # Re-pad with leading spaces to width w
  formatC(s, width = w, flag = " ")
}

# Shorthand formatters replicating Fortran Fw.d with leading-zero suppression

.cof_f73 <- function(v)  .cof_suppress_zero(sprintf("%7.3f", v),  7L)
.cof_f83 <- function(v)  .cof_suppress_zero(sprintf("%8.3f", v),  8L)
.cof_f63 <- function(v)  .cof_suppress_zero(sprintf("%6.3f", v),  6L)
.cof_f52 <- function(v)  .cof_suppress_zero(sprintf("%5.2f", v),  5L)
.cof_f72 <- function(v)  sprintf("%7.2f", v)   # values > 1 — no leading zero issue


# .resolve_series  --  character IDs or 1-based positions -> character IDs
#   pool: the IDs that positions index into (e.g. colnames(rwl)); unknown
#   names warn and are dropped, out-of-range positions are an error, the
#   order of 'sel' is kept.
.resolve_series <- function(sel, pool, what = "series") {
  if (is.null(sel)) return(pool)
  if (is.numeric(sel)) {
    idx <- as.integer(sel)
    bad <- idx[is.na(idx) | idx < 1L | idx > length(pool)]
    if (length(bad))
      stop(sprintf("'%s' positions out of range (1..%d): %s",
                   what, length(pool), paste(bad, collapse = ", ")))
    return(pool[idx])
  }
  sel  <- as.character(sel)
  miss <- setdiff(sel, pool)
  if (length(miss))
    warning(sprintf("'%s' not found and ignored: %s", what, paste(miss, collapse = ", ")))
  sel[sel %in% pool]
}

# .cof_nint  --  Fortran NINT: round half AWAY from zero
#   R's round() is IEC 60559 round-half-to-even, so round(80.5) == 80 while
#   Fortran NINT(80.5) == 81.  BARPL's decile indices hit exact .5 for some
#   series lengths (e.g. NY = 161 -> 161*11/22 = 80.5), so the two disagree.
.cof_nint <- function(x) sign(x) * floor(abs(x) + 0.5)

# .cof_barpl_cuts  --  the 10 decile cut-points of BARPL (Fortran RANKRI step)
#   z_norm: normalised values.  Returns Z(1..10).
.cof_barpl_cuts <- function(z_norm) {
  NY  <- length(z_norm)
  zs  <- sort(z_norm, decreasing = TRUE)
  DEC <- NY / 11.0
  Z   <- numeric(10L)
  for (j in 1L:10L) {
    j1   <- max(1L, min(NY, .cof_nint(DEC * (j - 0.5))))
    j2   <- max(1L, min(NY, .cof_nint(DEC * (j + 0.5))))
    Z[j] <- 0.5 * (zs[j1] + zs[j2])
  }
  Z
}

# .cof_barpl_car  --  one 16-char BARPL cell for a normalised value
.cof_barpl_car <- function(yr, yn, Z) {
  lb <- 16L
  for (k in 1L:10L) if (yn < Z[k]) lb <- 16L - k
  lb <- max(6L, min(16L, lb))
  lp <- .cof_nint(yn * 4)
  if (lp < 0) { lp <- 96 - lp; if (lp > 122) lp <- 60 }   # a-z, '<'
  else        { lp <- lp + 64; if (lp > 90)  lp <- 62 }   # A-Z, '>'
  formatC(paste0(formatC(yr, width = 5L), strrep("-", lb - 6L), intToUtf8(lp)),
          width = 16L, flag = "-")
}

# .cof_barpl_pages  --  BARPL page layout, verified line-for-line against the
#   DPL COFECHA 4.04P PEL benchmark (PELCOF.OUT):
#     * pages of 400 years anchored at (JYR/400)*400; 8 columns x 50 rows
#     * row i (0..49) holds years IA+i, IA+50+i, ..., IA+350+i
#     * column header directly followed by the first row (no blank line)
#     * decade separator after rows 9, 19, 29, 39 (not 49).  DPL prints
#       "  " + 8 x " ----"; cofechar deliberately prints a blank line instead,
#       matching the decade spacing used in Part 3 (the only intentional
#       departure from the DPL Part 4 layout).
#     * one blank line closing each page
#     * page_hdr (character vector) repeated at the top of every page
#   make_car(yr) must return a 16-char cell or 16 spaces outside the span.
.cof_barpl_pages <- function(make_car, jyr, lyr, page_hdr) {
  col_hdr <- paste0("  ", paste(rep(" Year Rel value ", 8L), collapse = ""))
  sep_row <- ""   # DPL: "  " + 8 x " ----"; see header comment
  out <- character(0)
  IA  <- (jyr %/% 400L) * 400L
  if (IA > jyr) IA <- IA - 400L
  repeat {
    out <- c(out, page_hdr, col_hdr)
    for (i in 0L:49L) {
      cells <- vapply(0L:7L, function(jj) make_car(IA + i + jj * 50L), character(1))
      out <- c(out, paste0("  ", paste(cells, collapse = "")))
      if (i %% 10L == 9L && i < 49L) out <- c(out, sep_row)
    }
    out <- c(out, "")
    IA <- IA + 400L
    if (IA > lyr) break
  }
  sub("\\s+$", "", out)
}

# -----------------------------------------------------------------------------
# .cof_normts  —  convert array to mean = 0, variance = 1  [NORMTS]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Normalises a numeric vector in-place to mean = 0, variance = 1.
#   Translates Holmes' NORMTS subroutine exactly, including the convention
#   of dividing by the population standard deviation (denominator N, not N-1)
#   when k = 0 (the default used throughout COFECHA).
#
# ARGUMENTS
#   x   Numeric vector to normalise (modified in place semantics via return).
#   k   Integer. 0 = divide by N (population SD, Holmes default);
#               1 = divide by N-1 (sample SD).
#
# RETURNS
#   List with:
#     $z    Numeric vector, normalised.
#     $mean Numeric scalar, original mean.
#     $sd   Numeric scalar, original standard deviation (using chosen k).
#
.cof_normts <- function(x, k = 0L) {
  n  <- length(x)
  xm <- mean(x)
  sd <- sqrt(max(0, (sum((x - xm)^2)) / (n - k)))
  z  <- if (sd > 0) (x - xm) / sd else x - xm
  list(z = z, mean = xm, sd = sd)
}


# -----------------------------------------------------------------------------
# .cof_logtr  —  log-transform with additive constant  [LOGTR]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Adds a small constant (c * mean) to all values, takes natural log, then
#   subtracts log(constant) — effectively log-transforming departures from
#   the mean while preserving the scale of variation.  Translates LOGTR.
#
# ARGUMENTS
#   x   Numeric vector. All values should be non-negative after the constant
#       is applied; negative values are first shifted to zero.
#   c   Numeric. Proportion of the mean added before the transform (Holmes
#       default: 1/3 = 0.3333).
#
# RETURNS
#   Numeric vector of the same length as x, log-transformed.
#
.cof_logtr <- function(x, c = 1/3) {
  # Shift if negative values are present
  ym  <- min(x)
  if (ym < 0) x <- x - ym
  yb  <- mean(x)
  cc  <- c * yb
  if (cc <= 1e-4) cc <- 1.0
  cl  <- log(cc)
  log(x + cc) - cl
}


# -----------------------------------------------------------------------------
# .cof_divser  —  divide one series by another  [DIVSER]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Divides x by y element-wise to produce indices.  If y contains zeros or
#   negative values a constant is added temporarily to keep y strictly
#   positive before dividing, then removed from x to restore the original
#   scale.  Translates Holmes' DIVSER subroutine.
#
# ARGUMENTS
#   x   Numeric vector (measurement or filtered series).
#   y   Numeric vector of the same length (spline curve).
#
# RETURNS
#   Numeric vector: x / y (adjusted as above where necessary).
#
.cof_divser <- function(x, y) {
  ym <- min(y)
  if (ym <= 0) {
    adj <- 0.25 - ym
    x   <- x + adj
    y   <- y + adj
  } else {
    adj <- 0
  }
  z <- x / y
  if (adj > 0) {
    x <- x - adj
  }
  z
}


# -----------------------------------------------------------------------------
# .cof_varsta  —  variance stabilisation by spline on absolute values  [VARSTA]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Stabilises variance by fitting a cubic smoothing spline to the absolute
#   values of the normalised series and dividing.  The original sign of each
#   departure is restored, and the original mean and SD are re-applied.
#   Translates Holmes' VARSTA subroutine.
#
# ARGUMENTS
#   x    Numeric vector (filtered, normalised series).
#   ls   Integer. Spline rigidity in years (>0) or negative percent of N.
#
# RETURNS
#   Numeric vector of same length as x, variance-stabilised.
#
.cof_varsta <- function(x, ls = 32L) {
  n   <- length(x)
  nr  <- .cof_normts(x, k = 0L)
  z   <- nr$z
  xm  <- nr$mean
  sd  <- nr$sd
  sgn <- ifelse(z < 0, -1L, 1L)
  az  <- abs(z)
  lsp <- if (ls > 0) ls else max(1L, round(n * abs(ls) * 0.01))
  sp  <- .cof_spline(az, lsp, 0.5)
  az2 <- .cof_divser(az, sp)
  z2  <- az2 * sgn
  z2 * sd + xm
}


# -----------------------------------------------------------------------------
# .cof_spline  —  cubic smoothing spline  [SPLINE]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Fits a cubic smoothing spline to a time series using the Lagrange
#   multiplier approach of Cook & Peters (1981) as implemented by Holmes in
#   SPLINE.  The rigidity parameter (zsp) specifies the 50%-frequency-
#   response wavelength in years; at that wavelength the spline transmits
#   exactly pvp (default 0.5) of the variance.  The resulting spline is
#   stored in a banded LDLT factorisation and solved via forward-back
#   substitution, matching the numerical behaviour of Holmes' Fortran code.
#
# ARGUMENTS
#   x    Numeric vector. Time series to be smoothed.
#   zsp  Numeric. Spline rigidity: 50%-response wavelength in years.
#        Negative values are interpreted as percent of series length.
#   pvp  Numeric in (0,1). Proportion of variance at wavelength zsp that
#        the spline should pass.  Holmes default: 0.5 (i.e. 50%).
#
# RETURNS
#   Numeric vector of same length as x, the fitted spline curve.
#   If the matrix is not positive-definite, returns a flat line at the mean.
#
.cof_spline <- function(x, zsp, pvp = 0.5) {
  n <- length(x)
  if (n < 4L) return(rep(mean(x), n))

  if (zsp < 0) zsp <- max(1, n * abs(zsp) * 0.01)

  # Lagrange multiplier: defines the frequency response
  pi2   <- 2 * pi
  omega <- pi2 / zsp
  pd    <- ((1 / (1 - pvp) - 1) * 6 * (cos(omega) - 1)^2) /
             (cos(omega) + 2)

  nm2 <- n - 2L
  c1  <- c(1, -4, 6, -2)
  c2  <- c(0, 1/3, 4/3)

  # Build banded symmetric matrix A (nm2 x 3 stored columns)
  A <- matrix(0, nrow = nm2, ncol = 4L)
  for (i in seq_len(nm2)) {
    A[i, 1L] <- c1[1L] + pd * c2[1L]
    A[i, 2L] <- c1[2L] + pd * c2[2L]
    A[i, 3L] <- c1[3L] + pd * c2[3L]
    A[i, 4L] <- x[i] + c1[4L] * x[i + 1L] + x[i + 2L]
  }
  A[1L, 1L] <- c2[1L]
  A[1L, 2L] <- c2[1L]
  if (nm2 >= 2L) A[2L, 1L] <- c2[1L]

  nc   <- 2L
  ncp1 <- nc + 1L
  rn   <- 1 / (nm2 * 16.0)

  # LDLT factorisation (translated from Holmes' LUDAPB / LUELPB)
  # -- forward sweep
  for (i in seq_len(nm2)) {
    i1 <- max(1L, 1L - (i - ncp1))
    for (j in i1:ncp1) {
      l   <- i - ncp1 + j
      i2  <- ncp1 - j
      s   <- A[i, j]
      jm1 <- j - 1L
      if (jm1 > 0L && l >= 1L) {
        for (kk in seq_len(jm1)) {
          m <- i2 + kk
          if (l <= nm2 && m <= ncp1)
            s <- s - A[i, kk] * A[l, m]
        }
      }
      if (j == ncp1) {
        if (A[i, j] + s * rn <= A[i, j]) {
          # not positive definite — return flat line
          return(rep(mean(x), n))
        }
        A[i, j] <- 1 / sqrt(s)
      } else {
        if (l >= 1L && l <= nm2)
          A[i, j] <- s * A[l, ncp1]
      }
    }
  }

  # -- forward substitution LY = b
  iw <- 0L; l <- 0L
  for (i in seq_len(nm2)) {
    s <- A[i, 4L]
    if (nc > 0L) {
      if (iw != 0L) {
        l  <- min(l + 1L, nc)
        kl <- i - l
        for (j in (ncp1 - l):nc) {
          s  <- s - A[kl, 4L] * A[i, j]
          kl <- kl + 1L   # Fortran: KL=KL+1 inside the J loop
        }
      } else {
        if (s != 0) iw <- 1L
      }
    }
    A[i, 4L] <- s * A[i, ncp1]
  }

  # -- back substitution UX = Y
  # Fortran DO loops with an empty range do not execute; R's a:b counts
  # backwards, so guard the two loops below for nm2 < 3 (series of 4 years).
  A[nm2, 4L] <- A[nm2, 4L] * A[nm2, ncp1]
  for (i in seq_len(nm2)[-1L]) {
    kk  <- nm2 + 1L - i
    s   <- A[kk, 4L]
    kl  <- kk + 1L
    k1  <- min(nm2, kk + nc)
    lv  <- 1L
    for (j in kl:k1) {
      s  <- s - A[j, 4L] * A[j, ncp1 - lv]
      lv <- lv + 1L
    }
    A[kk, 4L] <- s * A[kk, ncp1]
  }

  # Reconstruct spline from second differences
  f <- numeric(n)
  for (i in seq_len(nm2)[-(1:2)])
    f[i] <- A[i - 2L, 4L] + c1[4L] * A[i - 1L, 4L] + A[i, 4L]
  f[1L] <- A[1L, 4L]
  f[2L] <- c1[4L] * A[1L, 4L] + A[2L, 4L]
  f[n - 1L] <- A[nm2 - 1L, 4L] + c1[4L] * A[nm2, 4L]
  f[n]       <- A[nm2, 4L]

  x - f
}


# -----------------------------------------------------------------------------
# .cof_mempr  —  Burg autoregressive modelling, first-min AIC  [MEMPR]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Fits an autoregressive model to a zero-mean time series using the Burg
#   maximum-entropy method (MEM) and selects the model order by the first
#   minimum of the Akaike Information Criterion (IOPT = 1 in Holmes' code).
#   Returns the AR residuals (prewhitened series) with the fitted mean
#   restored (+1 to give indices near 1.0).
#
#   This is a direct translation of Holmes' MEMPR subroutine (Cook/Holmes
#   1984 algorithm), not R's ar().  The Burg recursion, AIC computation and
#   stopping rule are numerically equivalent to the Fortran original.
#
# ARGUMENTS
#   x    Numeric vector. The (filtered, normalised) time series to model.
#        Mean is removed internally before fitting and restored on output.
#   lg   Integer. Maximum AR order to try (Holmes default: 10).
#
# RETURNS
#   List with:
#     $residuals  Numeric vector of AR residuals (length n), mean restored.
#     $order      Integer. Selected AR order (0 if AR modelling not helpful).
#     $aic        Numeric vector of AIC values for orders 1..lg.
#     $phi        Numeric vector of autocorrelation coefficients.
#
.cof_mempr <- function(x, lg = 10L) {
  # Exact translation of Fortran MEMPR (Burg maximum-entropy AR modelling).
  # Uses the LATTICE recursion with PEF/PER arrays exactly as Fortran does,
  # NOT the direct computation — these are theoretically equivalent but give
  # different numerical results due to floating-point accumulation order.
  #
  # Fortran MEMPR signature:
  #   SUBROUTINE MEMPR(M, F, H, PEF, PER, LG, SVE, PHI, IP, AIC, PM, NT, IOPT, RES, MXY)
  # Called with IOPT=1 (first-minimum AIC stopping).
  # RES = WK (whole array), so WK(I,1) = residuals = RES(I).

  n   <- length(x)
  rn  <- as.double(n)   # Fortran: RN = NT = N
  xm  <- mean(x)
  z   <- x - xm         # Fortran: F(I) = F(I) - XM (demean in-place)

  # Initialise arrays  (Fortran: G(60), H, PEF, PER as local/passed arrays)
  g    <- numeric(lg + 1L)   # G(NN) reflection coefficients, 1-indexed from G(1)
  h    <- numeric(lg + 1L)
  pef  <- numeric(n)         # forward prediction errors
  per  <- numeric(n)         # backward prediction errors
  aic  <- numeric(lg)
  sve  <- numeric(12L)       # AR coefficients SVE(IJK) = -G(IJK+1)

  # Initial variance (Fortran lines 5824-5833):
  # PHI(1) = sum(z^2)/M;  PM = DM = PHI(1);  AIC(1) = RN*log(PM)+2
  ssum   <- sum(z^2)
  phi1   <- ssum / n
  pm     <- phi1
  dm     <- phi1
  aicm   <- rn * log(pm) + 2.0
  ip     <- 0L             # selected AR order

  for (nn in 2L:(lg + 1L)) {
    N <- nn - 2L          # Fortran N = NN-2

    # Initialise PEF/PER to zero on first iteration (Fortran lines 5839-5841)
    if (N == 0L) {
      pef[] <- 0.0
      per[] <- 0.0
    }

    # Compute SN, SD using lattice prediction errors (Fortran lines 5844-5847):
    # JJ = M - N - 1
    jj <- n - N - 1L
    sn <- 0.0; sd <- 0.0
    for (j in seq_len(jj)) {
      ef_j <- z[j + N + 1L] + pef[j]
      er_j <- z[j]           + per[j]
      sn   <- sn - 2.0 * ef_j * er_j
      sd   <- sd + ef_j^2 + er_j^2
    }
    g_nn <- if (abs(sd) > 0.0) sn / sd else 0.0
    g[nn] <- g_nn

    # Update G coefficients (Fortran lines 5853-5858, only if N != 0):
    if (N != 0L) {
      for (j in seq_len(N)) {
        k <- N - j + 2L
        h[j + 1L] <- g[j + 1L] + g_nn * g[k]
      }
      for (j in seq_len(N)) g[j + 1L] <- h[j + 1L]
      jj <- jj - 1L   # Fortran: JJ=JJ-1 before DO 10 (only when N!=0)
    }

    # Update PEF/PER lattice errors (Fortran lines 5860-5862, DO 10 J=1,JJ):
    # Note: PER is updated first, then PEF uses the UPDATED PER(J+1).
    # Loop runs from J=1 to JJ (already decremented if N>0).
    for (j in seq_len(jj)) {
      per_j_new <- per[j] + g_nn * pef[j] + g_nn * z[j + nn - 1L]
      pef[j]    <- pef[j + 1L] + g_nn * per[j + 1L] + g_nn * z[j + 1L]
      per[j]    <- per_j_new
    }

    # Update PHI (Fortran lines 5863-5866)
    phi_sum <- 0.0
    for (j in 2L:nn) phi_sum <- phi_sum - g[nn + 1L - j] * g[j]
    # Fortran: DO 14 J=2,NN: SUM=SUM-PHI(NN+1-J)*G(J)
    # PHI here is a local scratch array; for our purposes only G matters for SVE.
    # (PHI in Fortran is used for the final PHI(40)=PM and the autocorrelation
    #  extension in lines 5881-5885 which we don't need for residuals.)

    dm   <- (1.0 - g_nn^2) * dm
    pm   <- dm
    aic_nn <- rn * log(max(pm, 1e-30)) + 2.0 * nn
    if (nn <= lg) aic[nn - 1L] <- aic_nn

    # IOPT=1 stopping: first minimum AIC (Fortran lines 5871-5876)
    if (aic_nn > aicm && nn > 2L) break     # GOTO 98: stop, ip stays at prev value

    aicm <- aic_nn
    ip   <- min(nn - 1L, 10L)
    for (ijk in seq_len(ip)) sve[ijk] <- -g[ijk + 1L]
  }

  # Compute residuals for selected order ip (Fortran lines 5886-5901):
  # RES(I) = F(I) for I=1..IP (first IP values = demeaned input)
  # RES(I) = F(I) - sum_{j=1}^{IP} SVE(j)*F(I-j) for I=IP+1..M
  if (ip == 0L) {
    resid <- z
  } else {
    resid <- numeric(n)
    resid[1L:ip] <- z[1L:ip]
    for (i in (ip + 1L):n) {
      tmp <- z[i]
      for (j in seq_len(ip)) tmp <- tmp - sve[j] * z[i - j]
      resid[i] <- tmp
    }
  }

  list(
    residuals = resid,   # demeaned AR residuals, matching Fortran WK(I,1) = RES(I)
    order     = ip,
    aic       = aic[seq_len(max(ip, 1L))],
    phi       = -sve[seq_len(max(ip, 1L))]
  )
}


# -----------------------------------------------------------------------------
# .cof_crit99  —  critical correlation at 99% confidence  [CRIT99]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Returns the critical Pearson correlation coefficient at the 99% confidence
#   level for a given segment length, using Holmes' hardcoded lookup table
#   from CRIT99.  If the segment length falls between table entries the next
#   lower entry is used (conservative).
#
# ARGUMENTS
#   lsg   Integer. Segment length in years.
#
# RETURNS
#   Numeric scalar: critical correlation coefficient.
#
.cof_crit99 <- function(lsg) {
  ndf  <- c(10L, 15L, 20L, 25L, 30L, 35L, 40L, 50L, 60L,
             70L, 80L, 90L, 100L, 120L)
  cr99 <- c(0.7155, 0.5923, 0.5155, 0.4622, 0.4226, 0.3916,
             0.3665, 0.3281, 0.2997, 0.2776, 0.2597, 0.2449,
             0.2324, 0.2122)
  crt <- cr99[1L]
  for (i in seq_along(ndf)) {
    if (lsg >= ndf[i]) crt <- cr99[i] else break
  }
  crt
}


# -----------------------------------------------------------------------------
# .cof_qseg  —  find starting years of segments  [QSEG]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Computes the starting years of all segments of a series, given the series
#   span, segment length and lag.  Handles the irregular first and last
#   segments so that the full span of the series is covered without duplication.
#   Translates Holmes' QSEG subroutine exactly.
#
# ARGUMENTS
#   jyr   Integer. First year of the series.
#   lyr   Integer. Last year of the series.
#   ls    Integer. Segment length in years.
#   lag   Integer. Lag between successive segment starts.
#
# RETURNS
#   Integer vector of segment starting years.
#
.cof_qseg <- function(jyr, lyr, ls, lag) {
  n <- lyr - jyr + 1L
  if (n <= ls) return(jyr)

  iyr   <- (jyr %/% lag) * lag
  if (iyr > jyr) iyr <- iyr - lag
  jl    <- iyr - lag
  starts <- integer(0)

  j <- iyr
  while (j <= lyr) {
    ja <- max(j, jyr)
    ja <- min(ja, lyr - ls + 1L)
    if (ja != jl) {
      starts <- c(starts, ja)
      jl <- ja
    }
    j <- j + lag
  }
  starts
}


# -----------------------------------------------------------------------------
# .cof_correl  —  Pearson correlation between two vectors  [CORREL]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Computes the Pearson correlation coefficient between two numeric vectors
#   using the computational formula from Holmes' CORREL subroutine.
#   Returns -9.0 if n <= 1 or the denominator is zero (matching DPL behaviour).
#
# ARGUMENTS
#   x   Numeric vector.
#   y   Numeric vector of the same length as x.
#
# RETURNS
#   Numeric scalar: correlation coefficient in [-1, 1], or -9.0 on error.
#
.cof_correl <- function(x, y) {
  n <- length(x)
  if (n <= 1L) return(-9.0)
  an    <- 1.0 / n
  sumx  <- sum(x);  sumy  <- sum(y)
  sumx2 <- sum(x^2); sumy2 <- sum(y^2)
  sumxy <- sum(x * y)
  denom <- sqrt(abs(sumx2 - an * sumx^2)) * sqrt(abs(sumy2 - an * sumy^2))
  if (abs(denom) < 1e-8) return(0.0)
  (sumxy - an * sumx * sumy) / denom
}


# -----------------------------------------------------------------------------
# .cof_corrxn  —  leave-one-out correlation influence  [CORRXN]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   For a pair of series (z1, z2) of length n, computes the overall Pearson
#   correlation R and then the influence of each point on R: the change in
#   correlation when that point is omitted (R - R_minus_i).  Points are
#   ranked by influence (most lowering first, most raising last), matching
#   the output of Holmes' CORRXN subroutine.
#
# ARGUMENTS
#   z1   Numeric vector, master-minus-self segment.
#   z2   Numeric vector of same length, test series segment.
#   jyr  Integer. First calendar year of the segment (for labelling).
#
# RETURNS
#   List with:
#     $r        Numeric. Overall correlation.
#     $years    Integer vector (length n), calendar years.
#     $delta_r  Numeric vector (length n), influence of each year on r.
#     $ranked_years  Integer vector, years ranked most-lowering to most-raising.
#     $ranked_delta  Numeric vector, delta_r in the same order.
#
.cof_corrxn <- function(z1, z2, jyr) {
  # z1 = master-minus-self (YMSMA), z2 = test series (ZSERM)
  # Returns ranked delta_r with z_sign: '>' if z2 > z1 (series > master), '<' otherwise
  n     <- length(z1)
  years <- jyr + seq(0L, n - 1L)
  if (n < 8L) {
    return(list(r = NA_real_, years = years, delta_r = rep(NA_real_, n),
                ranked_years = years, ranked_delta = rep(NA_real_, n),
                z_sign = rep(">", n)))
  }
  sumx  <- sum(z1);  sumy  <- sum(z2)
  sumx2 <- sum(z1^2); sumy2 <- sum(z2^2)
  sumxy <- sum(z1 * z2)
  an    <- 1.0 / n
  denom <- sqrt(abs(sumx2 - an * sumx^2)) * sqrt(abs(sumy2 - an * sumy^2))
  R     <- if (abs(denom) < 1e-8) 0.0 else (sumxy - an * sumx * sumy) / denom

  # Leave-one-out
  ajab    <- 1.0 / (n - 1L)
  delta_r <- numeric(n)
  for (i in seq_len(n)) {
    xa  <- sumx  - z1[i];  ya  <- sumy  - z2[i]
    x2a <- sumx2 - z1[i]^2; y2a <- sumy2 - z2[i]^2
    xya <- sumxy - z1[i] * z2[i]
    d2  <- sqrt(abs(x2a - ajab * xa^2)) * sqrt(abs(y2a - ajab * ya^2))
    ri  <- if (abs(d2) < 1e-8) 0.0 else (xya - ajab * xa * ya) / d2
    delta_r[i] <- R - ri
  }

  # Sign: '>' if series (z2) > master (z1), '<' otherwise — matches Fortran CORRXN output
  z_sign_raw <- ifelse(z2 > z1, ">", "<")

  ord          <- order(delta_r)
  ranked_years <- years[ord]
  ranked_delta <- delta_r[ord]
  z_sign       <- z_sign_raw[ord]   # reordered to match ranked_years

  list(r = R, years = years, delta_r = delta_r,
       ranked_years = ranked_years, ranked_delta = ranked_delta,
       z_sign = z_sign)
}


# -----------------------------------------------------------------------------
# .cof_cofdif  —  year-to-year first-difference divergence  [COFDIF]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Computes normalised first-differences of the master-minus-self and the
#   test series, then identifies years where the two series diverge by more
#   than a threshold number of standard deviations.  Translates COFDIF.
#
# ARGUMENTS
#   ymsa   Numeric vector. Master-minus-self series, dated (subscripted by year).
#   zser   Numeric vector of same length. Test series.
#   jyr    Integer. First calendar year of the overlap segment.
#   lyr    Integer. Last calendar year of the overlap segment.
#   anot   Numeric. Divergence threshold in SDs (Holmes default: 4.0).
#
# RETURNS
#   Data.frame with columns:
#     year      Integer. Calendar year of divergence (year-1 to year pair).
#     delta     Numeric. Divergence in normalised SD units.
#   Returns an empty data.frame if no divergences exceed the threshold.
#
.cof_cofdif <- function(ymsa, zser, jyr, lyr, anot = 4.0) {
  idx  <- seq_along(ymsa)  # 1-based index
  n    <- length(ymsa)

  # Normalised first-differences of master and test series
  dm <- diff(ymsa)
  dt <- c(0, diff(zser))  # first element is zero (not used)

  nr_m <- .cof_normts(c(0, dm), k = 0L)$z
  nr_t <- .cof_normts(dt,        k = 0L)$z

  divg <- nr_t - nr_m
  hits <- which(abs(divg[-1L]) >= anot) + 1L  # skip first

  if (length(hits) == 0L)
    return(data.frame(year = integer(0), delta = numeric(0)))

  data.frame(
    year  = as.integer(jyr + hits - 1L),
    delta = divg[hits]
  )
}


# -----------------------------------------------------------------------------
# .cof_outabs  —  absent rings and statistical outliers  [OUTABS]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Scans a filtered test series for values that deviate from the master by
#   more than outp SDs (high) or outn SDs (low), and separately records
#   any years at which the raw measurement was zero (absent ring).
#   Translates Holmes' OUTABS subroutine.
#
# ARGUMENTS
#   zser   Numeric vector. Filtered, normalised test series (subscripted jyr:lyr).
#   ymsa   Numeric vector. Filtered master-minus-self, same subscript.
#   sdm    Numeric. Mean standard deviation of the master series across years.
#   jyr    Integer. First calendar year.
#   lyr    Integer. Last calendar year.
#   outp   Numeric. Upper outlier threshold in SDs (Holmes default:  3.0).
#   outn   Numeric. Lower outlier threshold in SDs (Holmes default: -4.5).
#
# RETURNS
#   Data.frame with columns:
#     year   Integer. Calendar year.
#     zsd    Numeric. Departure in SDs from master.
#   Only years exceeding either threshold are included.
#
.cof_outabs <- function(zser, ymsa, sdm, jyr, lyr, outp = 3.0, outn = -4.5) {
  years <- jyr:lyr
  hits  <- data.frame(year = integer(0), zsd = numeric(0))
  if (sdm == 0) return(hits)
  for (k in seq_along(years)) {
    zsd <- (zser[k] - ymsa[k]) / sdm
    if (zsd > outp || zsd < outn)
      hits <- rbind(hits, data.frame(year = years[k], zsd = zsd))
  }
  hits
}


# -----------------------------------------------------------------------------
# .cof_slsg  —  segment sliding-window correlation  [SLSG]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Slides a segment of the test series from -10 to +10 years relative to
#   its dated position and computes the Pearson correlation with the master-
#   minus-self at each lag.  Returns the dated correlation, the position of
#   the maximum, a flag indicating whether the maximum is at a non-zero lag
#   (flag "B") or the dated correlation is below the critical level (flag "A"),
#   and the full 21-point correlation vector for the segment table.
#   Translates Holmes' SLSG subroutine.
#
# ARGUMENTS
#   master_seg   Numeric vector. Master-minus-self over its full span.
#   test_seg     Numeric vector. Test series segment (already extracted,
#                length = lsg).
#   jyrs         Integer. Dated first year of the test segment.
#   jyrm         Integer. First year of the master series.
#   lyrm         Integer. Last year of the master series.
#   lsg          Integer. Segment length.
#   crt          Numeric. Critical correlation threshold.
#
# RETURNS
#   List with:
#     $r_dated   Numeric. Correlation at dated position (lag 0).
#     $r_max     Numeric. Maximum correlation across all lags.
#     $lag_max   Integer. Lag (years) at which r_max occurs; 0 = dated.
#     $flag      Character. "" = OK; "A" = below threshold; "B" = max elsewhere.
#     $corr21    Numeric vector of length 21, correlations at lags -10:+10.
#                Positions outside the master range are -9.99.
#
.cof_slsg <- function(master_seg, test_seg, jyrs, jyrm, lyrm, lsg, crt) {
  corr21 <- numeric(21L)
  nm     <- lyrm - jyrm + 1L

  for (k in 1L:21L) {
    lag <- k - 11L           # -10 to +10
    ia  <- jyrs - jyrm + lag + 1L          # 1-based start in master
    iz  <- lyrm - (jyrs + lsg - 1L) - lag + 1L  # clearance at end
    if (ia < 1L || iz < 1L || ia + lsg - 1L > nm) {
      corr21[k] <- -9.99
    } else {
      corr21[k] <- .cof_correl(master_seg[ia:(ia + lsg - 1L)], test_seg)
    }
  }

  r_dated <- corr21[11L]
  mxk     <- which.max(corr21)
  r_max   <- corr21[mxk]
  lag_max <- mxk - 11L

  flag <- ""
  if (lag_max != 0L) {
    flag <- "B"
  } else if (r_dated < crt) {
    flag <- "A"
  }

  list(r_dated = r_dated, r_max = r_max, lag_max = lag_max,
       flag = flag, corr21 = corr21)
}


# =============================================================================
# COF output formatters  (character-line generators for Parts 1-7)
# =============================================================================


# -----------------------------------------------------------------------------
# .cof_fmt_hisst  —  one series line for the time-span histogram  [HISST]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Generates a single character line in the style of Holmes' HISST subroutine
#   — a 103-character span bar from decade 1000 to 2000 using '<', '=', '>'
#   symbols, followed by the series ident, sequence number, first year, last
#   year and length.
#
# ARGUMENTS
#   id     Character. Series identifier (up to 8 chars).
#   iseq   Integer. Sequential number of this series in the file.
#   jyr    Integer. First year of the series.
#   n      Integer. Length of the series.
#
# RETURNS
#   Character scalar: one formatted line.
#
.cof_fmt_hisst <- function(id, iseq, jyr, n) {
  # Translation of Fortran HISST — formula calibrated against benchmark output.
  #
  # Fortran declares CHARACTER LN(98:200) (103 elements, indexed 98..200).
  # Background: ' ' everywhere; dots at j = 100, 105, ..., 200 (every 5th).
  # Bar: IDA to IDZ filled with '=', '<' at IDA, '>' at IDZ.
  # WRITE(IU,'(T2,103A,T106,A,I4,3I5)') (LN(J),J=98,200)
  #   → bar occupies output columns 2..104 (103 chars); label at col 106.
  #
  # NOTE: The benchmark (generated by a compiled DPL binary) uses the formula
  #   IDA = max((JYR - 50) %/% 10, 98)
  #   IDZ = min((LYR - 50) %/% 10, 200)
  # This differs from the Fortran source formula JYR%/%10 by -5 decades (-50yr),
  # which is the observed empirical offset confirmed across all 12 PUE series.
  # We replicate this empirical behaviour for benchmark-exact output.

  lyr <- jyr + n - 1L

  # Build LN(98:200) — R index i (1-based) = j - 97, i.e. ln[j-97] = LN(j)
  ln <- rep(" ", 103L)
  for (j in seq(100L, 200L, by = 5L)) ln[j - 97L] <- "."

  # Empirical bar extent formula (matches compiled DPL benchmark)
  IDA <- max((jyr - 50L) %/% 10L, 98L)
  IDA <- min(IDA, 200L)
  IDZ <- min((lyr - 50L) %/% 10L, 200L)
  IDZ <- max(IDZ, 98L)

  for (j in IDA:IDZ) ln[j - 97L] <- "="
  # Guard conditions match Fortran: IF(JYR/10 .GE. 98) and IF(LYR/10 .LE. 200)
  if ((jyr - 50L) %/% 10L >= 98L) ln[IDA - 97L] <- "<"
  if ((lyr - 50L) %/% 10L <= 200L) ln[IDZ - 97L] <- ">"

  bar <- paste(ln, collapse = "")
  # Leading space = Fortran CC col 1; bar at cols 2..104; gap; label at col 106
  sprintf(" %s %-8s %3d %4d %4d %4d", bar, id, iseq, jyr, lyr, n)
}


# -----------------------------------------------------------------------------
# .cof_fmt_part2  —  Part 2 time-span histogram  [HISST block]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Generates all lines for Part 2 (the time-span histogram) from a list of
#   series metadata.  Produces header, per-series bars, and scale footer,
#   matching the DPL layout.
#
# ARGUMENTS
#   series_meta   List of lists, each with: $id, $jyr, $n, $seq.
#   title         Character. Run title string for page heading.
#
# RETURNS
#   Character vector of formatted lines (no trailing newlines).
#
.cof_fmt_part2 <- function(series_meta, title) {
  # Exact translation of Fortran HISST header/footer.
  #
  # Header WRITE format (Fortran):
  #   WRITE(IU,'(I5,10I10,T119,"Beg   End"/T106,
  #   &"Ident    Seq year year  Yrs"/"   :",20("    :"),T106,
  #   &"-------- --- ---- ---- ----")')(J,J=1000,2000,100)
  #
  # Scale line: I5(1000) + 10×I10(1100..2000) = 105 chars total.
  # "Beg   End" at T119 (col 119, 1-indexed → 14 spaces padding after 105 chars).
  # Col-2 label at T106 = 105 leading spaces.
  # Tick row: '   :' + 20×'    :' (104 chars).
  # Separator '-------- --- ---- ---- ----' at T106 = col 106.
  #
  # Footer WRITE format:
  #   WRITE(IU,'(/"   :",20("    :")/I5,10I10)')(J,J=1000,2000,100)
  # → blank line, tick row, scale row (no extra label columns).

  # Scale row: I5 for 1000, then I10 for 1100..2000
  scale_row <- paste0(
    formatC(1000L, width = 5L),
    paste(formatC(seq(1100L, 2000L, by = 100L), width = 10L), collapse = "")
  )  # 105 chars

  # Scale header line: scale + "Beg   End" at col 119 (14 chars gap)
  scale_row_hdr <- sprintf("%-118s%s", scale_row, "Beg   End")

  # Tick row: '   :' + 20×'    :' = 104 chars
  tick_row <- paste0("   :", paste(rep("    :", 20L), collapse = ""))

  # Labels at T106 (col 106 = 105 leading chars from start of line)
  col106 <- strrep(" ", 105L)

  # Separator line: tick_row is 104 chars, separator starts at col 106
  # → pad tick_row to 105 chars, then append separator
  sep_line <- sprintf("%-105s%s", tick_row, "-------- --- ---- ---- ----")

  hdr <- c(
    sprintf("PART 2:  TIME PLOT OF TREE-RING SERIES: %s", title),
    strrep("-", 132L),
    scale_row_hdr,
    sprintf("%s%s", col106, "Ident    Seq year year  Yrs"),
    sep_line
  )

  bars <- vapply(series_meta, function(s)
    .cof_fmt_hisst(s$id, s$seq, s$jyr, s$n), character(1))

  # Footer: blank line, tick row, scale row
  foot <- c("", tick_row, scale_row)

  c(hdr, bars, foot)
}


# -----------------------------------------------------------------------------
# .cof_fmt_part3  —  Part 3 master dating series listing  [CRONV]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Generates the vertical listing of the master dating series (Part 3),
#   showing value, sample depth and absent-ring count for each year, formatted
#   in six columns of 50 years each — matching Holmes' CRONV subroutine output.
#
# ARGUMENTS
#   master   Named numeric vector (names = years as character).
#   depth    Named integer vector, sample depth per year (same names as master).
#   absent   Named integer vector, absent rings per year (same names as master).
#   title    Character. Run title string.
#
# RETURNS
#   Character vector of formatted lines.
#
.cof_fmt_part3 <- function(master, depth, absent, title) {
  # Exact translation of Fortran CRONV subroutine (lines 2097-2144).
  #
  # ESPA(J) is a 22-char field: WRITE(ESPA(J),'(I6,F7.3)') N, Y(N)
  #   cols  1- 6: year as I6 (right-justified)
  #   cols  7-13: value as F7.3 (leading zero suppressed: '  .116' not '  0.116')
  #   cols 14-17: NC(N) sample depth as I4
  #   cols 18-20: NA(N) absent rings as I3
  #   cols 21-22: '<<' if absent ring AND value > -0.4
  # Row output: WRITE(IO,'(6A)') ESPA → 6×22 = 132 chars
  # Blank line: IF(MOD(I,10).EQ.9 .AND. I.LT.IA+49) → blank every 10 rows, not on last
  #
  # IA formula (Fortran lines 2106-2110):
  #   IA = ((JYR+100)/300)*300 - 100   (integer division)
  #   WHILE IA > JYR: IA = IA - 100    (back up in 100-year steps until IA <= JYR)
  # Page covers 300 years (IA..IA+299), 6 columns × 50 rows.
  # Next page: IA = IA + 300 (line 2141).

  years <- as.integer(names(master))
  jyr   <- min(years);  lyr <- max(years)

  hdr <- c(
    sprintf("PART 3:  Master Dating Series: %s", title),
    strrep("-", 132L),
    paste(rep("  Year  Value  No Ab  ", 6L), collapse = ""),
    paste(rep("  ------------------  ", 6L), collapse = "")
  )

  # Correct IA formula from Fortran CRONV lines 2106-2110
  ia <- ((jyr + 100L) %/% 300L) * 300L - 100L
  while (ia > jyr) ia <- ia - 100L

  out_lines <- character(0)
  repeat {
    for (i in ia:(ia + 49L)) {
      espa <- vapply(1L:6L, function(j) {
        yr <- i + (j - 1L) * 50L
        if (yr < jyr || yr > lyr) return(strrep(" ", 22L))
        v    <- master[as.character(yr)]
        nd   <- depth[as.character(yr)]
        na_  <- absent[as.character(yr)]
        na_i <- if (is.na(na_)) 0L else as.integer(na_)
        nd_i <- if (is.na(nd))  0L else as.integer(nd)
        # Build 22-char ESPA cell: I6(6) + F7.3(7) + 9 spaces = 22 total
        # Fortran: CHARACTER ESPA(6)*22 — always exactly 22 chars per cell
        cell <- paste0(formatC(yr, width = 6L),
                       if (is.na(v)) "       " else .cof_f73(v),
                       strrep(" ", 9L))
        cell <- substring(cell, 1L, 22L)
        substr(cell, 14L, 17L) <- formatC(nd_i, width = 4L)
        if (na_i > 0L) {
          substr(cell, 18L, 20L) <- formatC(na_i, width = 3L)
          if (!is.na(v) && v > -0.4) substr(cell, 21L, 22L) <- "<<"
        }
        cell
      }, character(1))
      out_lines <- c(out_lines, paste(espa, collapse = ""))
      # Blank line every 10 rows within the block, but NOT on the last row
      if (i %% 10L == 9L && i < ia + 49L)
        out_lines <- c(out_lines, " ")
    }
    ia <- ia + 300L        # page increment: 300 years (Fortran line 2141)
    if (ia > lyr) break
    out_lines <- c(out_lines, "")
  }
  c(hdr, out_lines, strrep("-", 132L))
}


# -----------------------------------------------------------------------------
# .cof_fmt_part4  —  Part 4 ASCII bar plot of master series  [BARPL]
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Generates the ASCII bar plot of the (normalised) master dating series,
#   replicating Holmes' BARPL subroutine output exactly.  Each value is
#   represented by a horizontal bar whose length corresponds to its decile
#   rank; the terminal symbol encodes the standard-deviation class (A-Z for
#   positive, a-z for negative, @ for zero).
#
#   Verified line-for-line against the DPL COFECHA 4.04P PEL benchmark
#   (PELCOF.OUT, 161 years, 2 pages) and the PUE benchmark (67 bars).
#
# ALGORITHM  (direct translation of Fortran BARPL; shared helpers)
#   1. NORMTS(NY,Y,XM,SD,0) — normalise to mean=0, population SD (k=0).
#   2. .cof_barpl_cuts: RANKRI descending → 10 decile cut-points
#        Z(j) = 0.5*(sorted[J1]+sorted[J2])
#        J1 = NINT(NY/11*(j-0.5)),  J2 = NINT(NY/11*(j+0.5))
#      NINT = round half away from zero (.cof_nint), NOT R's round().
#   3. .cof_barpl_car: LB = 16; for k=1..10: if Y<Z(k) then LB=16-k;
#        LB=max(6,min(16,LB)).  LP = NINT(Y*4): negative → chr(96-LP)
#        capped '<'; positive → chr(64+LP) capped '>'.
#        CAR cell (16 chars): I5 year + (LB-6) dashes + symbol.
#   4. .cof_barpl_pages: IA=(JYR/400)*400; pages of 400 yr, 8 cols × 50 rows;
#        column header, rows, blank decade separator after rows 9/19/29/39
#        (DPL prints ' ----' there), blank line closing the page; page
#        header repeated per page.
#
# ARGUMENTS
#   master   Named numeric vector. The master dating series (named by year).
#            Must already be normalised (output of .cof_normts passed as
#            master_norm from dpl_cof) OR raw — normalisation is applied
#            internally, matching Fortran's second NORMTS call inside BARPL.
#   title    Character. Run title string.
#
# RETURNS
#   Character vector of formatted lines.
#
.cof_fmt_part4 <- function(master, title) {
  years <- as.integer(names(master))
  jyr   <- min(years);  lyr <- max(years)

  # NORMTS(NY, Y, XM, SD, k=0) -- population SD, as in Fortran BARPL
  Y <- setNames(.cof_normts(as.numeric(master), k = 0L)$z, names(master))
  Z <- .cof_barpl_cuts(Y)

  make_car <- function(yr) {
    if (yr < jyr || yr > lyr) return(strrep(" ", 16L))
    .cof_barpl_car(yr, Y[[as.character(yr)]], Z)
  }

  page_hdr <- c(sprintf("PART 4:  Master Bar Plot: %s", title), strrep("-", 132L))
  .cof_barpl_pages(make_car, jyr, lyr, page_hdr)
}


# -----------------------------------------------------------------------------
# .cof_fmt_seg_table  —  Part 5 segment correlation table for one series
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Generates the Part 5 segment correlation table lines for a single series,
#   in the style of Holmes' SGTBL subroutine.  Called once per series during
#   the second loop; the results are collected into the full Part 5 block by
#   the main engine.
#
# ARGUMENTS
#   id         Character. Series identifier.
#   seq_no     Integer. Sequential number of this series.
#   jyr        Integer. First year of the series.
#   lyr        Integer. Last year of the series.
#   n          Integer. Series length.
#   segs       List of segment results from .cof_slsg(), each with:
#                $start   Integer. Segment start year.
#                $r_dated Numeric. Correlation at dated position.
#                $flag    Character. "", "A", or "B".
#   iam        Integer. First year of master series.
#   izm        Integer. Last year of master series.
#   lag        Integer. Lag between successive segments.
#   crt        Numeric. Critical correlation threshold.
#
# RETURNS
#   List with:
#     $lines     Character vector of formatted text lines.
#     $n_segs    Integer. Total number of segments checked.
#     $n_below   Integer. Segments with flag "A".
#     $n_other   Integer. Segments with flag "B".
#
.cof_fmt_seg_table <- function(id, seq_no, jyr, lyr, n, segs,
                               iam, izm, ia2, iz2, lag, crt,
                               jyr_clipped) {
  # Exact translation of Fortran SGTBL (lines 1319-1400).
  #
  # SGTBL receives IA2/IZ2 as its IAM/IZM (Fortran line 822):
  #   CALL SGTBL(..., IA2, IZ2, ...)
  # so the display grid spans only the portion with ≥2 series.
  #
  # ISEGA formula (lines 1320-1321):
  #   ISR = IAX = ((JYR_clipped + LAG - 1) / LAG) * LAG
  #   ISEGA = (ISR - IA2) / LAG + 2
  #   IF (JYR_orig < ISR) ISEGA = ISEGA - 1
  #
  # Column mapping (lines 1371-1372):
  #   KSEG = KSEG + 1  (increments 1..NSP)
  #   KR   = KSEG - ISEGA + 1  → KR=1 when KSEG=ISEGA
  # So segment k (0-indexed) maps to master-grid column ISEGA + k (1-indexed).
  #
  # Overprint merge ('+' carriage-control, Fortran lines 1395-1396):
  #   Data at T25; underline '+' row at T26 (+1 offset).
  #   UND(j) = ' ___A' or ' ___B'; 5th char ('A'/'B') lands at start of next cell.
  #   Merge: skip '+' at index 1, for ui >= 2: data[ui] = under[ui] if non-sp/non-und.
  #
  # fmt_r — Fortran F5.2 leading-zero suppression (lines 1387-1390):
  #   ' 0.51' → '  .51'  (5 chars; char2 moves to char2's slot, char1 becomes space)
  #   '-0.16' → ' -.16'  (5 chars)

  n_segs  <- length(segs)
  n_below <- sum(vapply(segs, function(s) s$flag == "A", logical(1)))
  n_other <- sum(vapply(segs, function(s) s$flag == "B", logical(1)))

  # Master grid using IA2/IZ2
  ifa <- ((ia2 - 1L) %/% lag) * lag
  nsp <- (iz2 - ia2) %/% lag + 1L

  # ISEGA: column in master grid where this series' first segment sits
  iax   <- ((jyr_clipped + lag - 1L) %/% lag) * lag
  isega <- (iax - ia2) %/% lag + 2L
  if (jyr < iax) isega <- isega - 1L

  # Value format: F4.2 (4 chars) + 1 trailing space = 5-char cell.
  # Position 5 is always ' ', so the flag letter at UND[[col]][5] replaces it → ' .51A'.
  # F4.2 with leading-zero suppression: ' 0.51' (F5.2) → drop leading '0' → ' .51' (4 chars).
  fmt_r5 <- function(r) {
    if (is.na(r) || r <= -9.0) return("   - ")
    if (r < 0.0) {
      s <- sprintf("%.2f", -r)                    # "0.16"
      paste0("-", substr(s, 2L, 4L), " ")         # "-.16 " (5 chars)
    } else {
      s <- sprintf("%.2f", r)                     # "0.51"
      paste0(" ", substr(s, 2L, 4L), " ")         # " .51 " (5 chars)
    }
  }

  # Build RVL (value) and UND (flag underline) arrays — NSP elements each
  b   <- "   - "
  rvl <- rep(b,       nsp)
  und <- rep("     ", nsp)

  for (k in seq_along(segs)) {
    col <- isega + k      # 1-indexed master-grid column for this segment
    if (col < 1L || col > nsp) next
    r   <- segs[[k]]$r_dated
    flg <- segs[[k]]$flag
    if (!is.na(r) && r > -9.0) rvl[col] <- fmt_r5(r)
    # Flag letter at position 5 of UND[[col]]: replaces the trailing space of ' .51 ' → ' .51A'.
    if (flg %in% c("A", "B")) {
      und[[col]] <- paste0("    ", flg)   # 4 spaces + flag at pos 5
    }
  }

  # Data row: T25 → values start at col 25 of the line.
  # Under row: '+' at col 1, then 24 spaces to reach col 26, then UND values.
  # UND '    A' places 'A' at position 5 of the cell = the trailing space of ' .51 '.
  data_row  <- sprintf("%4d %-8s%5d%5d  %s",
                       seq_no, id, jyr, lyr, paste(rvl, collapse = ""))
  under_row <- paste0("+", strrep(" ", 24L), paste(und, collapse = ""))

  # Overprint merge: skip '+' at pos 1. Replace data char with UND char where UND != ' '.
  da <- strsplit(data_row,  "")[[1L]]
  ua <- strsplit(under_row, "")[[1L]]
  for (ui in seq(2L, length(ua))) {
    if (ua[ui] != " " && ui <= length(da)) da[ui] <- ua[ui]
  }
  merged_row <- paste(da, collapse = "")

  list(
    lines      = merged_row,
    n_segs     = n_segs,
    n_below    = n_below,
    n_other    = n_other,
    isega      = isega,
    nsp        = nsp,
    ifa        = ifa,
    seg_starts = vapply(segs, function(s) as.integer(s$start), integer(1))
  )
}


# -----------------------------------------------------------------------------
# .cof_fmt_part6_series  —  Part 6 potential-problems block for one series
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Generates the Part 6 diagnostic output block for a single series,
#   combining [A] segment lag tables, [B] influence values, [C] first-
#   difference divergences, [D] absent rings and [E] outliers into a single
#   character vector of formatted lines.  This is the function called by
#   both dpl_cof() (for the full output file) and dpl_cof_diag() (for
#   selective on-screen printing).
#
# ARGUMENTS
#   id          Character. Series identifier.
#   seq_no      Integer. Sequential number.
#   jyr         Integer. First year of the series.
#   lyr         Integer. Last year of the series.
#   n           Integer. Series length.
#   problems    List with elements from one series in cof$problems:
#                 $low_r    — list of flagged segments with $start, $r_dated,
#                             $r_max, $lag_max, $corr21, $flag
#                 $corrxn   — result of .cof_corrxn() for each flagged seg
#                 $cofdif   — data.frame from .cof_cofdif()
#                 $absent   — integer vector of absent ring years
#                 $outliers — data.frame from .cof_outabs()
#   crt         Numeric. Critical correlation threshold.
#   lsg         Integer. Segment length.
#
# RETURNS
#   Character vector of formatted lines for this series' Part 6 block.
#
.cof_fmt_part6_series <- function(id, seq_no, jyr, lyr, n,
                                  problems, crt, lsg, lag,
                                  master_ysd = NULL, depth = NULL,
                                  absent_by_yr = NULL) {
  # Exact translation of Fortran Part 6 output routines.
  #
  # Series separator + header (Fortran lines 689-690):
  #   WRITE(IU(26),'(1X,131(''='')//1X,A,I6,'' to'',I6,I8,'' years'',T123,''Series'',I4)')
  #   → ' ' + 131×'=' + blank line + ' ' + ID(A) + I6 jyr + ' to' + I6 lyr + I8 n + T123 + 'Series' + I4 seq
  #
  # [A] segment tables: the WRITE in Fortran SLSG for [A] header is COMMENTED OUT
  #   (lines 1547-1549 have 'c' prefix). [A] is NEVER emitted in standard COFECHA.
  #
  # [B] Fortran CORRXN:
  #   WRITE(IU,'(/'' [B] Entire series, effect on correlation ('',F6.3,'') is:'')')R
  #   WRITE(IU,'(7X,''Lower'',4(I7,F7.3),''  Higher'',4(I7,F7.3))')
  #       (JAHR(J),RC(J),J=N,N-3,-1),(JAHR(I),RC(I),I=1,4)
  #   → RANKRI sorts RC ascending (most-negative delta first).
  #   → J=N,N-3,-1: 4 values from the HIGH end (most raising = least negative delta).
  #   → I=1,4: 4 values from the LOW end (most lowering = most negative delta).
  #   → Labeled 'Lower' (4 values) and 'Higher' (4 values) — but benchmark shows
  #     'Lower' with 6 values and 'Higher' with 2 values and year</>marker.
  #   → The '<' or '>' after year: Z2(I) > Z1(I) means series > master → '>'
  #     Z2(I) < Z1(I) means series < master → '<'. This comes from comparing
  #     ZSERM and YMSMA at that year inside CORRXN (not stored in RC/JAHR directly).
  #     We replicate this using problems$corrxn$z1 and $z2 if available.
  #   → F6.3 for correlation: '  .515' (6 chars, 3 decimals, no leading zero)
  #   → I7 year (7 chars) + F7.3 delta (7 chars, no leading zero, with sign)

  # Series header: ' ' + 131×'=' + blank + series title line
  hdr_sep <- paste0(" ", strrep("=", 131L))
  # Fortran: 1X,A,I6,' to',I6,I8,' years',T123,'Series',I4
  # T123 means tab to column 123. Build by padding.
  title_body <- sprintf(" %-8s%6d to%6d%8d years", id, jyr, lyr, n)
  series_tag <- sprintf("Series%4d", seq_no)
  # Pad to col 123 (1-indexed) = 122 chars before 'Series'
  title_line <- formatC(title_body, width = 122L, flag = "-")
  title_line <- paste0(title_line, series_tag)

  lines <- c(hdr_sep, "", title_line, "")

  # [A] Segment lag-correlation table (Fortran SLSG lines 1543-1577)
  #
  # The [A] header line WAS commented out in the Fortran source (line 1548-1549 prefixed 'c'),
  # but the benchmark output shows it — the version that produced the benchmark had it active.
  # The individual data rows (lines 1572, 1577) ARE emitted for every flagged segment.
  #
  # Format (Fortran line 1572): I8 start, I5 end, 7X, 21×A5 correlations
  # Overprint row (line 1577):  '+', 13X, I4 lag_offset, 3X, 21×A5 markers
  #   After merge: '*' marks the MAXIMUM position, '|' marks the DATED position (lag 0).
  #   When MXCOR=11 (A-flag): '*' at lag 0 only (max = dated, no '|' needed).
  #   When MXCOR≠11 (B-flag): '|' at lag 0, '*' at max position.
  # Leading-zero suppression same as Part 5 (Fortran lines 1566-1569).
  # Separator line between non-consecutive segments (line 1553):
  #   if |JYRS - JYRL| > LAG: emit '   ' + 61×' -'

  # Shared 5-char F5.2 formatter with leading-zero suppression (Fortran lines 1563-1569)
  fmt_a <- function(r) {
    if (is.na(r) || r <= -9.0) return("   - ")
    rv <- sprintf("%5.2f", r)
    if (substr(rv, 2L, 2L) == "0")
      paste0(" ", substr(rv, 1L, 1L), substr(rv, 3L, 5L))  # 5 chars: ' -.16' or '  .51'
    else rv
  }

  if (length(problems$low_r) > 0L) {
    has_hdr  <- FALSE
    prev_start <- NA_integer_

    for (sl in problems$low_r) {
      if (sl$flag == "") next          # skip unflagged (shouldn't happen, but guard)
      c21  <- sl$corr21
      mxk  <- which.max(c21)          # 1-indexed, 1=lag-10 ... 11=lag0 ... 21=lag+10
      lag_off <- mxk - 11L            # lag offset of maximum (-10..+10)

      # [A] header printed once before the first flagged segment of each series
      if (!has_hdr) {
        lines <- c(lines,
          " [A] Segment   High   -10   -9   -8   -7   -6   -5   -4   -3   -2   -1   +0   +1   +2   +3   +4   +5   +6   +7   +8   +9  +10",
          "    ---------  ----   ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---"
        )
        has_hdr <- TRUE
      }

      # Separator between non-consecutive segments (Fortran line 1553: > LAG)
      if (!is.na(prev_start) && abs(sl$start - prev_start) > lag) {
        lines <- c(lines, paste0("   ", strrep(" -", 61L)))
      }
      prev_start <- sl$start

      # Format the 21 correlation values
      ln_vals <- vapply(c21, fmt_a, character(1))

      # Build marker arrays: blanks except at lag-0 ('|') and max position ('*')
      ln_und <- rep("     ", 21L)
      ln_und[11L] <- "    |"             # dated position marker
      ln_und[mxk] <- "    *"             # maximum position marker (overrides if mxk=11)

      # Data row: I8 start, I5 end, 7X, 21×A5  (Fortran line 1572)
      data_a  <- sprintf("%8d%5d       %s", sl$start, sl$start + lsg - 1L,
                         paste(ln_vals, collapse = ""))

      # Overprint row: '+', 13X, I4 lag_offset, 3X, 21×A5  (Fortran line 1577)
      # '+' is carriage-control (not printed). Content at col 14 onward.
      # Merge: skip index 1 ('+'), for ui>=2: data[ui]=under[ui] if non-sp/non-und
      und_a   <- sprintf("+%13s%4d   %s", "", lag_off,
                         paste(ln_und, collapse = ""))

      da <- strsplit(data_a, "")[[1L]]
      ua <- strsplit(und_a,  "")[[1L]]
      for (ui in seq(2L, length(ua))) {
        ch <- ua[ui]
        if (ch != " " && ch != "_" && ui <= length(da)) da[ui] <- ch
      }
      lines <- c(lines, paste(da, collapse = ""))
    }

    if (has_hdr) lines <- c(lines, "")
  }

  # [B] Influence values — emitted for EVERY series (called from SLSG for full overlap)
  cx <- problems$corrxn
  if (!is.null(cx) && !is.na(cx$r)) {
    nr <- length(cx$ranked_years)
    n_show <- min(4L, nr)

    # Low end (most lowering): ranked_years[1..4], ranked_delta[1..4]
    low_y <- cx$ranked_years[1L:n_show]
    low_d <- cx$ranked_delta[1L:n_show]
    # High end (most raising): ranked_years[N..N-3] reversed
    high_y <- rev(cx$ranked_years[(nr - n_show + 1L):nr])
    high_d <- rev(cx$ranked_delta[(nr - n_show + 1L):nr])

    # Sign indicator for each year: '>' if series value > master, '<' otherwise
    sign_for <- function(yrs, delta) {
      vapply(seq_along(yrs), function(k) {
        # Positive delta means this year LOWERED correlation by being removed →
        # the value was pulling r DOWN. Negative delta (raising removal) →
        # the value was pulling r UP. The '<'/'>' shows series vs master value.
        # We use the sign of the delta itself as a proxy (matches benchmark pattern):
        # delta > 0 → removing this year LOWERED r → year was outlier in one direction
        # The actual sign (< or >) compares ZSERM[i] vs YMSMA[i].
        # Since we don't store per-year z1/z2 in the result, use delta sign
        # as a proxy: positive delta → series < master for that year → '<'
        #             negative delta → series > master for that year → '>'
        # This matches the Fortran pattern shown in benchmark.
        if (!is.null(cx$z_sign) && length(cx$z_sign) >= k)
          cx$z_sign[k]
        else if (delta[k] > 0) "<" else ">"
      }, character(1))
    }

    low_sign  <- sign_for(low_y,  low_d)
    high_sign <- sign_for(high_y, high_d)

    # F6.3 correlation with no leading zero — use shared helper
    r_str <- .cof_f63(cx$r)

    # Format each year+sign+delta: I7 year + sign char + F7.3 delta
    fmt_yr_delta <- function(yr, sgn, d) {
      sprintf("%7d%s%s", yr, sgn, .cof_f73(d))
    }

    low_str  <- paste(mapply(fmt_yr_delta, low_y,  low_sign,  low_d),  collapse = "")
    high_str <- paste(mapply(fmt_yr_delta, high_y, high_sign, high_d), collapse = "")

    lines <- c(lines,
      sprintf("\n [B] Entire series, effect on correlation (%s) is:", r_str),
      sprintf("       Lower %s  Higher %s", low_str, high_str)
    )
  }

  # [C] Year-to-year first-difference divergences
  if (!is.null(problems$cofdif) && nrow(problems$cofdif) > 0L) {
    cd <- problems$cofdif
    lines <- c(lines,
      " [C] Year-to-year changes very different from the mean change in other series",
      paste(sprintf("%8d%+6.1f SD;", cd$year, cd$delta), collapse = "  ")
    )
  }

  # [D] Absent rings
  # Fortran (lines 843-854):
  #   WRITE(IU(26),'(/'' [D]'',I5,'' Absent rings:  Year   Master  N series Absent'')')JABB
  #   Each: WRITE(IU(26),'(I29,F9.3,2I8)') JP(K),YSD(JP(K)),JNS(JP(K)),JAB(JP(K))
  #   Warning overprint for wide ring (overprint → not emitted in plain text)
  if (length(problems$absent) > 0L) {
    ab_yrs <- problems$absent
    lines <- c(lines,
      sprintf("\n [D]%5d Absent rings:  Year   Master  N series Absent",
              length(ab_yrs))
    )
    for (ayr in ab_yrs) {
      yr_ch <- as.character(ayr)
      mv  <- if (!is.null(master_ysd) && yr_ch %in% names(master_ysd))
               master_ysd[yr_ch] else NA_real_
      nd  <- if (!is.null(depth) && yr_ch %in% names(depth))
               depth[yr_ch] else NA_integer_
      nab <- if (!is.null(absent_by_yr) && yr_ch %in% names(absent_by_yr))
               absent_by_yr[yr_ch] else NA_integer_
      # I29 year, F9.3 master, I8 N_series, I8 N_absent
      lines <- c(lines,
        sprintf("%29d%9.3f%8d%8d",
                ayr,
                ifelse(is.na(mv), 0.0, mv),
                ifelse(is.na(nd), 0L, as.integer(nd)),
                ifelse(is.na(nab), 0L, as.integer(nab)))
      )
    }
  }

  # [E] Outliers (OUTABS, Fortran lines 1463-1467)
  # Header prints fixed OUTP/OUTN thresholds — NOT the observed max/min ZSD.
  # Fortran: WRITE(IU,'(...'' [E] Outliers'',I6,F6.1,'' SD above or'',F5.1,'' SD below mean for year'')')
  #           NOUT, OUTP, OUTN
  if (!is.null(problems$outliers) && nrow(problems$outliers) > 0L) {
    ot     <- problems$outliers
    outp_v <- if (!is.null(problems$outp)) problems$outp else 3.0
    outn_v <- if (!is.null(problems$outn)) problems$outn else -4.5
    lines <- c(lines,
      sprintf("\n [E] Outliers%6d%5.1f SD above or%5.1f SD below mean for year",
              nrow(ot), outp_v, outn_v)
    )
    chunks <- split(seq_len(nrow(ot)), ceiling(seq_len(nrow(ot)) / 7L))
    for (ch in chunks) {
      lines <- c(lines,
        paste(sprintf("%8d%+5.1f SD;", ot$year[ch], ot$zsd[ch]), collapse = "  ")
      )
    }
  }

  lines
}


# -----------------------------------------------------------------------------
# .cof_fmt_part7_header  —  Part 7 table column header lines
# -----------------------------------------------------------------------------
#
# RETURNS   Character vector (3 lines) for the Part 7 statistics table header.
#
.cof_fmt_part7_header <- function() {
  # Exact Fortran WRITE (lines 810-816):
  # WRITE(IU(27),'(/T49,''Corr   //-------- Unfiltered --------\\  '',
  # &''//---- Filtered -----\\''/T24,3(''    No.''),
  # &T49,''with   Mean   Max     Std   Auto   Mean   Max     Std   Auto  AR''/
  # &'' Seq Series   Interval   Years  Segmt  Flags   Master  msmt   msmt    dev   corr   sens  '',
  # &''value    dev   corr  ()''/
  # &'' --- -------- ---------'',3(''  -----''),''   ------'',8('' ----- ''),'' --'')')
  # T49 = tab to col 49; T24 = tab to col 24
  c(
    "",
    sprintf("%s%s", strrep(" ", 48L),
            "Corr   //-------- Unfiltered --------\\\\  //---- Filtered -----\\\\"),
    sprintf("%s%s%s",
            strrep(" ", 23L), "    No.    No.    No.",
            sprintf("%s%s", strrep(" ", 48L - 23L - 21L),
                    "with   Mean   Max     Std   Auto   Mean   Max     Std   Auto  AR")),
    " Seq Series   Interval   Years  Segmt  Flags   Master  msmt   msmt    dev   corr   sens  value    dev   corr  ()",
    paste0(" --- -------- ---------",
           "  -----  -----  -----",
           "   ------",
           paste(rep(" ----- ", 8L), collapse = ""),
           " --")
  )
}


.cof_fmt_part7_row <- function(seq_no, id, jyr, lyr, n,
                                n_segs, n_flags, r_master,
                                stats_u, stats_f, ar_order) {
  # Exact Fortran format, assembled in two steps:
  # WRITE(LN,'(A,2F7.2,3F7.3,F7.2,2F7.3,I4)')
  #   ID, YMEAN, YMX, SD, ACOR, SEN, YMXF, SDF, ACORF, IAR
  # Then WRITE(QLN4(44:132),'(F8.3,A)') RAV, LN(9:68)
  # QLN4(1:43) from SGTBL: I3 seq, 1X, A8 id, I5 jyr, I5 lyr, I7 n, I7 nseg, I7 nund
  # Final: WRITE(IU(27),'(1X,A)') QLN4
  #
  # F8.3 for r_master: '   .515' (8 chars, no leading zero)
  # F7.2 for mean/max/filt_max: '   1.56' (7 chars, 2 decimals)
  # F7.3 for sd/acor/sen/filt_sd/filt_acor: '   .647' (no leading zero)
  # I4 for AR order
  # All format helpers use the shared .cof_fXX functions defined at module top.

  # QLN4(1:43): I3, 1X, A8, I5, I5, I7, I7, I7
  q1_43 <- sprintf("%3d %-8s%5d%5d%7d%7d%7d",
                   seq_no, id, jyr, lyr, n, n_segs, n_flags)
  # QLN4(44:51): F8.3 r_master
  q44 <- .cof_f83(r_master)
  # LN(9:68): 2F7.2, 3F7.3, F7.2, 2F7.3, I4
  ln9_68 <- paste0(
    .cof_f72(stats_u["mean"]), .cof_f72(stats_u["max"]),
    .cof_f73(stats_u["sd"]),   .cof_f73(stats_u["acor"]),  .cof_f73(stats_u["sen"]),
    .cof_f72(stats_f["max"]),
    .cof_f73(stats_f["sd"]),   .cof_f73(stats_f["acor"]),
    sprintf("%4d", ar_order)
  )
  paste0(" ", q1_43, q44, ln9_68)
}


#' Run COFECHA quality control and crossdating
#'
#' @description
#' Runs the complete COFECHA algorithm (Holmes 1983, 1994) on a set of
#' crossdated ring-width series, implementing the two-loop procedure:
#'
#' **First loop** --- applies a cubic smoothing spline, optional log-transform
#' and Burg AR prewhitening to each series, and accumulates the weighted mean
#' master dating chronology.
#'
#' **Second loop** --- removes each series from the master, segments the series,
#' slides each segment +/-10 years, flags low correlations and dating offsets,
#' identifies influential values, detects divergent year-to-year changes,
#' absent rings, and statistical outliers.
#'
#' All internal routines (spline, AR model, log-transform, normalisation) are
#' translated directly from Holmes' Fortran source to ensure numerical
#' equivalence with DPL reference output.
#'
#' @param rwl A dplR `rwl` data.frame with integer-coercible row names
#'   (calendar years) and at least two series columns.
#' @param seg_length Integer. Segment length in years correlated against the
#'   master. Holmes default: `50`.
#' @param seg_lag Integer. Lag between successive segment starts. Holmes
#'   default: `25` (half of `seg_length`).
#' @param spline_period Integer. 50\%-frequency-response wavelength of the
#'   cubic smoothing spline used to remove low-frequency variance before
#'   cross-correlation. Holmes default: `32`.
#' @param ar_model Logical. Apply Burg autoregressive prewhitening before
#'   computing master correlations. Default `TRUE` (Holmes `QM = 'A'`).
#' @param log_transform Logical. Log-transform the filtered series (additive
#'   constant = 1/3 of mean) before normalisation. Default `TRUE`
#'   (Holmes `QLGT = 'Y'`).
#' @param crit_level Numeric or `NULL`. Critical correlation threshold below
#'   which a segment is flagged. `NULL` (default) uses Holmes' CRIT99 table
#'   at 99\% confidence for `seg_length`.
#' @param outp Numeric. Upper outlier threshold in SDs. Holmes default: `3.0`.
#' @param outn Numeric. Lower outlier threshold (negative SD).
#'   Holmes default: `-4.5`.
#' @param parts Integer vector. Subset of `1:7` controlling which output parts
#'   appear in `$output` and `output_file`. Defaults to all parts. Part 6
#'   diagnostics are always computed and stored in `$problems` regardless.
#' @param output_file Character or `NULL`. If supplied, writes the full
#'   formatted COFECHA output to that file.
#' @param min_length Integer (default `10L`). Series with fewer measured years
#'   are excluded from the master and from segment testing, listed in Part 1
#'   and returned in `$short`. COFECHA's critical-correlation table starts at
#'   10 years; shorter series cannot be tested meaningfully with sliding
#'   segments. Use \code{\link{dpl_short}} to check them against the master.
#' @param verbose Logical. Print per-series progress and summary box. Default
#'   `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`$short`}{`data.frame` (`series`, `jyr`, `lyr`, `n`) of series
#'     excluded as shorter than `min_length` (empty if none).}
#'   \item{`$filtered`}{`rwl` data.frame of the fully transformed series
#'     (spline-detrended, log-transformed, AR-prewhitened, normalised) on the
#'     master's year axis: the values COFECHA correlates in Part 5. Use with
#'     \code{\link{dpl_cormat}} or your own analyses of the year-to-year
#'     signal.}
#'   \item{`$master`}{Named numeric vector. Globally normalised master dating
#'     series (names = character years). Used for Part 3/4 output and plotting.}
#'   \item{`$master_raw`}{Named numeric vector. Pre-normalisation mean of
#'     z-norm values. DPL-exact equivalent of Part 3 values; replaces
#'     DPL-exact normalisation.}
#'   \item{`$sample_depth`}{Named integer vector. Number of series per year.}
#'   \item{`$absent`}{`data.frame` (`series`, `year`): all absent rings
#'     (zero values) detected in the raw data.}
#'   \item{`$stats`}{`data.frame` (one row per series): `series`, `jyr`,
#'     `lyr`, `n`, `n_segs`, `n_flags`, `r_master`, `mean_u`, `max_u`,
#'     `sd_u`, `acor_u`, `sens_u`, `max_f`, `sd_f`, `acor_f`, `ar_order`.
#'     Use `n_flags > 0` or `r_master < crit` to gate quality control.}
#'   \item{`$segments`}{`data.frame` (one row per segment): `series`,
#'     `seq_no`, `seg_start`, `seg_end`, `r_dated`, `r_max`, `lag_max`,
#'     `flag`. \code{"A"} means correlation below threshold; \code{"B"} means
#'     correlation peaks at a non-zero lag (possible dating error).}
#'   \item{`$problems`}{Named list (one element per series) with `$low_r`,
#'     `$corrxn`, `$cofdif`, `$absent`, `$outliers`. Passed to
#'     \code{\link{dpl_cof_diag}} for on-demand Part 6 output.}
#'   \item{`$crit`}{Numeric. Critical correlation threshold used.}
#'   \item{`$options`}{List of all argument values used in this run.}
#'   \item{`$output`}{Character vector of formatted text for requested parts.}
#' }
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # Run COFECHA with Holmes defaults
#' cof <- dpl_cof(rwl)
#'
#' # Inspect results
#' cof$stats[order(cof$stats$r_master), ]       # ranked by master correlation
#' subset(cof$segments, flag != "")              # flagged segments
#'
#' # Write full 7-part output to file
#' cof <- dpl_cof(rwl, output_file = tempfile(fileext = ".out"))
#'
#' # Light run (Parts 1 and 7 only), then on-demand diagnostics
#' cof <- dpl_cof(rwl, parts = c(1L, 7L), verbose = FALSE)
#' dpl_cof_diag(cof, series = "ACC014A")
#'
#' # Iterate: remove a low-correlation series and re-run
#' drop <- cof$stats$series[cof$stats$r_master < cof$crit]
#' rwl2 <- dpl_edt(rwl,
#'   edits  = lapply(drop, function(id) list(series = id, op = "omit")),
#'   verbose = FALSE
#' )
#' cof2 <- dpl_cof(rwl2, verbose = FALSE)
#' range(cof2$stats$r_master)
#' }
#'
#' @references
#' Holmes, R. L. (1983). Computer-assisted quality control in tree-ring dating
#' and measurement. *Tree-Ring Bulletin* 43:69--78.
#'
#' Holmes, R. L. (1994). *Dendrochronology Program Library*. Laboratory of
#' Tree-Ring Research, University of Arizona, Tucson.
#'
#' @seealso \code{\link{dpl_edt}}, \code{\link{dpl_cof_diag}}, \code{\link{dpl_dateme}}, \code{\link{dpl_barplot}}
#' @export
dpl_cof <- function(rwl,
                     seg_length    = 50L,
                     seg_lag       = 25L,
                     spline_period = 32L,
                     ar_model      = TRUE,
                     log_transform = TRUE,
                     crit_level    = NULL,
                     outp          =  3.0,
                     outn          = -4.5,
                     parts         = 1:7,
                     min_length    = 10L,
                     output_file   = NULL,
                     verbose       = TRUE) {

  # ---- Input validation ----------------------------------------------------
  if (!.is_rwl(rwl)) stop("'rwl' must be a dplR rwl data.frame.")
  min_length <- max(1L, as.integer(min_length))
  if (ncol(rwl) < 2L)
    stop("'rwl' must contain at least two series for a meaningful master.")
  seg_lag <- min(seg_lag, seg_length %/% 2L)
  crt     <- if (!is.null(crit_level)) crit_level else .cof_crit99(seg_length)

  years_all <- as.integer(rownames(rwl))
  iayr      <- min(years_all)
  izyr      <- max(years_all)
  ser_ids   <- colnames(rwl)
  nser_tot  <- length(ser_ids)

  # ---- Storage initialisation ----------------------------------------------
  yms  <- setNames(rep(0.0, izyr - iayr + 1L), as.character(iayr:izyr))
  ymsm <- yms           # AR-modelled master
  ysd  <- yms           # sd accumulator (for normalised master)
  jns  <- setNames(rep(0L, izyr - iayr + 1L), as.character(iayr:izyr))
  jab_yr <- jns         # absent rings per year (across all series)

  absent_df   <- data.frame(series = character(0), year = integer(0))
  stats_list  <- vector("list", nser_tot)
  seg_rows    <- list()
  problems    <- vector("list", nser_tot)
  names(problems) <- ser_ids

  # Scratch: filtered series for second loop (list of lists)
  ser_store   <- vector("list", nser_tot)

  # ---- FIRST LOOP: filter series, accumulate master -----------------------
  n_dated <- 0L
  nrtot   <- 0L    # total rings across all series
  short_df <- data.frame(series = character(0), jyr = integer(0),
                         lyr = integer(0), n = integer(0),
                         stringsAsFactors = FALSE)
  for (seq_no in seq_len(nser_tot)) {
    id   <- ser_ids[seq_no]
    col  <- rwl[[id]]
    ok   <- which(!is.na(col))
    if (length(ok) == 0L) next

    # Series shorter than min_length cannot be tested by segments (CRIT99 is
    # tabulated from n = 10) and would only add noise to the master: exclude
    # them from the run and report them in Part 1 and $short.
    n_ok <- length(ok)
    if (n_ok < min_length) {
      short_df <- rbind(short_df, data.frame(
        series = id, jyr = years_all[ok[1L]], lyr = years_all[ok[n_ok]],
        n = n_ok, stringsAsFactors = FALSE))
      if (verbose)
        message(sprintf("%4d  %-8s  %d - %d  (%d yr)  -- too short (< %d), not used",
                        seq_no, id, years_all[ok[1L]], years_all[ok[n_ok]],
                        n_ok, min_length))
      next
    }

    jyr  <- years_all[ok[1L]]
    lyr  <- years_all[ok[length(ok)]]
    n    <- lyr - jyr + 1L
    z    <- col[ok[1L]:ok[length(ok)]]
    z[is.na(z)] <- 0.0   # treat embedded NAs as absent rings

    # Raw statistics (unfiltered)
    zmean  <- mean(z, na.rm = TRUE)
    zmax   <- max(z, na.rm = TRUE)
    zsd_u  <- sd(z, na.rm = TRUE)
    zacor  <- if (n > 2L) .cof_correl(z[-n], z[-1L]) else NA_real_
    zsen   <- if (zmean > 0 && n > 1L) {
      mean(abs(2 * diff(z) / (z[-n] + z[-1L])), na.rm = TRUE)
    } else NA_real_

    # Absent rings
    ab_years <- jyr + which(z == 0.0 | z < -.5) - 1L
    if (length(ab_years) > 0L) {
      absent_df <- rbind(absent_df,
                         data.frame(series = id, year = as.integer(ab_years)))
      for (k in ab_years) {
        jab_yr[as.character(k)] <- jab_yr[as.character(k)] + 1L
        # replace zero with small positive for spline stability
        pos <- k - jyr + 1L
        if (z[pos] == 0.0) z[pos] <- -0.006
      }
    }

    nrtot <- nrtot + n
    n_dated <- n_dated + 1L

    # Cubic smoothing spline and variance stabilisation
    if (spline_period > 0L) {
      sp_curve <- .cof_spline(z, spline_period, 0.5)
      zf       <- .cof_divser(z, sp_curve)
      zf       <- .cof_varsta(zf, 32L)
    } else {
      zf <- .cof_normts(z, k = 0L)$z
    }

    # Fortran LOGTR, NORMTS, YMS, YMSM are all REAL (32-bit single precision).
    # Simulate this by casting through single precision at each step,
    # matching the numerical behaviour of the compiled Fortran binary exactly.
    sp <- function(v) as.double(as.single(v))

    # Log-transform for standard master accumulation
    # Fortran: LOGTR operates on REAL → result is REAL
    z_log  <- sp(if (log_transform) .cof_logtr(sp(zf), 1/3) else zf)
    # NORMTS on REAL input → REAL output
    z_norm <- sp(.cof_normts(z_log, k = 0L)$z)

    # Accumulate into master (standard = log-normalised); YMS is REAL
    for (j in seq_len(n)) {
      yr_ch <- as.character(jyr + j - 1L)
      jns[yr_ch]  <- jns[yr_ch] + 1L
      yms[yr_ch]  <- sp(yms[yr_ch] + z_norm[j])
      ysd[yr_ch]  <- sp(ysd[yr_ch] + z_norm[j]^2)
    }

    # AR modelling — MEMPR operates on REAL Z (variance-stabilised, single precision)
    ar_res <- if (ar_model && n >= 8L) .cof_mempr(sp(zf), lg = 10L) else NULL
    iar    <- if (!is.null(ar_res)) ar_res$order else 0L
    # Fortran: Z = WK(I,1) = demeaned AR residuals (REAL, single precision from MEMPR)
    # Then LOGTR (REAL) + NORMTS (REAL) applied. YMSM is REAL.
    z_ar   <- if (!is.null(ar_res)) ar_res$residuals else sp(zf)

    # Log-transform and normalise for AR master — both REAL in Fortran
    z_ar_log  <- sp(if (log_transform) .cof_logtr(z_ar, 1/3) else z_ar)
    z_ar_norm <- sp(.cof_normts(z_ar_log, k = 0L)$z)

    for (j in seq_len(n)) {
      yr_ch <- as.character(jyr + j - 1L)
      ymsm[yr_ch] <- sp(ymsm[yr_ch] + z_ar_norm[j])
    }

    # Filtered stats (from z_ar_log, in REAL space)
    zf_max  <- max(z_ar_log, na.rm = TRUE)
    zf_sd   <- sd(z_ar_log, na.rm = TRUE)
    zf_acor <- if (n > 2L) .cof_correl(z_ar_log[-n], z_ar_log[-1L]) else NA_real_

    # Store for second loop.
    # Fortran writes T (variance-stabilised series) to a direct-access file using
    # format F8.4 (4 decimal places), then reads it back in the second loop.
    # This rounding means ZSER/ZSERM in the second loop differ from what was
    # accumulated into YMS/YMSM in the first loop.
    # We simulate this by rounding zf to 4 decimal places, then recomputing
    # z_norm and z_ar_norm — these match Fortran's second-loop ZSER and ZSERM exactly.
    zf_f84 <- round(sp(zf), 4L)   # F8.4: 4 decimal places, single precision first

    # Recompute ZSER (standard) from rounded zf — matches Fortran second-loop
    z_log_f84  <- sp(if (log_transform) .cof_logtr(zf_f84, 1/3) else zf_f84)
    z_norm_f84 <- sp(.cof_normts(z_log_f84, k = 0L)$z)

    # Recompute ZSERM (AR) from rounded zf — matches Fortran second-loop
    ar_f84     <- if (ar_model && n >= 8L) .cof_mempr(zf_f84, lg = 10L) else NULL
    z_ar_f84   <- if (!is.null(ar_f84)) ar_f84$residuals else sp(zf_f84)
    z_arl_f84  <- sp(if (log_transform) .cof_logtr(z_ar_f84, 1/3) else z_ar_f84)
    z_arn_f84  <- sp(.cof_normts(z_arl_f84, k = 0L)$z)

    ser_store[[seq_no]] <- list(
      id     = id, jyr = jyr, lyr = lyr, n = n,
      z_norm    = z_norm_f84,   # ZSER:  from F8.4-rounded T, matches Fortran second loop
      z_ar_norm = z_arn_f84,    # ZSERM: from F8.4-rounded T, matches Fortran second loop
      iar    = iar,
      stats_u = c(mean = zmean, max = zmax, sd = zsd_u,
                   acor = zacor, sen = zsen),
      stats_f = c(max = zf_max, sd = zf_sd, acor = zf_acor)
    )

    if (verbose)
      message(sprintf("%4d  %-8s  %d - %d  (%d yr)", seq_no, id, jyr, lyr, n))
  }

  # Finalise master series (mean)
  ok_yrs <- which(jns > 0L)
  for (k in ok_yrs) {
    # Fortran: YMS(I)=YMS(I)/AJNS  (REAL division)
    yms[k]  <- as.double(as.single(yms[k]  / jns[k]))
    ymsm[k] <- as.double(as.single(ymsm[k] / jns[k]))
    if (jns[k] > 1L)
      # Fortran: YSD(I)=SQRT(ABS((YSD(I)-AJNS*YMS(I)*YMS(I))/(AJNS-1.)))  REAL
      ysd[k] <- as.double(as.single(
        sqrt(abs((ysd[k] - jns[k] * yms[k]^2) / (jns[k] - 1L)))))
    else
      ysd[k] <- 0.0
  }
  # SDM = mean SD of master — Fortran lines 607, 615-616:
  #   SDM = SDM + YSD(I)   [for years with JNS>1]
  #   N   = IZM - IAM + 1  [FULL master span, including JNS=1 years]
  #   SDM = SDM / FLOAT(N)
  # Must divide by full master length, NOT by count(JNS>1).
  n_master_full <- if (length(ok_yrs) > 0L) max(ok_yrs) - min(ok_yrs) + 1L else 1L
  sdm <- if (sum(jns > 1L) > 0L) sum(ysd[jns > 1L]) / n_master_full else 0.0

  # Master span with >= 1 series; span with >= 2
  iam <- min(as.integer(names(jns)[jns >= 1L]))
  izm <- max(as.integer(names(jns)[jns >= 1L]))
  ia2 <- if (any(jns >= 2L)) min(as.integer(names(jns)[jns >= 2L])) else iam
  iz2 <- if (any(jns >= 2L)) max(as.integer(names(jns)[jns >= 2L])) else izm

  if (verbose) {
    message(sprintf("\n Time span of Master dating series: %d to %d", iam, izm))
    message(sprintf(" Portion with two or more series:   %d to %d", ia2, iz2))
  }

  # Trim master to its actual span
  master_yrs <- as.character(iam:izm)
  master     <- yms[master_yrs]      # mean of z_norm per series (REAL)
  master_ar  <- ymsm[master_yrs]
  depth      <- jns[master_yrs]
  absent_by_yr <- jab_yr[master_yrs]

  # Fortran: YSD(J)=YMS(J) then CALL NORMTS(N,YSD(IAM),YNF,XSD,0)
  # Both YSD and NORMTS are REAL — simulate single precision on the global normalisation.
  ysd_sp  <- as.double(as.single(as.numeric(master)))  # YSD = YMS copy (REAL)
  nr_sp   <- .cof_normts(ysd_sp, k = 0L)               # NORMTS on REAL input
  master_norm <- setNames(as.double(as.single(nr_sp$z)), master_yrs)

  # Series meta for Part 2
  series_meta <- lapply(seq_len(nser_tot), function(i) {
    s <- ser_store[[i]]
    if (is.null(s)) return(NULL)
    list(id = s$id, seq = i, jyr = s$jyr, n = s$n)
  })
  series_meta <- Filter(Negate(is.null), series_meta)

  # ---- SECOND LOOP: segment correlations and diagnostics ------------------
  nrchq   <- 0L    # total rings in checked portions
  ravt    <- 0.0   # running sum for mean intercorrelation

  for (seq_no in seq_len(nser_tot)) {
    s <- ser_store[[seq_no]]
    if (is.null(s)) next

    id  <- s$id;  jyr <- s$jyr;  lyr <- s$lyr;  n <- s$n
    z_norm   <- s$z_norm      # log-normalised filtered series (ZSER in Fortran)
    z_ar_nm  <- s$z_ar_norm   # AR normalised series (ZSERM in Fortran)

    # Remove this series from master.
    # Fortran: YMSA(L)  = (YMS(L)*JNS(L) - ZSER(L))  / (JNS(L)-1)
    #          YMSMA(L) = (YMSM(L)*JNS(L) - ZSERM(L)) / (JNS(L)-1)
    # All arrays are REAL → single-precision arithmetic.
    ymsa  <- master      # master-minus-self, standard (REAL)
    ymsma <- master_ar   # master-minus-self, AR (REAL)

    for (j in jyr:lyr) {
      yr_ch <- as.character(j)
      if (jns[yr_ch] > 1L) {
        n_j <- jns[yr_ch]
        ymsa[yr_ch]  <- as.double(as.single(
          (yms[yr_ch] * n_j - z_norm[j - jyr + 1L]) / (n_j - 1L)))
        ymsma[yr_ch] <- as.double(as.single(
          (ymsm[yr_ch] * n_j - z_ar_nm[j - jyr + 1L]) / (n_j - 1L)))
      }
    }

    # Determine checkable portion
    jyr1 <- max(jyr, ia2);  lyr1 <- min(lyr, iz2)
    n_check <- max(0L, lyr1 - jyr1 + 1L)
    nrchq   <- nrchq + n_check

    # Segment starting years
    segs_starts <- .cof_qseg(jyr1, lyr1, seg_length, seg_lag)

    # Run .cof_slsg for each segment
    seg_results <- vector("list", length(segs_starts))
    n_flagged   <- 0L
    low_r_segs  <- list()

    for (ki in seq_along(segs_starts)) {
      js   <- segs_starts[ki]
      je   <- min(js + seg_length - 1L, lyr1)
      lsk  <- je - js + 1L

      # Extract master-minus-self window and test series window
      ma_vec   <- as.numeric(ymsma[as.character(iam:izm)])
      test_vec <- z_ar_nm[js - jyr + 1L : (je - jyr + 1L)]

      sl <- .cof_slsg(
        master_seg = ma_vec,
        test_seg   = z_ar_nm[(js - jyr + 1L):(je - jyr + 1L)],
        jyrs       = js,
        jyrm       = iam,
        lyrm       = izm,
        lsg        = lsk,
        crt        = crt
      )
      sl$start <- as.integer(js)
      sl$end   <- as.integer(je)
      seg_results[[ki]] <- sl

      if (sl$flag != "") {
        n_flagged <- n_flagged + 1L
        low_r_segs <- c(low_r_segs, list(sl))
      }

      seg_rows <- c(seg_rows, list(data.frame(
        series    = id,
        seq_no    = as.integer(seq_no),
        seg_start = as.integer(js),
        seg_end   = as.integer(je),
        r_dated   = sl$r_dated,
        r_max     = sl$r_max,
        lag_max   = as.integer(sl$lag_max),
        flag      = sl$flag,
        stringsAsFactors = FALSE
      )))
    }

    # Overall correlation with master-minus-self (Fortran line 824):
    #   CORREL(N, YMSMA(JYR), ZSERM(JYR), RAV)
    # where N and JYR are the CLIPPED values (after lines 731-736).
    # Must use jyr1/lyr1, not max(jyr,iam)/min(lyr,izm).
    if (n_check >= 2L) {
      ma_ov <- as.numeric(ymsma[as.character(jyr1:lyr1)])
      ts_ov <- z_ar_nm[(jyr1 - jyr + 1L):(lyr1 - jyr + 1L)]
      rav   <- .cof_correl(ma_ov, ts_ov)
      ravt  <- ravt + rav * n_check
    } else {
      rav   <- NA_real_
      ma_ov <- numeric(0L)
      ts_ov <- numeric(0L)
    }

    # Influence analysis (Fortran line 760: CORRXN uses same clipped JYR/N)
    cx_full <- if (n_check >= 8L && length(ma_ov) >= 8L) {
      .cof_corrxn(ma_ov, ts_ov, jyr1)
    } else NULL

    # First-difference divergence (Fortran: COFDIF uses YMSA and ZSER = z_norm)
    cd <- if (n_check > 2L) {
      .cof_cofdif(
        as.numeric(ymsa[as.character(jyr1:lyr1)]),
        z_norm[(jyr1 - jyr + 1L):(lyr1 - jyr + 1L)],
        jyr1, lyr1
      )
    } else data.frame(year = integer(0), delta = numeric(0))

    # Outliers (Fortran: OUTABS uses ZSER = z_norm and YMSA)
    ot <- if (n_check > 0L) {
      .cof_outabs(
        z_norm[(jyr1 - jyr + 1L):(lyr1 - jyr + 1L)],
        as.numeric(ymsa[as.character(jyr1:lyr1)]),
        sdm, jyr1, lyr1, outp, outn
      )
    } else data.frame(year = integer(0), zsd = numeric(0))

    # Absent rings for this series
    abs_this <- absent_df$year[absent_df$series == id]

    problems[[id]] <- list(
      low_r    = low_r_segs,
      corrxn   = cx_full,
      cofdif   = cd,
      absent   = as.integer(abs_this),
      outliers = ot,
      outp     = outp,    # fixed threshold for [E] header
      outn     = outn     # fixed threshold for [E] header
    )

    # Stats row
    stats_list[[seq_no]] <- data.frame(
      series   = id,
      jyr      = jyr,
      lyr      = lyr,
      n        = n,
      n_segs   = length(segs_starts),
      n_flags  = n_flagged,
      r_master = rav,
      mean_u   = s$stats_u["mean"],
      max_u    = s$stats_u["max"],
      sd_u     = s$stats_u["sd"],
      acor_u   = s$stats_u["acor"],
      sens_u   = s$stats_u["sen"],
      max_f    = s$stats_f["max"],
      sd_f     = s$stats_f["sd"],
      acor_f   = s$stats_f["acor"],
      ar_order = s$iar,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      lrav  <- max(0L, round(rav * 20))
      message(sprintf("%4d  %-8s  %d - %d  %d yr  r=%.3f  %s",
                      seq_no, id, jyr, lyr, n,
                      ifelse(is.na(rav), -9.99, rav),
                      strrep("]", lrav)))
    }
  }

  # Final summary statistics
  mean_r   <- if (nrchq > 0L) ravt / nrchq else NA_real_
  mean_sen <- if (nrow(absent_df) > 0L || n_dated > 0L) {
    mean(vapply(stats_list, function(r) if (!is.null(r)) r$sens_u else NA_real_,
                numeric(1)), na.rm = TRUE)
  } else NA_real_
  n_prob <- sum(vapply(stats_list,
                       function(r) if (!is.null(r)) r$n_flags else 0L, integer(1)))

  if (verbose) {
    message(sprintf(
      "\n %s\n *C* Number of dated series%10d *C*\n *O* Master series %4d %4d %4d yrs *O*\n *F* Total rings in all series%7d *F*\n *E* Total dated rings checked%7d *E*\n *C* Series intercorrelation%9.3f *C*\n *H* Average mean sensitivity%8.3f *H*\n *A* Segments, possible problems%5d *A*\n %s",
      strrep("*", 40L),
      n_dated, iam, izm, izm - iam + 1L,
      nrtot, nrchq,
      ifelse(is.na(mean_r), -9.99, mean_r),
      ifelse(is.na(mean_sen), 0.0, mean_sen),
      n_prob,
      strrep("*", 40L)
    ))
  }

  # ---- Build stats and segments data.frames --------------------------------
  stats_df <- do.call(rbind, Filter(Negate(is.null), stats_list))
  rownames(stats_df) <- NULL

  segs_df <- if (length(seg_rows) > 0L) {
    do.call(rbind, seg_rows)
  } else {
    data.frame(series = character(0), seq_no = integer(0),
               seg_start = integer(0), seg_end = integer(0),
               r_dated = numeric(0), r_max = numeric(0),
               lag_max = integer(0), flag = character(0))
  }
  rownames(segs_df) <- NULL

  # ---- Formatted text output -----------------------------------------------
  out_lines <- character(0)
  title_str <- sprintf("xDPL COFECHA  %s", format(Sys.Date(), "%Y-%m-%d"))

  if (1L %in% parts) {
    mean_len <- mean(vapply(Filter(Negate(is.null), ser_store),
                            function(s) s$n, numeric(1)))
    p1 <- c(
      sprintf("PART 1:  %s", title_str),
      strrep("-", 132L),
      "",
      " QUALITY CONTROL AND DATING CHECK OF TREE-RING MEASUREMENTS",
      "",
      " RUN CONTROL OPTIONS SELECTED                             VALUE",
      "",
      sprintf("         1  Cubic smoothing spline 50%% wavelength cutoff for filtering"),
      sprintf("                                                            %d years", spline_period),
      sprintf("         2  Segments examined are                           %d years lagged successively by %3d years",
              seg_length, seg_lag),
      sprintf("         3  Autoregressive model %s",
              if (ar_model) "applied                     A  Residuals are used in master dating series and testing"
              else "not applied                  N"),
      sprintf("         4  Series %stransformed to logarithms%s%s",
              if (log_transform) "" else "not ",
              if (log_transform) "                 " else "              ",
              if (log_transform) " Y  Each series log-transformed for master dating series and testing" else " N"),
      sprintf("         5  Critical correlation, 99%% confidence level   %.4f", crt),
      "",
      sprintf("                                        %s", strrep("*", 40L)),
      sprintf("                                        *C* Number of dated series%10d *C*", n_dated),
      sprintf("                                        *O* Master series%5d%5d%5d yrs *O*", iam, izm, izm - iam + 1L),
      sprintf("                                        *F* Total rings in all series%7d *F*", nrtot),
      sprintf("                                        *E* Total dated rings checked%7d *E*", nrchq),
      sprintf("                                        *C* Series intercorrelation%9.3f *C*",
              ifelse(is.na(mean_r), -9.99, mean_r)),
      sprintf("                                        *H* Average mean sensitivity%8.3f *H*",
              ifelse(is.na(mean_sen), 0.0, mean_sen)),
      sprintf("                                        *A* Segments, possible problems%5d *A*", n_prob),
      sprintf("                                        *** Mean length of series%10.1f ***", mean_len),
      sprintf("                                        %s", strrep("*", 40L))
    )
    if (nrow(absent_df) > 0L) {
      total_ab <- nrow(absent_df)
      pct_ab   <- 100 * total_ab / nrtot
      p1 <- c(p1, "",
              " ABSENT RINGS listed by SERIES:            (See Master Dating Series for absent rings listed by year)",
              "")
      for (sid in unique(absent_df$series)) {
        ab <- absent_df$year[absent_df$series == sid]
        p1 <- c(p1,
          sprintf(" %-8s  %d absent rings:  %s",
                  sid, length(ab),
                  paste(sprintf("%5d", ab), collapse = "")))
      }
      p1 <- c(p1, "",
        sprintf("%14d absent rings%8.3f%%", total_ab, pct_ab))
    }
    if (nrow(short_df) > 0L) {
      p1 <- c(p1, "",
        sprintf(" SERIES NOT USED: shorter than %d years (not included in master, not tested by segments)",
                min_length), "")
      for (k in seq_len(nrow(short_df)))
        p1 <- c(p1, sprintf(" %-8s  %4d to %4d  %4d years",
                            short_df$series[k], short_df$jyr[k],
                            short_df$lyr[k], short_df$n[k]))
    }
    out_lines <- c(out_lines, p1)
  }

  if (2L %in% parts)
    out_lines <- c(out_lines, .cof_fmt_part2(series_meta, title_str))

  if (3L %in% parts)
    # Fortran passes YSD which is YMS(J) normalized globally — this is master_norm
    out_lines <- c(out_lines,
      .cof_fmt_part3(master_norm, depth, absent_by_yr, title_str))

  if (4L %in% parts)
    # BARPL also receives YSD (same normalised master); NORMTS applied internally
    out_lines <- c(out_lines, .cof_fmt_part4(master_norm, title_str))

  # ---- Part 5: build full table then emit ----------------------------------
  # Collect all per-series rows first (needed for column header and avg row)
  # Pre-sized so that series skipped in the first loop (empty, or shorter
  # than min_length) leave NULL slots; indexing p5_rows[[seq_no]] past the
  # last stored series would otherwise be out of bounds.
  p5_rows     <- vector("list", nser_tot)   # .cof_fmt_seg_table() results
  p5_all_segs <- vector("list", nser_tot)   # all segment results (for avg)

  for (seq_no in seq_len(nser_tot)) {
    s <- ser_store[[seq_no]]
    if (is.null(s)) next
    id  <- s$id
    # Rebuild full segment list from segs_df (all segments, flagged or not)
    this_segs_df <- segs_df[segs_df$series == id, ]
    full_segs <- lapply(seq_len(nrow(this_segs_df)), function(k) {
      list(start   = this_segs_df$seg_start[k],
           end     = this_segs_df$seg_end[k],
           r_dated = this_segs_df$r_dated[k],
           r_max   = this_segs_df$r_max[k],
           lag_max = this_segs_df$lag_max[k],
           flag    = this_segs_df$flag[k])
    })
    # SGTBL receives IA2/IZ2 (Fortran line 822), not IAM/IZM.
    # jyr_clipped = max(orig_jyr, ia2) = the clipped start used in QSEG.
    jyr_clipped_s <- max(as.integer(s$jyr), as.integer(ia2))
    p5r <- .cof_fmt_seg_table(id, seq_no, s$jyr, s$lyr, s$n,
                               full_segs, iam, izm, ia2, iz2,
                               seg_lag, crt, jyr_clipped_s)
    p5_rows[[seq_no]]    <- p5r
    p5_all_segs[[seq_no]] <- full_segs
  }

  if (5L %in% parts) {
    # SGTBL uses IA2/IZ2 for the display grid (Fortran line 822).
    # IFA = ((IA2-1)/LAG)*LAG;  NSP = (IZ2-IA2)/LAG + 1
    NSP <- (iz2 - ia2) %/% seg_lag + 1L
    IFA <- ((ia2 - 1L) %/% seg_lag) * seg_lag
    seg_hdr_yrs <- IFA + seg_lag * seq(0L, NSP - 1L)
    seg_hdr_end <- seg_hdr_yrs + seg_length - 1L

    out_lines <- c(out_lines,
      sprintf("PART 5:  CORRELATION OF SERIES BY SEGMENTS: %s", title_str),
      strrep("-", 132L),
      sprintf(" Correlations of%4d-year dated segments, lagged%4d years",
              seg_length, seg_lag),
      sprintf(" Flags:  __A = correlation under%8.4f but highest as dated;  __B = correlation higher at other than dated position",
              crt),
      "",
      # Fortran line 1355: ' Seq Series   Interval',T25,20I5
      # ' Seq Series   Interval' = 22 chars; T25 pads to col 25 → 3 trailing spaces
      sprintf(" Seq Series   Interval   %s",
              paste(sprintf("%5d", seg_hdr_yrs), collapse = "")),
      # Interval end-years: T25,20I5 → 25-char blank prefix then years
      sprintf("%25s%s",
              "",
              paste(sprintf("%5d", seg_hdr_end), collapse = "")),
      # Fortran line 1357-1358: ' --- -------- ---------',T25,20(' ----')
      # ' --- -------- ---------' = 23 chars; T25 pads to 25 → 2 trailing spaces; then 20×' ----'
      sprintf(" --- -------- ---------  %s",
              paste(rep(" ----", NSP), collapse = ""))
    )

    for (seq_no in seq_len(nser_tot)) {
      pr <- p5_rows[[seq_no]]
      if (is.null(pr)) next
      out_lines <- c(out_lines, pr$lines)
    }

    # Average segment correlation — accumulate using each series' ISEGA offset
    avg_r <- numeric(NSP)
    n_avg <- integer(NSP)
    for (seq_no in seq_len(nser_tot)) {
      pr        <- p5_rows[[seq_no]]
      segs_this <- p5_all_segs[[seq_no]]
      if (is.null(pr) || is.null(segs_this)) next
      isega_s <- pr$isega
      for (k in seq_along(segs_this)) {
        col <- isega_s + k
        if (col < 1L || col > NSP) next
        r <- segs_this[[k]]$r_dated
        if (!is.na(r) && r > -9.0) {
          avg_r[col] <- avg_r[col] + r
          n_avg[col] <- n_avg[col] + 1L
        }
      }
    }
    avg_r <- ifelse(n_avg > 0L, avg_r / n_avg, NA_real_)
    avg_str <- vapply(avg_r, function(r) {
      if (is.na(r)) return("   - ")
      if (r < 0.0) {
        s <- sprintf("%.2f", -r); paste0("-", substr(s, 2L, 4L), " ")
      } else {
        s <- sprintf("%.2f",  r); paste0(" ", substr(s, 2L, 4L), " ")
      }
    }, character(1))
    # ' Av segment correlation' = 23 chars; pad to 25 with 2 trailing spaces
    out_lines <- c(out_lines,
      sprintf(" Av segment correlation  %s",
              paste(avg_str, collapse = "")))
  }

  # ---- Part 6 -----------------------------------------------------------
  if (6L %in% parts) {
    # Fortran Part 6 header (lines 664-676): note '[B]' line also mentions symbol key
    out_lines <- c(out_lines,
      sprintf("PART 6:  POTENTIAL PROBLEMS: %s", title_str),
      strrep("-", 132L),
      "",
      " For each series with potential problems the following diagnostics may appear:",
      "",
      sprintf(" [A] Correlations with master dating series of flagged%4d-year segments of series filtered with%4d-year spline,",
              seg_length, spline_period),
      "     at every point from ten years earlier (-10) to ten years later (+10) than dated",
      "",
      " [B] Effect of those data values which most lower or raise correlation with master series",
      "     Symbol following year indicates value in series is greater (>) or lesser (<) than master series value",
      "",
      " [C] Year-to-year changes very different from the mean change in other series",
      "",
      " [D] Absent rings (zero values)",
      "",
      " [E] Values which are statistical outliers from mean for the year",
      paste0(" ", strrep("=", 131L)))

    for (seq_no in seq_len(nser_tot)) {
      s <- ser_store[[seq_no]]
      if (is.null(s)) next
      id  <- s$id
      p6 <- .cof_fmt_part6_series(
        id, seq_no, s$jyr, s$lyr, s$n,
        problems[[id]], crt, seg_length, seg_lag,
        master_ysd   = master_norm,
        depth        = depth,
        absent_by_yr = absent_by_yr
      )
      out_lines <- c(out_lines, p6)
    }

    # Final [*] message if no segments had problems
    if (n_prob == 0L)
      out_lines <- c(out_lines,
        sprintf("\n [*] All segments correlate highest as dated, with correlation with master series over%8.4f",
                crt))
  }

  # ---- Part 7 -----------------------------------------------------------
  if (7L %in% parts) {
    out_lines <- c(out_lines,
      sprintf("PART 7:  DESCRIPTIVE STATISTICS: %s", title_str),
      strrep("-", 132L),
      .cof_fmt_part7_header())

    tot_n <- 0L; tot_segs <- 0L; tot_flags <- 0L
    sum_r <- 0.0; sum_mean <- 0.0; max_max <- 0.0
    sum_sd <- 0.0; sum_acor <- 0.0; sum_sen <- 0.0
    sum_maxf <- 0.0; sum_sdf <- 0.0; sum_acorf <- 0.0

    for (seq_no in seq_len(nser_tot)) {
      r <- stats_list[[seq_no]]
      if (is.null(r)) next
      # n_flags in Part 7 = NLO + NBE (flags A + flags B, total problems)
      nund <- r$n_flags  # already NLO+NBE from second loop
      tot_n     <- tot_n     + r$n
      tot_segs  <- tot_segs  + r$n_segs
      tot_flags <- tot_flags + nund
      sum_r     <- sum_r     + r$r_master * r$n
      sum_mean  <- sum_mean  + r$mean_u   * r$n
      max_max   <- max(max_max, r$max_u)
      sum_sd    <- sum_sd    + r$sd_u     * r$n
      sum_acor  <- sum_acor  + r$acor_u   * r$n
      sum_sen   <- sum_sen   + r$sens_u   * r$n
      sum_maxf  <- sum_maxf  + r$max_f    * r$n
      sum_sdf   <- sum_sdf   + r$sd_f     * r$n
      sum_acorf <- sum_acorf + r$acor_f   * r$n
      out_lines <- c(out_lines,
        .cof_fmt_part7_row(
          seq_no, r$series, r$jyr, r$lyr, r$n,
          r$n_segs, nund, r$r_master,
          c(mean = r$mean_u, max = r$max_u, sd = r$sd_u,
            acor = r$acor_u, sen = r$sens_u),
          c(max  = r$max_f,  sd = r$sd_f,  acor = r$acor_f),
          r$ar_order
        )
      )
    }

    # Separator + Total row — use shared .cof_fXX format helpers (no local closures)

    # Separator row from Fortran: ' --- -------- ---------', 3x'  -----', '   ------', 8x' ----- ', ' --'
    sep_row <- paste0(" --- -------- ---------",
                      "  -----  -----  -----",
                      "   ------",
                      paste(rep(" ----- ", 8L), collapse = ""),
                      " --")
    r_tot <- if (tot_n > 0L) sum_r / tot_n else 0.0
    out_lines <- c(out_lines,
      sep_row,
      paste0(
        sprintf(" Total or mean:        %7d%7d%7d",
                tot_n, tot_segs, tot_flags),
        .cof_f83(r_tot),
        .cof_f72(if (tot_n > 0L) sum_mean  / tot_n else 0.0),
        .cof_f72(max_max),
        .cof_f73(if (tot_n > 0L) sum_sd    / tot_n else 0.0),
        .cof_f73(if (tot_n > 0L) sum_acor  / tot_n else 0.0),
        .cof_f73(if (tot_n > 0L) sum_sen   / tot_n else 0.0),
        .cof_f72(if (tot_n > 0L) sum_maxf  / tot_n else 0.0),
        .cof_f73(if (tot_n > 0L) sum_sdf   / tot_n else 0.0),
        .cof_f73(if (tot_n > 0L) sum_acorf / tot_n else 0.0)
      )
    )
  }

  out_lines <- c(out_lines, "",
    sprintf("                                              - = [ COFECHA ] = -"))

  # Write to file if requested
  if (!is.null(output_file))
    writeLines(out_lines, output_file)

  # ---- Return --------------------------------------------------------------
  # ---- Filtered series as an rwl (what Part 5 correlates: spline, log, AR,
  # normalised, i.e. ZSERM) on the master's year axis ------------------------
  filt <- as.data.frame(matrix(NA_real_, length(master_yrs), 0L,
                               dimnames = list(master_yrs, NULL)))
  for (seq_no in seq_len(nser_tot)) {
    s <- ser_store[[seq_no]]
    if (is.null(s)) next
    col <- rep(NA_real_, length(master_yrs))
    col[match(as.character(s$jyr:s$lyr), master_yrs)] <- s$z_ar_norm
    filt[[s$id]] <- col
  }
  class(filt) <- c("rwl", "data.frame")

  result <- list(
    short        = short_df,      # series excluded as shorter than min_length
    filtered     = filt,          # COFECHA-transformed series (Part 5 input)
    master       = master_norm,   # globally normalised master (for plotting, Part 3/4)
    master_raw   = master,        # raw mean of z_norm per series (pre global normts)
    sample_depth = depth,
    absent       = absent_df,
    stats        = stats_df,
    segments     = segs_df,
    problems     = problems,
    crit         = crt,
    options      = list(
      seg_length    = seg_length,
      seg_lag       = seg_lag,
      spline_period = spline_period,
      ar_model      = ar_model,
      log_transform = log_transform,
      crit_level    = crit_level,
      outp          = outp,
      outn          = outn,
      min_length    = min_length
    ),
    output       = out_lines
  )

  if (!is.null(output_file)) invisible(result) else result
}


#' Print Part 5/6 diagnostics for selected series
#'
#' @description
#' Prints the Part 6 potential-problems diagnostic block for one or more named
#' series from a completed \code{\link{dpl_cof}} run. Output is generated on-demand from
#' `$problems` --- no recomputation required. The same formatter is used
#' internally by \code{\link{dpl_cof}} for file output, ensuring consistency.
#'
#' @param cof_result The list returned by \code{\link{dpl_cof}} (or \code{\link{dpl_dateme}}).
#'   Must contain `$problems`, `$stats`, `$crit`, and `$options`.
#' @param series Character vector of series IDs, integer vector of sequence
#'   numbers, or `NULL` (default = all series). Exact matching; no partial
#'   matching applied.
#' @param parts Integer vector. `5` = segment correlation table (+/-10-year lag
#'   table); `6` = potential problems (\[A\]--\[E\]). Default `6`.
#'   Request both with `c(5, 6)`.
#' @param to_file Character or `NULL`. If supplied, formatted lines are also
#'   appended to that file. Default `NULL` (console only).
#'
#' @return Invisibly returns a character vector of all lines printed.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#' cof <- dpl_cof(rwl, parts = c(1L, 7L), verbose = FALSE)
#'
#' # Part 6 for the lowest-correlating series
#' worst <- cof$stats$series[which.min(cof$stats$r_master)]
#' dpl_cof_diag(cof, series = worst)
#'
#' # Part 5 + 6 for all flagged series
#' flagged <- unique(subset(cof$segments, flag != "")$series)
#' dpl_cof_diag(cof, series = flagged, parts = c(5L, 6L))
#'
#' # Save diagnostics to file
#' dpl_cof_diag(cof, series = worst,
#'              to_file = tempfile(fileext = "_diag.txt"))
#' }
#'
#' @seealso \code{\link{dpl_cof}}, \code{\link{dpl_dateme}}
#' @export
dpl_cof_diag <- function(cof_result,
                           series  = NULL,
                           parts   = 6L,
                           to_file = NULL) {

  if (is.null(cof_result$problems))
    stop("'cof_result' does not contain $problems. Was it produced by dpl_cof()?")

  all_ids <- names(cof_result$problems)

  # Resolve series selection
  if (is.null(series)) {
    sel_ids <- all_ids
  } else if (is.numeric(series)) {
    sel_ids <- all_ids[as.integer(series)]
  } else {
    missing_ids <- setdiff(series, all_ids)
    if (length(missing_ids) > 0L)
      warning("Series not found in cof_result: ",
              paste(missing_ids, collapse = ", "))
    sel_ids <- intersect(series, all_ids)
  }

  if (length(sel_ids) == 0L) {
    message("No matching series found.")
    return(invisible(character(0)))
  }

  crt    <- cof_result$crit
  lsg    <- cof_result$options$seg_length
  stats  <- cof_result$stats

  all_lines <- character(0)

  for (id in sel_ids) {
    sr  <- stats[stats$series == id, ]
    prb <- cof_result$problems[[id]]

    if (5L %in% parts) {
      p5 <- .cof_fmt_seg_table(
        id, sr$seq_no[1L], sr$jyr[1L], sr$lyr[1L], sr$n[1L],
        prb$low_r, NA_integer_, NA_integer_, cof_result$options$seg_lag, crt
      )
      all_lines <- c(all_lines, p5$lines)
    }

    if (6L %in% parts) {
      p6 <- .cof_fmt_part6_series(
        id, which(all_ids == id), sr$jyr[1L], sr$lyr[1L], sr$n[1L],
        prb, crt, lsg, cof_result$options$seg_lag
      )
      all_lines <- c(all_lines, p6)
    }
  }

  cat(paste(all_lines, collapse = "\n"), "\n")
  if (!is.null(to_file))
    cat(paste(all_lines, collapse = "\n"), "\n", file = to_file, append = TRUE)

  invisible(all_lines)
}


#' Date floating (undated) series against a COFECHA master
#'
#' @description
#' Implements Holmes' UDATE algorithm: slides segments of each undated series
#' along the dated master chronology from a completed \code{\link{dpl_cof}} run and
#' finds positions with the highest correlation. For each segment the 11
#' best-correlating positions are reported with their year adjustments and
#' correlations. A tally across segments gives a consensus best-fit date.
#'
#' @param rwl_undated A dplR `rwl` data.frame of undated (floating) series.
#'   Row names are nominal years defining series length only.
#' @param cof_result The list returned by \code{\link{dpl_cof}} on the dated series.
#'   Must contain `$master`, `$options`, and `$crit`.
#' @param seg_length Integer or `NULL`. Segment length for sliding-window
#'   correlation. `NULL` inherits from `cof_result$options$seg_length`.
#' @param seg_lag Integer or `NULL`. Lag between successive segments. `NULL`
#'   inherits from `cof_result$options$seg_lag`.
#' @param n_best Integer. Best-correlating positions to report per segment.
#'   Holmes default: `11`.
#' @param parts Integer vector. `8` = date-adjustment table (Part 8 in DPL
#'   terminology). Default `8`.
#' @param output_file Character or `NULL`. If supplied, writes Part 8 to file.
#' @param verbose Logical. Print per-series progress. Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`$summary`}{`data.frame` (one row per series): `series`, `n`,
#'     `best_adj` (year adjustment with most segment support),
#'     `n_segs_supporting`, `mean_r_at_best`.}
#'   \item{`$problems`}{Named list compatible with \code{\link{dpl_cof_diag}}.}
#'   \item{`$output`}{Character vector of Part 8 formatted text.}
#' }
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # Build a dated master from the full set
#' cof <- dpl_cof(rwl, parts = integer(0), verbose = FALSE)
#'
#' # Simulate an undated set by stripping the dates from two series
#' rwl_undated <- rwl[, c("ACC013B", "ACC013C")]
#' rownames(rwl_undated) <- as.character(seq_len(nrow(rwl_undated)))
#'
#' dat <- dpl_dateme(rwl_undated, cof_result = cof)
#' dat$summary    # best_adj = year adjustment with most segment support
#' }
#'
#' @seealso \code{\link{dpl_cof}}, \code{\link{dpl_cof_diag}}
#' @export
dpl_dateme <- function(rwl_undated,
                       cof_result,
                       seg_length  = NULL,
                       seg_lag     = NULL,
                       n_best      = 11L,
                       parts       = 8L,
                       output_file = NULL,
                       verbose     = TRUE) {

  if (!.is_rwl(rwl_undated))
    stop("'rwl_undated' must be a dplR rwl data.frame.")
  if (is.null(cof_result$master))
    stop("'cof_result' must be the output of dpl_cof().")

  ls  <- as.integer(if (!is.null(seg_length)) seg_length
                    else cof_result$options$seg_length)
  lag <- as.integer(if (!is.null(seg_lag))   seg_lag
                    else cof_result$options$seg_lag)
  n_best <- as.integer(n_best)

  # Master AR chronology (YMSM in Fortran)
  master_yrs <- as.integer(names(cof_result$master))
  iam  <- min(master_yrs)   # IA
  izm  <- max(master_yrs)   # IZ
  mast <- as.numeric(cof_result$master)
  nm   <- length(mast)

  years_u <- as.integer(rownames(rwl_undated))
  ser_ids <- colnames(rwl_undated)

  # ---- Helpers ---------------------------------------------------------------

  # Fortran F4.2 with leading-zero suppression for correlation values (lines 1768-1771)
  # WRITE(LN(K),'(F4.2)') RMX(K); IF(LN(K)(1:1).EQ.'0') LN(K)(1:1)=' '
  fmt_ln <- function(r) {
    if (is.na(r) || r <= -9.0) return("   -")   # -9.99 → missing
    s <- sprintf("%.2f", min(r, 9.99))           # e.g. "0.73"
    if (substr(s, 1L, 1L) == "0") substr(s, 1L, 1L) <- " "   # " .73"
    substr(s, 1L, 4L)                             # exactly 4 chars
  }

  # ---- Output header (printed once for the first series) --------------------
  # Fortran lines 1707-1710: PAGE call + time-span line
  title_str <- if (is.null(cof_result$options$title)) "" else cof_result$options$title
  out_lines  <- character(0)
  nsr        <- 0L   # series counter (NSR)
  summary_list <- list()

  # ---- Per-series loop (Fortran: GOTO 1 at line 1812) -----------------------
  for (sid in ser_ids) {

    col <- rwl_undated[[sid]]
    ok  <- which(!is.na(col))
    if (length(ok) == 0L) next

    # Series span and filtered values
    i0 <- ok[1L];  i1 <- ok[length(ok)]
    n  <- i1 - i0 + 1L
    z  <- col[i0:i1];  z[is.na(z)] <- 0.0
    jyr <- years_u[i0]   # JYR (nominal first year)

    nsr <- nsr + 1L
    if (verbose)
      message(sprintf("%6d Undated %s%6d years", nsr, sid, n))

    # Filter (same pipeline as dated series in dpl_cof)
    sp_curve <- .cof_spline(z, cof_result$options$spline_period, 0.5)
    zf <- .cof_divser(z, sp_curve)
    zf <- .cof_varsta(zf, 32L)
    if (cof_result$options$log_transform) zf <- .cof_logtr(zf, 1/3)
    zf <- .cof_normts(zf, k = 0L)$z
    if (cof_result$options$ar_model && n >= 8L) {
      ar_res <- .cof_mempr(zf, lg = 10L)
      zf <- .cof_normts(ar_res$residuals, k = 0L)$z
    }

    # Page header — written once for the first series (Fortran lines 1706-1710)
    if (nsr == 1L && 8L %in% parts) {
      out_lines <- c(out_lines,
        sprintf("PART 8:  DATE ADJUSTMENT FOR UNKNOWN SERIES: %s", title_str),
        strrep("-", 132L),
        "",
        # '(/'' Time span'',2I6,'', best matches for'',I4,''-year segments lagged'',I4,'' years''/)'
        sprintf(" Time span%6d%6d, best matches for%4d-year segments lagged%4d years",
                iam, izm, ls, lag),
        "")
    }

    # Per-series 3-line column header (Fortran lines 1712-1714):
    # '(11X,''Counted '',11(6X,''Corr'')/
    #  '' Series'',4X,''Segment '',11(''  Add  #'',I2)/
    #  '' -------- ---------'',11(''  --------''))' (J,J=1,11)
    if (8L %in% parts) {
      out_lines <- c(out_lines,
        paste0(strrep(" ", 11L), "Counted ",
               paste(rep(paste0(strrep(" ", 6L), "Corr"), n_best), collapse = "")),
        paste0(" Series    Segment ",
               paste(sprintf("  Add  #%2d", seq_len(n_best)), collapse = "")),
        paste0(" -------- ---------",
               paste(rep("  --------", n_best), collapse = ""))
      )
    }

    # N1 / A1 tally arrays (Fortran: N1(-MX:MX), A1(-MX:MX))
    # Use named list keyed by character(adj) for sparse storage
    n1 <- list()   # occurrence counts
    a1 <- list()   # correlation sums

    nsg <- 0L        # segment counter (NSG)
    ipr <- -32000L   # previous accepted segment index (IPR)

    # Segment loop: DO II=1,N,LAG (Fortran lines 1724-1781)
    ii <- 1L
    while (ii <= n) {

      # I = MIN(II, N-LS+1)  (Fortran line 1725)
      i <- min(ii, n - ls + 1L)
      if (i < 1L) {
        lsa <- min(ls, n)
        i   <- 1L
      } else {
        lsa <- ls
      }

      # JA, JZ: nominal calendar year range of this segment (lines 1735-1736)
      ja <- jyr + i - 1L
      jz <- ja  + lsa - 1L

      # Near-identical segment check (Fortran lines 1738-1742):
      # IF(I-IPR < LAG/2 .OR. I == IPR) → skip; if I-IPR > 0 print message
      if ((i - ipr) < (lag %/% 2L) || i == ipr) {
        if ((i - ipr) > 0L && 8L %in% parts)
          out_lines <- c(out_lines,
            sprintf(" %-8s%5d%5d  Lag from prior segment%4d years - insufficient",
                    sid, ja, jz, i - ipr))
        ii <- ii + lag
        next
      }
      ipr <- i
      nsg <- nsg + 1L

      # Slide segment along master (Fortran lines 1747-1751):
      # MS = IZ - LSA + 1 (last valid start position in master, 1-indexed)
      # DO J=IA,MS: CORREL(LSA, YMSM(J), ZSER(I), R(NR))
      seg  <- zf[i:min(i + lsa - 1L, n)]
      lsa  <- length(seg)    # actual length after clipping at series end
      ms   <- nm - lsa + 1L  # number of positions in master
      if (ms < 1L) { ii <- ii + lag; next }

      r_all   <- numeric(ms)
      adj_all <- integer(ms)
      for (jj in seq_len(ms)) {
        r_all[jj]   <- .cof_correl(mast[jj:(jj + lsa - 1L)], seg)
        # MXD(J) = IA + K - JA - 1  where K is 1-based position index
        adj_all[jj] <- iam + jj - 1L - ja
      }

      # Find 11 highest correlations sequentially (Fortran lines 1754-1764):
      # Initialise RMX=-9.99, find max, set that position to -9.99, repeat
      rmx <- rep(-9.99, n_best)
      mxd <- integer(n_best)
      rv  <- r_all    # work copy
      av  <- adj_all
      for (jj in seq_len(n_best)) {
        kmx      <- which.max(rv)
        rmx[jj]  <- rv[kmx]
        mxd[jj]  <- av[kmx]
        rv[kmx]  <- -9.99   # eliminate (Fortran line 1763)
      }

      # Tally N1 and A1 (Fortran lines 1774-1777)
      for (jj in seq_len(n_best)) {
        ak <- as.character(mxd[jj])
        n1[[ak]] <- (if (is.null(n1[[ak]])) 0L   else n1[[ak]]) + 1L
        a1[[ak]] <- (if (is.null(a1[[ak]])) 0.0  else a1[[ak]]) + rmx[jj]
      }

      # Data row (Fortran lines 1779-1780):
      # WRITE(IUW,'(1X,A,2I5,11(I6,A))') ID, JA, JZ, (MXD(J),LN(J),J=1,11)
      if (8L %in% parts) {
        pairs <- paste(sprintf("%6d%s", mxd, vapply(rmx, fmt_ln, character(1))),
                       collapse = "")
        out_lines <- c(out_lines,
          sprintf(" %-8s%5d%5d%s", sid, ja, jz, pairs))
      }

      ii <- ii + lag
    }   # end DO II loop; falls through to label 10

    # ---- Tally section (Fortran label 10, lines 1784-1811) ------------------
    # Build JP, N2 arrays: adjustments with N1 >= 3
    all_adjs <- as.integer(names(n1))
    all_cnts <- vapply(names(n1), function(k) n1[[k]], integer(1))
    all_rsums <- vapply(names(n1), function(k) a1[[k]], numeric(1))

    keep <- all_cnts >= 3L
    jp   <- all_adjs[keep]
    n2   <- all_cnts[keep]
    rsum <- all_rsums[keep]
    nk   <- length(jp)

    # RANKII(NK, N2, JP, 'D') — sort JP by N2 descending (Fortran line 1792)
    if (nk > 0L) {
      ord_d <- order(n2, decreasing = TRUE)
      jp    <- jp[ord_d];  n2 <- n2[ord_d];  rsum <- rsum[ord_d]
    }
    # A2(J) = A1(JP(J)) / N2(J)  (Fortran line 1793-1795)
    a2 <- if (nk > 0L) rsum / n2 else numeric(0)

    if (8L %in% parts) {
      if (nk > 0L) {
        # Line 1799: I7(NSG) + ' segments' + 38×'  -'
        out_lines <- c(out_lines,
          sprintf("%7d segments%s", nsg, strrep("  -", 38L)),
          # Line 1799 (cont): '/ Number of segments'
          " Number of segments",
          # Line 1800: '4X,8(''   Add No R_av'')'
          paste0("    ", paste(rep("   Add No R_av", 8L), collapse = "")))

        # Line 1801: (4X,(8(SP,I6,S,I3,F5.2)))  — rows of 8
        # SP+I6: mandatory-sign 6-char integer → f'{adj:+6d}'
        # S+I3:  no-sign 3-char integer        → f'{cnt:3d}'
        # F5.2:  5-char float                  → f'{r:5.2f}'
        for (rs in seq(1L, nk, by = 8L)) {
          re    <- min(rs + 7L, nk)
          items <- paste(sprintf("%+6d%3d%5.2f", jp[rs:re], n2[rs:re], a2[rs:re]),
                         collapse = "")
          out_lines <- c(out_lines, paste0("    ", items))
        }

        # Chronological order if NK > 2 (Fortran lines 1802-1805)
        if (nk > 2L) {
          # RANKII(NK, JP, N2, 'A') — sort JP ascending by JP value
          ord_a <- order(jp)
          jp_c  <- jp[ord_a];  n2_c <- n2[ord_a]
          out_lines <- c(out_lines,
            " Chronological order",
            # '4X,14(''   Add No'')'
            paste0("    ", paste(rep("   Add No", 14L), collapse = "")))
          # '(4X,(14(SP,I6,S,I3)))'
          for (rs in seq(1L, nk, by = 14L)) {
            re    <- min(rs + 13L, nk)
            items <- paste(sprintf("%+6d%3d", jp_c[rs:re], n2_c[rs:re]),
                           collapse = "")
            out_lines <- c(out_lines, paste0("    ", items))
          }
        }

      } else {
        # No pattern found (Fortran line 1808)
        out_lines <- c(out_lines,
          sprintf("%7d segments; no pattern found", nsg))
      }

      # Separator (Fortran line 1811): '(T2,18(''=''),11(2X,8(''='')))'
      # T2 = 1 leading space; 18×'='; 11×(2sp + 8×'=')
      out_lines <- c(out_lines,
        paste0(" ", strrep("=", 18L),
               paste(rep(paste0("  ", strrep("=", 8L)), 11L), collapse = "")))
    }

    # Summary entry
    best_adj <- if (nk > 0L) jp[1L] else NA_integer_
    n_supp   <- if (nk > 0L) n2[1L] else 0L
    mean_r_b <- if (nk > 0L) a2[1L] else NA_real_

    summary_list <- c(summary_list, list(data.frame(
      series            = sid,
      n                 = n,
      best_adj          = best_adj,
      n_segs_supporting = n_supp,
      mean_r_at_best    = mean_r_b,
      stringsAsFactors  = FALSE
    )))
  }   # end series loop (Fortran: GOTO 1)

  # Footer (Fortran lines 1698-1699): IF(NSR > 0) WRITE(I9,' undated series') NSR
  if (8L %in% parts && nsr > 0L)
    out_lines <- c(out_lines, sprintf("%9d undated series", nsr))

  # ---- Assemble result -------------------------------------------------------
  summary_df <- if (length(summary_list) > 0L)
    do.call(rbind, summary_list)
  else
    data.frame(series = character(0), n = integer(0),
               best_adj = integer(0), n_segs_supporting = integer(0),
               mean_r_at_best = numeric(0))

  result <- list(
    summary  = summary_df,
    crit     = cof_result$crit,
    options  = cof_result$options,
    output   = out_lines
  )

  if (!is.null(output_file)) {
    writeLines(out_lines, output_file)
    invisible(result)
  } else {
    result
  }
}


#' ASCII bar plot of one or more ring-width series
#'
#' @description
#' Produces the COFECHA-style ASCII bar plot (Fortran BARPL) for any set of
#' ring-width series. Each annual value is normalised and displayed as a
#' horizontal bar whose length encodes its decile rank and whose terminal
#' symbol encodes its SD class (A--Z positive, a--z negative, `@` for zero).
#' Verified against PUE benchmark: all bars match DPL reference output.
#'
#' Two layout modes:
#' - **\code{"page"}** (default) --- 8 columns x 50 rows per centennial block,
#'   identical to COFECHA Part 4.
#' - **\code{"column"}** --- single-column chronological listing, easier for comparing
#'   multiple short series side by side.
#'
#' @param rwl A dplR `rwl` data.frame **or** the list returned by \code{\link{dpl_cof}}
#'   (uses `$master`) **or** a named numeric vector (single series,
#'   names = years).
#' @param series Series to plot: a character vector of IDs or an integer
#'   vector of column positions in `rwl` (e.g. `1:7`). `NULL` (default)
#'   plots all columns. For a `cof_result` input, \code{"master"} is always
#'   available.
#' @param years Integer vector `c(first, last)` defining the display period.
#'   `NULL` (default) uses each series' full span. Normalisation always uses
#'   the full span by default (see `norm_scope`).
#' @param layout Character. \code{"page"} (default) or \code{"column"}.
#' @param norm_scope Character. \code{"full"} (default) --- normalise over the entire
#'   series span; \code{"period"} --- normalise within the display window only.
#' @param output_file Character or `NULL`. Write output to file. Default `NULL`.
#' @param overwrite Logical. `TRUE` (default) truncates; `FALSE` appends.
#' @param quiet Logical. Suppress console output. Default `FALSE`.
#' @param normalise Logical. Normalise to mean = 0, population SD = 1 before
#'   plotting. Default `TRUE`.
#'
#' @return Invisibly, a named list of character vectors holding the formatted
#'   lines: one element per series for \code{layout = "page"}; a single
#'   element \code{all} for \code{layout = "column"}, where series are tiled
#'   side by side. Use \code{writeLines()} on an element to reprint it.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#'
#' # ASCII bar plot of a single long series
#' dpl_barplot(rwl, series = "ACC026B")
#'
#' # Bar plot of the COFECHA master (Part 4 equivalent)
#' cof <- dpl_cof(rwl, parts = integer(0), verbose = FALSE)
#' dpl_barplot(cof)
#'
#' # Period display with normalisation over that window only
#' dpl_barplot(rwl, series = "ACC026B",
#'             years = c(1850, 2002), norm_scope = "period")
#'
#' # Write to file silently
#' dpl_barplot(cof, output_file = tempfile(fileext = ".txt"), quiet = TRUE)
#' }
#'
#' @seealso \code{\link{dpl_cof}}, \code{\link{plot_cof_barplot}}
#' @export
dpl_barplot <- function(rwl,
                          series      = NULL,
                          years       = NULL,
                          layout      = c("page", "column"),
                          norm_scope  = c("full", "period"),
                          output_file = NULL,
                          overwrite   = TRUE,
                          quiet       = FALSE,
                          normalise   = TRUE) {

  layout     <- match.arg(layout)
  norm_scope <- match.arg(norm_scope)

  # ---- 1. Extract series data -----------------------------------------------
  # ser_list holds the FULL span of every series (never clipped here).
  ser_list <- list()

  if (is.list(rwl) && !is.data.frame(rwl) && !is.null(rwl$master)) {
    yrs_m <- as.integer(names(rwl$master))
    ser_list[["master"]] <- setNames(as.numeric(rwl$master),
                                     as.character(yrs_m))
  } else if (is.data.frame(rwl)) {
    all_ids <- colnames(rwl)
    ids_use <- .resolve_series(series, all_ids)
    if (length(ids_use) == 0L)
      stop("No series found matching 'series' argument.")
    yrs_all <- as.integer(rownames(rwl))
    for (sid in ids_use) {
      col <- rwl[[sid]]
      ok  <- which(!is.na(col))
      if (length(ok) == 0L) next
      j1  <- ok[1L]; j2 <- ok[length(ok)]
      v   <- col[j1:j2];  v[is.na(v)] <- 0.0
      ser_list[[sid]] <- setNames(as.numeric(v), as.character(yrs_all[j1:j2]))
    }
  } else if (is.numeric(rwl) && !is.null(names(rwl))) {
    sid <- if (is.character(series)) series[1L] else "series"
    ser_list[[sid]] <- rwl
  } else {
    stop("'rwl' must be a dplR rwl data.frame, a cof_result list, ",
         "or a named numeric vector.")
  }

  if (length(ser_list) == 0L) stop("No valid series to plot.")

  # ---- 2. Validate period ---------------------------------------------------
  if (!is.null(years)) {
    if (length(years) != 2L || years[1L] > years[2L])
      stop("'years' must be c(first_year, last_year) with first_year <= last_year.")
  }

  # ---- 3. Generate ASCII lines per series -----------------------------------
  # Pipeline per series:
  #   (a) normalise + compute decile cut-points over norm_scope
  #       norm_scope = "full"   -> entire series span  (default)
  #       norm_scope = "period" -> display window only (years argument)
  #   (b) bar positions use the mean/SD from (a) applied to the FULL span,
  #       so a year outside the display window would plot correctly if shown
  #   (c) display is clipped to 'years' — bars outside the window are blank

  all_lines <- vector("list", length(ser_list))
  names(all_lines) <- names(ser_list)

  for (sid in names(ser_list)) {

    v_full   <- ser_list[[sid]]
    yrs_full <- as.integer(names(v_full))

    # Subset used for normalisation and decile ranking
    if (!is.null(years) && norm_scope == "period") {
      keep_norm <- yrs_full >= years[1L] & yrs_full <= years[2L]
      v_norm    <- v_full[keep_norm]
    } else {
      v_norm    <- v_full
    }
    if (length(v_norm) == 0L) {
      all_lines[[sid]] <- NULL
      next
    }

    # Normalise v_norm; apply the same (mean, SD) to the full span
    if (normalise) {
      nr     <- .cof_normts(as.numeric(v_norm), k = 0L)
      xm_use <- nr$mean;  sd_use <- nr$sd
      if (sd_use > 0) {
        Y_full <- setNames((as.numeric(v_full) - xm_use) / sd_use,
                           names(v_full))
        y_dec  <- (as.numeric(v_norm) - xm_use) / sd_use
      } else {
        Y_full <- setNames(as.numeric(v_full) - xm_use, names(v_full))
        y_dec  <- as.numeric(v_norm) - xm_use
      }
    } else {
      Y_full <- v_full
      y_dec  <- as.numeric(v_norm)
    }

    # Decile cut-points from the norm-scope values (Fortran RANKRI / NINT)
    Z <- .cof_barpl_cuts(as.numeric(y_dec))

    # Display window for this series
    disp_jyr <- if (!is.null(years)) max(years[1L], min(yrs_full)) else min(yrs_full)
    disp_lyr <- if (!is.null(years)) min(years[2L], max(yrs_full)) else max(yrs_full)

    # make_car: 16-char cell (Fortran BARPL convention, shared helper);
    # blank outside the display window and for NA values
    make_car_fn <- local({
      Y_ <- Y_full;  Z_ <- Z
      dj <- as.integer(disp_jyr);  dl <- as.integer(disp_lyr)
      function(yr) {
        if (yr < dj || yr > dl) return(strrep(" ", 16L))
        yn <- Y_[as.character(yr)]
        if (length(yn) == 0L || is.na(yn)) return(strrep(" ", 16L))
        .cof_barpl_car(yr, as.numeric(yn), Z_)
      }
    })

    # Header notes
    norm_note <- if (normalise) {
      if (!is.null(years) && norm_scope == "period")
        sprintf("  [norm: period %d-%d]", years[1L], years[2L])
      else
        sprintf("  [norm: full span %d-%d]", min(yrs_full), max(yrs_full))
    } else "  [raw values]"
    disp_note <- if (!is.null(years))
      sprintf("  [display: %d-%d]", disp_jyr, disp_lyr) else ""

    all_lines[[sid]] <- list(
      make_car = make_car_fn,
      disp_jyr = as.integer(disp_jyr),
      disp_lyr = as.integer(disp_lyr),
      full_jyr = as.integer(min(yrs_full)),
      full_lyr = as.integer(max(yrs_full)),
      full_n   = as.integer(length(v_full)),
      norm_note = norm_note,
      disp_note = disp_note
    )
  }

  all_lines <- Filter(Negate(is.null), all_lines)
  if (length(all_lines) == 0L)
    stop("No series have data in the requested 'years' range.")

  sids <- names(all_lines)   # ordered series IDs

  # ---- 4. Assemble output lines per layout ----------------------------------

  out_lines  <- character(0)
  per_series <- list()

  if (layout == "page") {
    # PAGE layout: one series per block, 8 columns x 50 rows -----------------
    for (sid in sids) {
      s       <- all_lines[[sid]]
      make_c  <- s$make_car
      djyr    <- s$disp_jyr;  dlyr <- s$disp_lyr

      page_hdr <- c(
        "",
        strrep("=", 132L),
        sprintf(" Series: %-8s   %d to %d   (%d yr)%s%s",
                sid, s$full_jyr, s$full_lyr, s$full_n,
                s$norm_note, s$disp_note),
        strrep("-", 132L)
      )
      # Same 400-year page engine as COFECHA Part 4 (see .cof_barpl_pages)
      body   <- .cof_barpl_pages(make_c, djyr, dlyr, page_hdr)
      footer <- c(paste0(" ", strrep("-", 131L)),
                  sprintf("%8d years displayed", dlyr - djyr + 1L))
      per_series[[sid]] <- c(body, footer)
      out_lines <- c(out_lines, body, footer)
    }

  } else {
    # COLUMN layout: multiple series side-by-side, tiled across pages ---------
    #
    # Layout dimensions:
    #   YEAR_W = 6   (" YYYY  ")
    #   BAR_W  = 12  (bar right-padded to 11 chars + 1 space separator)
    #   PAGE_W = 132
    #   COLS_PER_PAGE = floor((PAGE_W - YEAR_W) / BAR_W) = 10
    #
    # Year axis: union of all display windows; all years in the full range
    # shown even if blank for some series (keeps a clean calendar grid).

    YEAR_W        <- 6L
    BAR_W         <- 12L
    COLS_PER_PAGE <- (132L - YEAR_W) %/% BAR_W   # = 10

    # Global year range: span all series' display windows
    all_djyr <- vapply(all_lines, `[[`, integer(1), "disp_jyr")
    all_dlyr <- vapply(all_lines, `[[`, integer(1), "disp_lyr")
    g_jyr    <- min(all_djyr)
    g_lyr    <- max(all_dlyr)
    all_yrs  <- g_jyr:g_lyr

    # Tile series into pages of COLS_PER_PAGE columns each
    n_ser    <- length(sids)
    n_pages  <- ceiling(n_ser / COLS_PER_PAGE)

    for (pg in seq_len(n_pages)) {
      idx_first <- (pg - 1L) * COLS_PER_PAGE + 1L
      idx_last  <- min(pg * COLS_PER_PAGE, n_ser)
      page_sids <- sids[idx_first:idx_last]
      n_cols    <- length(page_sids)

      # Column header: year label + one ID per series column
      id_hdr <- vapply(page_sids, function(s)
        formatC(s, width = BAR_W - 1L, flag = "-"), character(1))
      hdr_line <- paste0(strrep(" ", YEAR_W),
                         paste(id_hdr, collapse = " "))

      # Span line: show each series' own span under its ID
      span_hdr <- vapply(page_sids, function(s) {
        info <- all_lines[[s]]
        sp   <- sprintf("%d-%d", info$full_jyr, info$full_lyr)
        formatC(sp, width = BAR_W - 1L, flag = "-")
      }, character(1))
      span_line <- paste0(strrep(" ", YEAR_W),
                          paste(span_hdr, collapse = " "))

      # Norm note: one per series column (abbreviated)
      norm_hdr <- vapply(page_sids, function(s) {
        nt <- all_lines[[s]]$norm_note
        # Trim to fit in BAR_W-1 chars
        nt_short <- sub("^  \\[", "[", nt)
        formatC(substr(nt_short, 1L, BAR_W - 1L), width = BAR_W - 1L,
                flag = "-")
      }, character(1))
      norm_line <- paste0(strrep(" ", YEAR_W),
                          paste(norm_hdr, collapse = " "))

      page_hdr <- c(
        "",
        strrep("=", 132L),
        if (n_pages > 1L)
          sprintf(" Column layout -- page %d of %d  (years %d to %d)",
                  pg, n_pages, g_jyr, g_lyr)
        else
          sprintf(" Column layout  (years %d to %d)", g_jyr, g_lyr),
        strrep("-", 132L),
        hdr_line,
        span_line,
        norm_line,
        strrep("-", min(132L, YEAR_W + n_cols * BAR_W))
      )
      out_lines <- c(out_lines, page_hdr)

      # Data rows: one row per year
      for (idx in seq_along(all_yrs)) {
        yr <- all_yrs[idx]

        # Year label: right-justified in 4 chars + 2 spaces
        yr_lbl <- sprintf("%4d  ", yr)

        # Bar cell for each series: strip the year prefix, right-pad to BAR_W-1
        bar_cells <- vapply(page_sids, function(s) {
          car <- all_lines[[s]]$make_car(yr)
          # car is I5(yr) + bar, padded to 16. We want just the bar symbol part
          # (strip the year prefix already in yr_lbl, keep only the bar chars).
          bar_only <- substr(car, 6L, 16L)   # chars 6-16 = the bar (max 11 chars)
          bar_only <- trimws(bar_only, which = "right")
          formatC(bar_only, width = BAR_W - 1L, flag = "-")
        }, character(1))

        out_lines <- c(out_lines,
                       paste0(yr_lbl, paste(bar_cells, collapse = " ")))

        # Blank line every 10 years (between groups, not after last)
        if (idx %% 10L == 0L && idx < length(all_yrs)) out_lines <- c(out_lines, "")
      }

      footer <- c(
        "",
        strrep("-", min(132L, YEAR_W + n_cols * BAR_W)),
        sprintf("%4d years  |  %d series shown on this page", length(all_yrs), n_cols)
      )
      out_lines <- c(out_lines, footer)
    }
  }

  # ---- 5. Output ------------------------------------------------------------
  flat <- out_lines

  if (!quiet)
    cat(paste(flat, collapse = "\n"), "\n", sep = "")

  if (!is.null(output_file)) {
    mode <- if (overwrite) "w" else "a"
    con  <- file(output_file, open = mode, encoding = "native.enc")
    on.exit(close(con), add = TRUE)
    writeLines(flat, con)
    if (!quiet)
      message(sprintf("ASCII bar plot written to: %s", output_file))
  }

  # Page layout: one element per series.  Column layout tiles series side by
  # side, so the lines cannot be split per series: return them under "all".
  if (layout == "page") invisible(per_series) else invisible(list(all = flat))
}



#' Graphical bar plot of the COFECHA master dating series
#'
#' @description
#' Produces a base-R graphical bar plot of the master dating series from a
#' completed \code{\link{dpl_cof}} run. Positive departures are coloured `col_pos`
#' and negative departures `col_neg`. An optional secondary axis overlays
#' sample depth as a grey line. The layout mirrors the information content of
#' Holmes' Part 4 ASCII bar plot but is designed for screen display and
#' publication export.
#'
#' @param cof_result The list returned by \code{\link{dpl_cof}}. Must contain `$master`
#'   and (if `add_sample_depth = TRUE`) `$sample_depth`.
#' @param add_sample_depth Logical. If `TRUE` (default), overlays sample depth
#'   on a secondary right axis as a grey line.
#' @param col_pos Character. Colour for positive departures. Default
#'   \code{"steelblue"}.
#' @param col_neg Character. Colour for negative departures. Default \code{"tomato"}.
#' @param main Character or `NULL`. Plot title. `NULL` (default) auto-constructs
#'   a title from the master span.
#' @param ... Additional arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @return Invisibly returns the normalised master series plotted.
#'
#' @examples
#' \dontrun{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#' cof <- dpl_cof(rwl, parts = integer(0), verbose = FALSE)
#'
#' # Interactive / screen plot
#' plot_cof_barplot(cof,
#'   main = "Fitzroya cupressoides --- Cerro Mirador, Alerce Costero")
#'
#' # Publication export
#' pdf(tempfile(fileext = ".pdf"), width = 14, height = 4)
#' plot_cof_barplot(cof,
#'   main = "Cerro Mirador master dating series (Barichivich 2005)",
#'   col_pos = "steelblue", col_neg = "tomato")
#' dev.off()
#' }
#'
#' @seealso \code{\link{dpl_cof}}, \code{\link{dpl_barplot}}
#' @export
plot_cof_barplot <- function(cof_result,
                              add_sample_depth = TRUE,
                              col_pos  = "steelblue",
                              col_neg  = "tomato",
                              main     = NULL,
                              ...) {

  if (is.null(cof_result$master))
    stop("'cof_result' must be the output of dpl_cof().")

  master <- as.numeric(cof_result$master)
  yrs    <- as.integer(names(cof_result$master))

  # Normalise to mean = 0, sd = 1
  nr  <- .cof_normts(master, k = 0L)
  y   <- nr$z

  cols <- ifelse(y >= 0, col_pos, col_neg)
  ttl  <- if (!is.null(main)) main else
    sprintf("Master Dating Series  %d \u2013 %d", min(yrs), max(yrs))

  if (add_sample_depth && !is.null(cof_result$sample_depth)) {
    nd  <- as.integer(cof_result$sample_depth)
    op  <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))
    graphics::par(mar = c(4, 4, 3, 4))
  }

  graphics::barplot(y, names.arg = yrs, col = cols, border = NA,
                    ylab = "Normalised index", xlab = "Year",
                    main = ttl, las = 2L, ...)
  graphics::abline(h = 0, col = "black", lwd = 0.5)

  if (add_sample_depth && !is.null(cof_result$sample_depth)) {
    nd_sc <- nd / max(nd) * graphics::par("usr")[4L] * 0.35
    graphics::par(new = TRUE)
    graphics::plot(yrs, nd_sc, type = "l", col = "grey50", lwd = 1.5,
                   axes = FALSE, xlab = "", ylab = "",
                   ylim = c(0, graphics::par("usr")[4L]))
    graphics::axis(4, at = pretty(nd_sc),
                   labels = round(pretty(nd_sc) / graphics::par("usr")[4L]
                                  * max(nd) / 0.35),
                   col.axis = "grey50", col = "grey50", las = 1L)
    graphics::mtext("Sample depth (n series)", side = 4L, line = 2.5,
                    col = "grey50", cex = 0.85)
  }

  invisible(y)
}
