#' cofechar: COFECHA and EDT --- R Port of the DPL Dendrochronology Program Library
#'
#' @description
#' A faithful R port of two core modules from Holmes' Dendrochronology Program
#' Library (DPL): **EDT** (ring-measurement editing) and **COFECHA** (quality
#' control and crossdating). All internal statistical routines --- cubic smoothing
#' spline, Burg autoregressive modelling, log-transform, normalisation --- are
#' translated directly from the Fortran source (Holmes 1982--1994) to ensure
#' numerical equivalence with DPL reference output.
#'
#' @details
#' ## EDT functions
#' Reading, writing, editing, and merging ring-measurement files:
#' \describe{
#'   \item{\code{\link{dpl_read}}}{Read a DPL file (compact `~` or Tucson `.rwl`) into R}
#'   \item{\code{\link{dpl_read_dec}}}{Read a decadal-format file with any label length}
#'   \item{\code{\link{dpl_write}}}{Write a dplR `rwl` or internal series list to a DPL file}
#'   \item{\code{\link{dpl_merge}}}{Merge two or more `rwl` objects onto a common year axis}
#'   \item{\code{\link{dpl_edt}}}{Apply a list of named edits to a `rwl` or series list}
#'   \item{\code{\link{dpl_display}}}{Print a ~50-year inspection window around a given year}
#'   \item{\code{\link{dpl_edit_file}}}{Convenience wrapper: read, edit, write in one call}
#' }
#'
#' ## COF functions
#' COFECHA quality control, diagnostics and floating-series dating:
#' \describe{
#'   \item{\code{\link{dpl_cof}}}{Run COFECHA on a dated `rwl`; return structured results}
#'   \item{\code{\link{dpl_cof_diag}}}{Print Part 5/6 diagnostics for selected series}
#'   \item{\code{\link{dpl_dateme}}}{Date floating (undated) series against a COFECHA master}
#'   \item{\code{\link{dpl_barplot}}}{ASCII bar plot of any `rwl` series or the COF master}
#'   \item{\code{\link{plot_cof_barplot}}}{Graphical bar plot of the master dating series}
#' }
#'
#' ## dplR integration
#' Every public function accepts dplR `rwl` data.frames and returns them by
#' default, so results slot directly into \code{\link[dplR]{detrend}}, \code{\link[dplR]{chron}},
#' \code{\link[dplR]{rwl.stats}}, and related functions without further conversion.
#'
#' ## Supported file formats
#'
#' | Format        | Description                                             |
#' |---------------|---------------------------------------------------------|
#' | Tucson `.rwl` | ITRDB standard; decades-per-line, F6.2, sentinel 9.99  |
#' | DPL compact   | Holmes `~` format; fixed-width integer, auto-scaled     |
#' | Holmes decadal | 8-char label + 4-digit decade year + integer values    |
#'
#' Format is auto-detected from file content when `format = "auto"`.
#'
#' ## Example data
#'
#' The package includes \code{CL-MIR.rwl}: ring-width measurements from 36
#' \emph{Fitzroya cupressoides} cores collected at Cerro Mirador,
#' Alerce Costero National Park, Chile (Barichivich 2005). The chronology
#' spans 1406--2002 and serves as the example dataset throughout the
#' documentation. Access it with:
#'
#' \preformatted{
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#' }
#'
#' ## References
#'
#' Holmes, R. L. (1983). Computer-assisted quality control in tree-ring dating
#' and measurement. *Tree-Ring Bulletin* 43:69--78.
#'
#' Holmes, R. L. (1994). *Dendrochronology Program Library*. Laboratory of
#' Tree-Ring Research, University of Arizona, Tucson.
#'
#' @seealso \code{\link[dplR]{read.rwl}}, \code{\link[dplR]{detrend}}, \code{\link[dplR]{chron}}
#'
#' @keywords internal
"_PACKAGE"
