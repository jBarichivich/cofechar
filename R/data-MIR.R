#' Fitzroya cupressoides ring-width data --- Cerro Mirador, Alerce Costero
#'
#' @description
#' Ring-width measurements from 36 \emph{Fitzroya cupressoides} cores collected in semi-open stands at the summit of Cerro Mirador, Alerce
#' Costero National Park, Cordillera Pelada, southern Chile (ca. 40S).
#' The file is stored as \code{inst/extdata/CL-MIR.rwl} in the Holmes decadal
#' format (8-character label, 4-digit decade year, up to 10 space-separated
#' integer values per line, sentinel \code{-9999} marking end of series).
#' Values are in 1/100 mm.
#'
#' @details
#' The chronology spans 1406--2002, making it one of the longest tree-ring
#' records from South American temperate rainforests. Series identifiers follow
#' the convention: tree code (ACC0XX / ACC2X) + core letter (A/B/C or TA/TB).
#'
#' The file is accessed via \code{system.file("extdata", "CL-MIR.rwl",
#' package = "cofechar")} and read with \code{\link{dpl_read_dec}}.
#'
#' @format A plain-text file in Holmes decadal format:
#' \describe{
#'   \item{Series}{36 cores from individual trees (IDs: ACC013B to ACC27TB)}
#'   \item{Span}{1406--2002 (maximum series length 597 years)}
#'   \item{Units}{1/100 mm (integer; divide by 100 for mm)}
#'   \item{Sentinel}{\code{-9999} marks end of each series}
#'   \item{Missing}{No within-series missing values; \code{-9999} is end-of-series only}
#' }
#'
#' @source
#' Barichivich, J. (2005). Dendroclimatologia de la Cordillera Pelada, X
#' Region, Chile. Undergraduate thesis (Licenciatura en Ciencias Forestales),
#' Universidad Austral de Chile, Valdivia.
#'
#' Site: Cerro Mirador, Alerce Costero National Park, Los Rios Region, Chile.
#' Approximate coordinates: 40 deg 10' S, 73 deg 50' W, ~1000 m a.s.l.
#' Collector: Jonathan Barichivich.
#'
#' @name CL-MIR
#' @examples
#' mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
#' rwl <- dpl_read_dec(mir_file, label_length = NULL, stop_val = -9999L,
#'                     unit = "0.001mm")
#' head(rwl[, 1:4])
#' range(as.integer(rownames(rwl)))   # 1406 to 2002
NULL
