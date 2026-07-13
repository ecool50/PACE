#' PACEFit: a fitted PACE model
#'
#' An S4 container for the result of a PACE fit. It holds the raw penalised
#' quasi-likelihood solver state together with the three interpretable output
#' layers: shrunken per-(gene, focal, neighbour) proximity slopes, the per-gene
#' variance decomposition, and the per-pair driver tables. Construct one with
#' [paceFit()]; read it out with [neighbourSlopes()], [varianceDecomposition()],
#' and [topDrivers()].
#'
#' @slot fit The raw streaming PQL solver state (fixed and random effects,
#'   per-cell contamination loadings, gene-wise overdispersions, variance
#'   components).
#' @slot neighbourSlopes A data frame of shrunken proximity coefficients, one row
#'   per (gene, focal cell type, neighbour cell type[, condition term]), with the
#'   raw estimate, shrunken estimate, and local false sign rate (lfsr).
#' @slot varianceDecomposition A list of per-gene, per-focal variance
#'   decomposition tables (cell-type identity, spatial cell state, contamination,
#'   residual; and, for condition cohorts, a responder spatial block).
#' @slot topDrivers A list of per-pair driver-score tables ranking the genes that
#'   mediate each focal-neighbour spatial relationship.
#' @slot cellTypes The cell types modelled, in fitting order.
#' @slot params The fitting parameters actually used (bandwidths, contamination
#'   and dispersion settings, condition column, and so on).
#' @slot context The working frame and fixed-effect design retained from the fit
#'   so that the downstream stages ([paceShrink()], [paceDecompose()],
#'   [paceDrivers()]) can run without refitting.
#'
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' fit
#' cellTypes <- fit@cellTypes
#' @name PACEFit
#' @rdname PACEFit-class
#' @exportClass PACEFit
setClass(
  "PACEFit",
  slots = c(
    fit                   = "list",
    neighbourSlopes       = "data.frame",
    varianceDecomposition = "list",
    topDrivers            = "list",
    cellTypes             = "character",
    params                = "list",
    context               = "list"
  ),
  prototype = list(
    neighbourSlopes       = data.frame(),
    varianceDecomposition = list(),
    topDrivers            = list(),
    context               = list()
  )
)
