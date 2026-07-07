#' PACE: Proximity-Associated Changes in Expression
#'
#' PACE quantifies how a cell's gene expression changes with proximity to specific
#' neighbouring cell types in imaging-based spatial transcriptomics. It fits
#' hierarchical negative binomial mixed models with partial pooling across cell
#' types, corrects for transcript contamination between adjacent cells with a
#' per-cell ambient term, and decomposes expression variance into cell-type
#' identity, spatial cell state, contamination, and residual components.
#'
#' The user-facing entry point is [paceFit()], which takes a
#' [SpatialExperiment::SpatialExperiment] and returns a [PACEFit] object. The
#' results are read out with [neighbourSlopes()], [varianceDecomposition()], and
#' [topDrivers()].
#'
#' @keywords internal
#' @name PACE-package
#' @aliases PACE
#'
#' @import methods
#' @import dplyr
#' @import tibble
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom mashr mash mash_set_data cov_canonical cov_ed cov_pca get_significant_results
#' @importFrom ashr get_lfsr get_pm get_psd
#' @importFrom stats var sd median quantile pnorm qnorm optimize coef predict
#'   as.formula model.matrix setNames rnorm
#' @importFrom methods new validObject is
"_PACKAGE"
