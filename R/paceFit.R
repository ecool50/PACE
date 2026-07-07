#' Fit a PACE model to a SpatialExperiment
#'
#' `paceFit()` runs the full PACE pipeline on a
#' [SpatialExperiment::SpatialExperiment]: it builds the biological and technical
#' neighbourhood kernels from the spatial coordinates and cell-type labels, fits
#' the hierarchical negative binomial mixed model by streaming penalised
#' quasi-likelihood (with the per-cell contamination correction), shrinks the
#' neighbour slopes with multivariate adaptive shrinkage, and computes the
#' variance decomposition and per-pair driver tables.
#'
#' Cell-type labels, spatial coordinates, and (optionally) an image/sample grouping
#' and a condition are read from the object following SpatialExperiment
#' conventions. The defaults reproduce the manuscript recipe (30 um biological and
#' 5 um technical bandwidths, per-cell homotypic-core contamination correction).
#'
#' @param object A [SpatialExperiment::SpatialExperiment] with a counts assay and
#'   two-dimensional spatial coordinates.
#' @param celltype_col Name of the [SummarizedExperiment::colData] column holding
#'   the discrete cell-type annotation.
#' @param image_col Optional name of the colData column grouping cells into images
#'   or samples. If `NULL` (default) the whole object is treated as one section.
#' @param condition_col Optional name of a colData column holding a binary
#'   condition (for example treatment response); enables the responder spatial
#'   block for multi-sample cohorts.
#' @param assay_name Name of the counts assay to use (default `"counts"`).
#' @param h_bio,h_tech Biological and technical kernel bandwidths in micrometres
#'   (defaults 30 and 5).
#' @param contamination Contamination correction: `"percell_hc"` (default,
#'   per-cell homotypic-core ambient term) or `"none"`.
#' @param dispersion Negative binomial parameterisation, `"nb1"` (default) or
#'   `"nb2"`.
#' @param kernel_per_image If `TRUE`, neighbourhoods are built within each image
#'   (use when images are separate samples/patients). Default `FALSE`.
#' @param image_re Second image-level random-effect block: one of `"none"`,
#'   `"intercept"`, `"slopes"`, `"condition_slopes"`.
#' @param resp_term For condition cohorts, the responder interaction term prefix
#'   used when shrinking and decomposing (for example `"ResponderPD"`). If `NULL`
#'   it is derived from `condition_col`.
#' @param pairs Optional list of focal-neighbour cell-type pairs for the driver
#'   tables; if `NULL`, all pairs are scored.
#' @param verbose Whether to print progress.
#' @param ... Further arguments passed to the streaming fitter (for example
#'   `n_iter`, `threads`, `chunk_size`).
#'
#' @return A [PACEFit] object.
#'
#' @examples
#' \dontrun{
#' library(SpatialExperiment)
#' fit <- paceFit(spe, celltype_col = "cellType")
#' head(neighbourSlopes(fit))
#' }
#'
#' @rdname paceFit
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom Matrix rowSums
#' @export
setMethod(
  "paceFit", "SpatialExperiment",
  function(object, celltype_col, image_col = NULL, condition_col = NULL,
           assay_name = "counts", h_bio = 30, h_tech = 5,
           contamination = c("percell_hc", "none"),
           dispersion = c("nb1", "nb2"),
           kernel_per_image = FALSE,
           image_re = c("none", "intercept", "slopes", "condition_slopes"),
           resp_term = NULL, pairs = NULL, verbose = TRUE, ...) {

    contamination <- match.arg(contamination)
    dispersion    <- match.arg(dispersion)
    image_re      <- match.arg(image_re)

    ## ---- pull counts, coordinates, and cell metadata from the SPE ----
    if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
      stop("assay '", assay_name, "' not found in the object.", call. = FALSE)
    }
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    if (!celltype_col %in% colnames(cd)) {
      stop("celltype_col '", celltype_col, "' not found in colData.", call. = FALSE)
    }

    ## counts are genes x cells in an SPE; PACE wants cells x genes.
    Y <- t(as.matrix(SummarizedExperiment::assay(object, assay_name)))

    coords <- SpatialExperiment::spatialCoords(object)
    if (ncol(coords) < 2L) stop("spatialCoords must have two columns.", call. = FALSE)
    df <- cd
    df[["x"]] <- coords[, 1]
    df[["y"]] <- coords[, 2]

    ## single-section objects get a constant image grouping.
    if (is.null(image_col)) {
      df[[".image"]] <- factor("all")
      image_col <- ".image"
    }

    ## derive the responder term prefix for condition cohorts.
    if (!is.null(condition_col) && is.null(resp_term)) {
      lev <- levels(factor(df[[condition_col]]))
      resp_term <- paste0("Responder", lev[length(lev)])
    }

    ## ---- 1. streaming PQL fit (the method) ----
    res <- pace_fit_streaming(
      Y, df, celltype_col = celltype_col, image_col = image_col,
      coord_cols = c("x", "y"), h_bio = h_bio, h_tech = h_tech,
      contamination = contamination, dispersion = dispersion,
      condition_col = condition_col, kernel_per_image = kernel_per_image,
      image_re = image_re, verbose = verbose, ...)
    types <- res$types

    ## ---- 2. mash shrinkage of the neighbour slopes ----
    slopes <- pace_shrink(res$fit, types, resp_term = resp_term)

    ## ---- 3. variance decomposition (link-scale blocks + observed single-frame) ----
    dec <- pace_decompose(res$fit, res$df, res$Y, types, res$X_fixed,
                          resp_term = resp_term)
    single_frame <- .pace_single_frame(res, dec, celltype_col)
    decomposition <- list(perGene = single_frame, blocks = dec)

    ## ---- 4. per-pair driver tables ----
    drivers <- tryCatch(
      pace_top_drivers(res$fit, slopes, dec, types, pairs = pairs),
      error = function(e) {
        warning("topDrivers could not be computed: ", conditionMessage(e),
                call. = FALSE)
        list()
      })

    methods::new(
      "PACEFit",
      fit                   = res$fit,
      neighbourSlopes       = as.data.frame(slopes),
      varianceDecomposition = decomposition,
      topDrivers            = as.list(drivers),
      cellTypes             = as.character(types),
      params = list(h_bio = h_bio, h_tech = h_tech,
                    contamination = contamination, dispersion = dispersion,
                    celltype_col = celltype_col, image_col = image_col,
                    condition_col = condition_col, resp_term = resp_term,
                    kernel_per_image = kernel_per_image, image_re = image_re))
  })

## Observed single-frame (log1p CP10k) decomposition from the fitted block
## proportions; the manuscript's headline per-gene frame. Internal.
.pace_single_frame <- function(res, dec, celltype_col) {
  block <- if (!is.null(dec$gene_focal_5block)) dec$gene_focal_5block
           else dec$gene_focal_4block
  if (is.null(block)) return(NULL)
  nCount <- as.numeric(Matrix::rowSums(res$Y))
  celltype <- res$df[[celltype_col]]
  tryCatch(
    single_frame_decomp_obs(res$Y, celltype, nCount, block),
    error = function(e) block)
}
