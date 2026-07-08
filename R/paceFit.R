## ---------------------------------------------------------------------------
## Shared SPE -> (counts, cell frame) extraction used by the fit entry points.
## ---------------------------------------------------------------------------
.pace_prepare <- function(object, celltype_col, image_col, condition_col,
                          assay_name, resp_term) {
  if (!assay_name %in% SummarizedExperiment::assayNames(object))
    stop("assay '", assay_name, "' not found in the object.", call. = FALSE)
  cd <- as.data.frame(SummarizedExperiment::colData(object))
  if (!celltype_col %in% colnames(cd))
    stop("celltype_col '", celltype_col, "' not found in colData.", call. = FALSE)

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
  list(Y = Y, df = df, image_col = image_col, resp_term = resp_term)
}

## Observed single-frame (log1p CP10k) decomposition from the fitted block
## proportions; the manuscript's headline per-gene frame. Internal.
.pace_single_frame <- function(Y, df, celltype_col, dec) {
  block <- if (!is.null(dec$gene_focal_5block)) dec$gene_focal_5block
           else dec$gene_focal_4block
  if (is.null(block)) return(NULL)
  nCount   <- as.numeric(Matrix::rowSums(Y))
  celltype <- df[[celltype_col]]
  tryCatch(single_frame_decomp_obs(Y, celltype, nCount, block),
           error = function(e) block)
}

## ===========================================================================
## Stage 1: fit the model
## ===========================================================================

#' Fit the PACE model (stage 1 of the pipeline)
#'
#' `paceModel()` builds the biological and technical neighbourhood kernels from
#' the spatial coordinates and cell-type labels and fits the hierarchical
#' negative binomial mixed model by streaming penalised quasi-likelihood (with
#' the per-cell contamination correction). It returns a [PACEFit] carrying only
#' the fitted model; the reporting layers are added by [paceShrink()],
#' [paceDecompose()], and [paceDrivers()] (or all at once by [paceFit()]).
#'
#' @param object A [SpatialExperiment::SpatialExperiment] with a counts assay and
#'   two-dimensional spatial coordinates.
#' @param celltype_col colData column with the discrete cell-type annotation.
#' @param image_col Optional colData column grouping cells into images/samples.
#'   `NULL` (default) treats the object as one section.
#' @param condition_col Optional colData column with a binary condition; enables
#'   the responder spatial block for multi-sample cohorts.
#' @param assay_name Counts assay name (default `"counts"`).
#' @param h_bio,h_tech Biological and technical kernel bandwidths in micrometres
#'   (defaults 30 and 5).
#' @param contamination `"percell_hc"` (default) or `"none"`.
#' @param dispersion `"nb1"` (default) or `"nb2"`.
#' @param kernel_per_image If `TRUE`, neighbourhoods are built within each image.
#' @param image_re Second image-level random-effect block: `"none"`,
#'   `"intercept"`, `"slopes"`, or `"condition_slopes"`.
#' @param resp_term Responder interaction term prefix for condition cohorts; if
#'   `NULL` it is derived from `condition_col`.
#' @param verbose Whether to print progress.
#' @param ... Further arguments passed to the streaming fitter (`n_iter`,
#'   `threads`, `drop_sparse_neff`, `within_image`, `edge_correct`,
#'   `data_informed_tau`, `tau_shrinkage`, ...).
#' @return A [PACEFit] with the fitted model (reporting layers empty).
#' @examples
#' \dontrun{
#' fit <- paceModel(spe, celltype_col = "cellType") |>
#'          paceShrink() |> paceDecompose(spe) |> paceDrivers()
#' }
#' @rdname paceModel
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom Matrix rowSums
#' @export
setMethod(
  "paceModel", "SpatialExperiment",
  function(object, celltype_col, image_col = NULL, condition_col = NULL,
           assay_name = "counts", h_bio = 30, h_tech = 5,
           contamination = c("percell_hc", "none"),
           dispersion = c("nb1", "nb2"), kernel_per_image = FALSE,
           image_re = c("none", "intercept", "slopes", "condition_slopes"),
           resp_term = NULL, verbose = TRUE, ...) {
    contamination <- match.arg(contamination)
    dispersion    <- match.arg(dispersion)
    image_re      <- match.arg(image_re)

    prep <- .pace_prepare(object, celltype_col, image_col, condition_col,
                          assay_name, resp_term)
    res <- pace_fit_streaming(
      prep$Y, prep$df, celltype_col = celltype_col, image_col = prep$image_col,
      coord_cols = c("x", "y"), h_bio = h_bio, h_tech = h_tech,
      contamination = contamination, dispersion = dispersion,
      condition_col = condition_col, kernel_per_image = kernel_per_image,
      image_re = image_re, verbose = verbose, ...)

    methods::new(
      "PACEFit",
      fit       = res$fit,
      cellTypes = as.character(res$types),
      context   = list(df = res$df, X_fixed = res$X_fixed,
                       genes = colnames(res$fit$U)),
      params = list(h_bio = h_bio, h_tech = h_tech,
                    contamination = contamination, dispersion = dispersion,
                    celltype_col = celltype_col, image_col = prep$image_col,
                    condition_col = condition_col, resp_term = prep$resp_term,
                    kernel_per_image = kernel_per_image, image_re = image_re,
                    assay_name = assay_name))
  })

## ===========================================================================
## Stage 2: shrink the neighbour slopes
## ===========================================================================

#' Shrink the neighbour slopes (stage 2)
#'
#' Stabilises the fitted per-(gene, focal, neighbour) proximity slopes with
#' multivariate adaptive shrinkage, populating [neighbourSlopes()].
#'
#' @param object A [PACEFit] from [paceModel()].
#' @param ... Unused.
#' @return The `PACEFit` with the shrunken neighbour slopes added.
#' @rdname paceShrink
#' @export
setMethod("paceShrink", "PACEFit", function(object, ...) {
  slopes <- pace_shrink(object@fit, object@cellTypes,
                        resp_term = object@params$resp_term)
  object@neighbourSlopes <- as.data.frame(slopes)
  object
})

## ===========================================================================
## Stage 3: variance decomposition
## ===========================================================================

#' Variance decomposition (stage 3)
#'
#' Partitions per-gene expression variance into cell-type identity, spatial cell
#' state, contamination, and residual (plus a responder block for condition
#' cohorts), populating [varianceDecomposition()]. Needs the fitted object plus
#' the same `SpatialExperiment` used for [paceModel()] (to read the counts).
#'
#' @param object A [PACEFit] from [paceModel()].
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param ... Unused.
#' @return The `PACEFit` with the variance decomposition added.
#' @rdname paceDecompose
#' @export
setMethod("paceDecompose", "PACEFit", function(object, spe, ...) {
  genes <- object@context$genes
  Y  <- t(as.matrix(SummarizedExperiment::assay(spe, object@params$assay_name)))
  Y  <- Y[, genes, drop = FALSE]
  df <- object@context$df
  if (nrow(Y) != nrow(df))
    stop("`spe` has ", nrow(Y), " cells but the fit has ", nrow(df),
         "; pass the same object used for paceModel().", call. = FALSE)
  dec <- pace_decompose(object@fit, df, Y, object@cellTypes,
                        object@context$X_fixed, resp_term = object@params$resp_term)
  sf  <- .pace_single_frame(Y, df, object@params$celltype_col, dec)
  object@varianceDecomposition <- list(perGene = sf, blocks = dec)
  object
})

## ===========================================================================
## Stage 4: per-pair driver tables
## ===========================================================================

#' Per-pair driver scores (stage 4)
#'
#' Ranks the genes mediating each focal-neighbour relationship by driver score,
#' populating [topDrivers()]. Requires [paceShrink()] and [paceDecompose()] to
#' have run first.
#'
#' @param object A [PACEFit] with shrunken slopes and a decomposition.
#' @param pairs Optional list of focal-neighbour pairs to score; `NULL` scores
#'   all pairs.
#' @param ... Unused.
#' @return The `PACEFit` with the driver tables added.
#' @rdname paceDrivers
#' @export
setMethod("paceDrivers", "PACEFit", function(object, pairs = NULL, ...) {
  if (nrow(object@neighbourSlopes) == 0L || length(object@varianceDecomposition) == 0L)
    stop("Run paceShrink() and paceDecompose() before paceDrivers().", call. = FALSE)
  drivers <- tryCatch(
    pace_top_drivers(object@fit, object@neighbourSlopes,
                     object@varianceDecomposition$blocks, object@cellTypes,
                     pairs = pairs),
    error = function(e) {
      warning("topDrivers could not be computed: ", conditionMessage(e), call. = FALSE)
      list()
    })
  object@topDrivers <- as.list(drivers)
  object
})

## ===========================================================================
## One-shot convenience: run the whole pipeline
## ===========================================================================

#' Fit a PACE model to a SpatialExperiment
#'
#' `paceFit()` runs the full pipeline in one call:
#' [paceModel()] -> [paceShrink()] -> [paceDecompose()] -> [paceDrivers()]. For
#' step-by-step control (inspecting the fitted model, re-running the downstream
#' without refitting), call those stages directly.
#'
#' @param object A [SpatialExperiment::SpatialExperiment].
#' @param pairs Optional list of focal-neighbour pairs for the driver tables.
#' @param ... Arguments passed to [paceModel()] (`celltype_col`, `image_col`,
#'   `contamination`, `dispersion`, ...).
#' @return A fully populated [PACEFit].
#' @examples
#' \dontrun{
#' fit <- paceFit(spe, celltype_col = "cellType")
#' neighbourSlopes(fit)
#' }
#' @rdname paceFit
#' @export
setMethod("paceFit", "SpatialExperiment", function(object, ..., pairs = NULL) {
  fit <- paceModel(object, ...)
  fit <- paceShrink(fit)
  fit <- paceDecompose(fit, object)
  paceDrivers(fit, pairs = pairs)
})

## ===========================================================================
## Fit-construction primitives (inspect the neighbourhood / contamination field)
## ===========================================================================

#' Build the neighbourhood kernels
#'
#' Returns the per-cell kernel-weighted neighbour abundances: the Gaussian
#' biological kernel `K_bio` (used for the proximity coefficients) and the
#' short-range exponential technical kernel `K_tech`. Useful for inspecting the
#' neighbourhood before or independently of a fit.
#'
#' @param object A [SpatialExperiment::SpatialExperiment].
#' @param celltype_col colData column with the cell-type annotation.
#' @param h_bio,h_tech Bandwidths in micrometres (defaults 30, 5).
#' @param eps Truncation radius; defaults to `3 * h_bio`.
#' @param ... Passed to the kernel builder.
#' @return A list with `K_bio` and `K_tech` (cells x cell types).
#' @rdname buildNeighbourhood
#' @export
setMethod("buildNeighbourhood", "SpatialExperiment",
  function(object, celltype_col, h_bio = 30, h_tech = 5, eps = NULL, ...) {
    ct     <- as.character(SummarizedExperiment::colData(object)[[celltype_col]])
    types  <- sort(unique(ct))
    coords <- SpatialExperiment::spatialCoords(object)
    if (is.null(eps)) eps <- 3 * h_bio
    pace_neighbour_kernel(coords, ct, types, h_bio = h_bio, h_tech = h_tech,
                          eps = eps, ...)
  })

#' Build the ambient contamination field
#'
#' Returns the sparse cross-cell-type ambient weight matrix that defines the
#' per-cell contamination term (the short-range technical field a cell receives
#' from its heterotypic neighbours).
#'
#' @param object A [SpatialExperiment::SpatialExperiment].
#' @param celltype_col colData column with the cell-type annotation.
#' @param image_col Optional colData column grouping cells into images.
#' @param h_tech Technical bandwidth in micrometres (default 5).
#' @param assay_name Counts assay name (default `"counts"`).
#' @param ... Passed to the ambient-field builder.
#' @return The ambient-field object (sparse weight matrix and image index).
#' @rdname ambientField
#' @export
setMethod("ambientField", "SpatialExperiment",
  function(object, celltype_col, image_col = NULL, h_tech = 5,
           assay_name = "counts", ...) {
    ct     <- as.character(SummarizedExperiment::colData(object)[[celltype_col]])
    types  <- sort(unique(ct))
    coords <- SpatialExperiment::spatialCoords(object)
    Y      <- t(as.matrix(SummarizedExperiment::assay(object, assay_name)))
    image  <- factor(if (!is.null(image_col))
                       as.character(SummarizedExperiment::colData(object)[[image_col]])
                     else rep("all", nrow(Y)))
    pace_ambient_field(coords, Y, ct, image, types, h_tech = h_tech, ...)
  })
