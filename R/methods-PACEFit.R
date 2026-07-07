#' Shrunken neighbour proximity slopes
#'
#' Returns the shrunken per-(gene, focal, neighbour) proximity coefficients with
#' their local false sign rate (lfsr). For condition cohorts the responder
#' interaction terms are included.
#'
#' @param object A [PACEFit].
#' @param ... Unused.
#' @return A data frame of shrunken slopes.
#' @examples
#' \dontrun{ neighbourSlopes(fit) }
#' @rdname neighbourSlopes
#' @export
setMethod("neighbourSlopes", "PACEFit", function(object, ...) {
  object@neighbourSlopes
})

#' Per-gene variance decomposition
#'
#' Returns the per-gene, per-focal variance decomposition. By default the
#' observed single-frame decomposition (cell-type identity, spatial cell state,
#' contamination, residual; and a responder spatial block for condition cohorts)
#' is returned; set `which = "blocks"` for the underlying link-scale block tables.
#'
#' @param object A [PACEFit].
#' @param which Either `"perGene"` (default, the observed single-frame table) or
#'   `"blocks"` (the raw link-scale decomposition).
#' @param ... Unused.
#' @return A data frame (for `"perGene"`) or a list (for `"blocks"`).
#' @examples
#' \dontrun{ varianceDecomposition(fit) }
#' @rdname varianceDecomposition
#' @export
setMethod("varianceDecomposition", "PACEFit",
  function(object, which = c("perGene", "blocks"), ...) {
    which <- match.arg(which)
    if (which == "perGene") object@varianceDecomposition$perGene
    else object@varianceDecomposition$blocks
  })

#' Per-pair driver tables
#'
#' Returns the driver-score tables ranking the genes that mediate each
#' focal-neighbour spatial relationship.
#'
#' @param object A [PACEFit].
#' @param ... Unused.
#' @return A list of per-pair driver tables.
#' @examples
#' \dontrun{ topDrivers(fit) }
#' @rdname topDrivers
#' @export
setMethod("topDrivers", "PACEFit", function(object, ...) {
  object@topDrivers
})

#' @importFrom methods show
#' @rdname PACEFit-class
#' @export
setMethod("show", "PACEFit", function(object) {
  p <- object@params
  cat("class: PACEFit\n")
  cat("cell types (", length(object@cellTypes), "): ",
      paste(object@cellTypes, collapse = ", "), "\n", sep = "")
  cat("neighbour slopes: ", nrow(object@neighbourSlopes),
      " (gene, focal, neighbour) rows\n", sep = "")
  n_sig <- sum(object@neighbourSlopes$lfsr < 0.05, na.rm = TRUE)
  cat("  significant (lfsr < 0.05): ", n_sig, "\n", sep = "")
  cat("kernels: h_bio = ", p$h_bio, " um, h_tech = ", p$h_tech, " um\n", sep = "")
  cat("contamination: ", p$contamination, "; dispersion: ", p$dispersion, "\n", sep = "")
  if (!is.null(p$condition_col)) {
    cat("condition: ", p$condition_col, " (", p$resp_term, ")\n", sep = "")
  }
  invisible(NULL)
})
