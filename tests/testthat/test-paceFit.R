# A small synthetic SpatialExperiment so the smoke test runs without the
# manuscript data. Three cell types with type-specific genes, laid out in space.
make_toy_spe <- function(n = 600L, seed = 1L) {
  set.seed(seed)
  types <- c("A", "B", "C")
  ct <- sample(types, n, replace = TRUE)
  coords <- cbind(x = runif(n, 0, 500), y = runif(n, 0, 500))

  # 12 genes per type: higher in the owning type than the others.
  markers <- do.call(rbind, lapply(types, function(t) {
    lam <- ifelse(ct == t, 8, 1)
    t(vapply(seq_len(12L), function(i) rpois(n, lam), numeric(n)))
  }))
  rownames(markers) <- paste0("g", seq_len(nrow(markers)))
  colnames(markers) <- paste0("c", seq_len(n))

  SpatialExperiment::SpatialExperiment(
    assays        = list(counts = markers),
    colData       = S4Vectors::DataFrame(cellType = ct),
    spatialCoords = coords)
}

test_that("paceFit returns a populated PACEFit from a SpatialExperiment", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_toy_spe()

  fit <- paceFit(spe, celltype_col = "cellType",
                 contamination = "none",     # no anchors needed on toy data
                 n_iter = 5L, threads = 1L, verbose = FALSE)

  expect_s4_class(fit, "PACEFit")
  expect_setequal(fit@cellTypes, c("A", "B", "C"))

  slopes <- neighbourSlopes(fit)
  expect_s3_class(slopes, "data.frame")
  expect_true(nrow(slopes) > 0)
  expect_true(all(c("gene", "focal", "neighbour", "lfsr") %in% colnames(slopes)))

  vd <- varianceDecomposition(fit)
  expect_true(!is.null(vd))

  expect_output(show(fit), "class: PACEFit")
})
