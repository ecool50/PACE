# Regression tests for the 0.99.1 correctness fixes. Each one pins a defect
# that produced wrong numbers with no error raised, so the failure mode cannot
# come back silently.

make_condition_spe <- function(n = 600L, seed = 2L, condition_name = "Outcome") {
  set.seed(seed)
  types <- c("A", "B", "C")
  ct <- sample(types, n, replace = TRUE)
  coords <- cbind(x = runif(n, 0, 500), y = runif(n, 0, 500))

  markers <- do.call(rbind, lapply(types, function(t) {
    lam <- ifelse(ct == t, 8, 1)
    t(vapply(seq_len(12L), function(i) rpois(n, lam), numeric(n)))
  }))
  rownames(markers) <- paste0("g", seq_len(nrow(markers)))
  colnames(markers) <- paste0("c", seq_len(n))

  # two images per condition arm so the condition varies between images
  image <- rep(paste0("img", 1:4), length.out = n)
  cond  <- ifelse(image %in% c("img1", "img2"), "Good", "Bad")

  cd <- S4Vectors::DataFrame(cellType = ct, image = image)
  cd[[condition_name]] <- cond

  SpatialExperiment::SpatialExperiment(
    assays        = list(counts = markers),
    colData       = cd,
    spatialCoords = coords)
}

test_that("an NA in a random-effect variable errors instead of recycling covariates", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_condition_spe()
  # One unlabelled cell, which is routine in real cohorts. model.matrix drops
  # the row; without the guard sparseMatrix recycles x and every later cell
  # receives another cell's covariate.
  SummarizedExperiment::colData(spe)$Outcome[1L] <- NA

  expect_error(
    paceModel(spe, celltype_col = "cellType", condition_col = "Outcome",
            image_col = "image", contamination = "none",
            n_iter = 3L, threads = 1L, verbose = FALSE),
    "dropped|NA")
})

test_that("a condition column not named 'Responder' still gets a live responder block", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_condition_spe(condition_name = "Outcome")

  fit <- paceModel(spe, celltype_col = "cellType", condition_col = "Outcome",
                 image_col = "image", contamination = "none",
                 n_iter = 5L, threads = 1L, verbose = FALSE)

  # The term is paste0(condition_col, level), never the literal "Responder".
  expect_match(fit@params$resp_term, "^Outcome")

  # The indicator the decomposition and driver scoring read must not be all
  # zeros, which is how the hardcoded prefix silently zeroed the whole block.
  dummy <- fit@context$df$.resp_dummy
  expect_true(!is.null(dummy))
  expect_true(any(dummy == 1L))
  expect_true(any(dummy == 0L))
})

test_that("resp_term that names no level of condition_col is refused", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_condition_spe(condition_name = "Outcome")

  expect_error(
    paceModel(spe, celltype_col = "cellType", condition_col = "Outcome",
            image_col = "image", contamination = "none",
            resp_term = "ResponderPD",
            n_iter = 3L, threads = 1L, verbose = FALSE),
    "does not name a level")
})

test_that("a converged fit carries no non-finite coefficients and reports the count", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_condition_spe()

  fit <- paceModel(spe, celltype_col = "cellType", condition_col = "Outcome",
                 image_col = "image", contamination = "none",
                 n_iter = 5L, threads = 1L, verbose = FALSE)

  eng <- fit@fit
  # The convergence metric used to mask non-finite values out of both the mean
  # and the max, so a fit could lose genes and still report success. Every
  # engine now records what the mask hid.
  expect_true(!is.null(eng$history$n_nonfinite))
  expect_equal(sum(eng$history$n_nonfinite, na.rm = TRUE), 0)
  if (isTRUE(eng$converged)) {
    expect_true(all(is.finite(eng$B)))
    expect_true(all(is.finite(eng$U)))
  }
  # Early exit at convergence used to leave these entirely NA.
  expect_true(all(is.finite(eng$se_B)))
  expect_true(all(is.finite(eng$se_U)))
})
