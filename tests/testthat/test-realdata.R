# Regression test on the shipped breast cancer Xenium subset: PACE must recover
# the macrophage-tumour signal (MRC1 down, APOC1 up near tumour). Guards against
# engine/packaging regressions in the contamination-correcting streaming path.
test_that("paceFit recovers the macrophage-tumour signal on the BC subset", {
  skip_if_not_installed("SpatialExperiment")
  f <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(f == "", "example dataset not installed")

  spe <- readRDS(f)
  fit <- paceFit(spe, celltype_col = "cellType",
                 contamination = "percell_hc", dispersion = "nb1",
                 n_iter = 32L, threads = 2L, verbose = FALSE)

  ns <- neighbourSlopes(fit)
  mt <- ns[ns$focal == "Macrophage" & ns$neighbour == "Tumour", ]
  get <- function(g) mt$estimate_shrunk[mt$gene == g]

  expect_lt(get("MRC1"),  0)      # tissue-resident marker reduced near tumour
  expect_gt(get("APOC1"), 0)      # lipid-associated marker elevated near tumour
  expect_lt(mt$lfsr[mt$gene == "MRC1"],  0.05)
  expect_lt(mt$lfsr[mt$gene == "APOC1"], 0.05)
})
