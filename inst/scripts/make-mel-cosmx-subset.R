## How inst/extdata/mel_cosmx_subset.rds was made.
##
## SOURCE. The CosMx melanoma cohort of Dong, Su, Kluger, Fan and Kluger,
## "SIMVI disentangles intrinsic and spatial-induced cellular states in spatial
## omics data", deposited at Zenodo, doi:10.5281/zenodo.14708000, file
## Melanoma_5612.h5ad. The deposit is licensed CC-BY-4.0, which permits
## redistribution of this derived subset with attribution; the subset is a
## modification of that data and is noted as such in the vignette.
##
## DERIVATION. The published preparation is applied first, exactly as in the
## accompanying paper: keep cells with a recorded best response, recode to PD
## versus non-PD (non-PD pools SD, PR and CR), convert pixel coordinates to
## micrometres by 0.12028, pool the six Tumor_* sub-clusters into Tumour and the
## T-cell subtypes into T_Cell, keep the six modelled types, drop genes detected
## in fewer than 5% of the cells of every type, and drop zero-count cells. That
## gives 56,274 cells and 927 genes across 26 sections, one per patient.
##
## The subset then keeps ALL 26 patients and crops each section spatially: the
## cells nearest that section's centroid, 15% of it or 120 cells, whichever is
## larger. Cropping contiguously rather than sampling at random is deliberate,
## because the model reads local neighbourhoods and a random subsample would
## thin them. Keeping every patient is also deliberate: the condition effect is
## a between-patient contrast, so dropping patients would weaken the very thing
## the vignette illustrates. The result is 8,454 cells, 927 genes, 26 sections.
##
## The subset is sized for a package vignette, not for inference. It does NOT
## reproduce the cohort-level macrophage SPP1 result; see the vignette, which
## says so and explains why.

library(SpatialExperiment)
library(zellkonverter)

## ... published preparation of Melanoma_5612.h5ad, yielding `Y` (cells x genes
## integer counts) and `df` (celltype, imageID, Responder, x, y) ...

crop <- function(df, frac = 0.15, min_cells = 120L) {
  sort(unlist(lapply(split(seq_len(nrow(df)), df$imageID), function(idx) {
    xy <- as.matrix(df[idx, c("x", "y")])
    d  <- sqrt((xy[, 1] - mean(xy[, 1]))^2 + (xy[, 2] - mean(xy[, 2]))^2)
    idx[order(d)][seq_len(max(min_cells, ceiling(frac * length(idx))))]
  }), use.names = FALSE))
}

idx <- crop(df)
spe <- SpatialExperiment(
  assays  = list(counts = as(t(Y[idx, , drop = FALSE]), "dgCMatrix")),
  colData = DataFrame(
    cellType  = as.character(df$celltype)[idx],
    Responder = factor(as.character(df$Responder)[idx], levels = c("nonPD", "PD")),
    image     = as.character(df$imageID)[idx]),
  spatialCoords = as.matrix(df[idx, c("x", "y")]))

saveRDS(spe, "inst/extdata/mel_cosmx_subset.rds", compress = "xz")
