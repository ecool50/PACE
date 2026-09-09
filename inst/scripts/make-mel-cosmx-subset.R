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
## The subset keeps EVERY cell and EVERY patient and reduces only the panel, to
## the 180 most widely detected genes by the same max-per-cell-type criterion
## already used for QC. File size is cells x genes, and of the two the cells are
## what must be preserved: the model reads local neighbourhoods, and the
## condition effect is a between-patient contrast, so cropping cells or dropping
## patients destroys exactly what this vignette demonstrates while cutting the
## panel does not. Cropping cells was tried first and abandoned; at 15% of each
## section nothing reached lfsr < 0.05 at all, and no cell-based crop under
## about 5 MB recovered the macrophage SPP1 result, including crops targeted at
## the macrophages themselves, which keep most of the tissue because the
## macrophages are spread through all of it.
##
## SPP1 sits at rank 173 by detection, so it is forced into the panel rather
## than arriving on merit; every other gene is taken in detection order.
##
## The panel cannot go much lower. At 120 genes some focal cell type has nothing
## left to decompose and mvpql_variance_decomposition_multi() fails with
## "Column `focal` is not found".
##
## The reduced panel leaves estimates alone but shifts mash's calibration, since
## it calibrates against the genes it is given: this subset yields 76 responder
## calls where the full panel yields 46. SPP1 Macrophage <- Tumour is -0.0663
## here against -0.0668 on the full panel.

library(SpatialExperiment)
library(zellkonverter)

## ... published preparation of Melanoma_5612.h5ad, yielding `Y` (cells x genes
## integer counts) and `df` (celltype, imageID, Responder, x, y) ...

## Panel reduction: all cells, all patients, 180 most-detected genes.
ct   <- as.character(df$celltype)
cbt  <- lapply(cell_types, function(t) which(ct == t))
maxdet <- apply(Y, 2, function(g) max(vapply(cbt, function(i) mean(g[i] > 0), numeric(1))))
keep <- names(sort(maxdet, decreasing = TRUE))[seq_len(180L)]
if (!"SPP1" %in% keep) keep <- c(keep[-180L], "SPP1")
keep <- sort(keep)

spe <- SpatialExperiment(
  assays  = list(counts = as(t(Y[, keep, drop = FALSE]), "dgCMatrix")),
  colData = DataFrame(
    cellType  = ct,
    Responder = factor(as.character(df$Responder), levels = c("nonPD", "PD")),
    image     = as.character(df$imageID)),
  spatialCoords = as.matrix(df[, c("x", "y")]))

saveRDS(spe, "inst/extdata/mel_cosmx_subset.rds", compress = "xz")
