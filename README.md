# PACE

**Proximity-Associated Changes in Expression** — a hierarchical empirical Bayes
framework for imaging-based spatial transcriptomics (Xenium, CosMx).

PACE quantifies how a cell's gene expression changes with proximity to specific
neighbouring cell types. It:

- fits hierarchical **negative binomial mixed models** with partial pooling
  across cell types;
- corrects for **transcript contamination** between adjacent cells with a
  per-cell ambient term;
- **decomposes expression variance** into cell-type identity, spatial cell state,
  contamination, and residual components;
- ranks the **genes that drive** each focal–neighbour spatial relationship.

Input and output follow [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment)
conventions.

## Installation

```r
# install.packages("BiocManager")
BiocManager::install("remotes")
remotes::install_github("ecool50/PACE")
```

## Quick start

```r
library(PACE)
library(SpatialExperiment)

fit <- paceFit(spe, celltype_col = "cellType")

neighbourSlopes(fit)          # shrunken (gene, focal, neighbour) slopes + lfsr
varianceDecomposition(fit)    # per-gene variance blocks
topDrivers(fit)               # per-pair driver-score tables
```

For a multi-sample, condition-stratified cohort:

```r
fit <- paceFit(spe,
               celltype_col     = "celltype",
               image_col        = "sampleID",
               condition_col     = "response",
               kernel_per_image = TRUE,
               image_re         = "condition_slopes")
```

## Status

Early release (`0.99.x`), targeting Bioconductor. The statistical engine is the
locked implementation used in the PACE manuscript.

## License

GPL-3.
