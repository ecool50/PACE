# Contamination anchor genes per cell type

Returns, for each cell type, the genes used as negative-control anchors
when identifying the per-cell contamination loading: genes owned by
another cell type (high mean expression there) that are near-absent in
the focal type. Anchors are recomputed with the same selection as the
fit, and ranked by expression in their owning cell type. This is the
source of the manuscript's supplementary anchor-gene tables.

## Usage

``` r
anchorGenes(object, spe, top_n = 10L, assay_name = "counts")
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  (supplies the cell-type and image columns used at fit time).

- spe:

  The
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  that was fitted.

- top_n:

  Number of top anchors to list per cell type (default 10).

- assay_name:

  Counts assay name (default `"counts"`).

## Value

A data frame with one row per cell type: `cellType`, `n_anchors` (total
anchors), and `top_anchors` (the `top_n` strongest, comma-separated).

## Examples

``` r
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
anchorGenes(fit, spe)
#>     anchors: homotypic-core cells 5168/7898; anchors/type median=136 range=[129,149]
#> Warning: NAs introduced by coercion
#>         cellType n_anchors
#> 1         B_Cell       129
#> 2 Dendritic_Cell       141
#> 3    Endothelial       129
#> 4     Macrophage       136
#> 5  Myoepithelial       136
#> 6        Stromal       131
#> 7         T_Cell       149
#> 8         Tumour       141
#>                                                            top_anchors
#> 1  ERBB2, KRT7, POSTN, ANKRD30A, FOXA1, EPCAM, SCD, GATA3, CCND1, MLPH
#> 2   ERBB2, KRT7, ANKRD30A, FOXA1, EPCAM, SCD, GATA3, CCND1, MLPH, FASN
#> 3 ERBB2, KRT7, ANKRD30A, FOXA1, EPCAM, SCD, GATA3, MLPH, FASN, TACSTD2
#> 4   ERBB2, KRT7, ANKRD30A, FOXA1, EPCAM, SCD, GATA3, CCND1, MLPH, FASN
#> 5       POSTN, ANKRD30A, FOXA1, LUM, SCD, MLPH, FASN, CXCR4, TCIM, LYZ
#> 6 KRT7, ANKRD30A, FOXA1, EPCAM, SCD, GATA3, CCND1, MLPH, FASN, TACSTD2
#> 7  ERBB2, KRT7, POSTN, ANKRD30A, FOXA1, EPCAM, SCD, GATA3, CCND1, MLPH
#> 8     POSTN, LUM, MYLK, SERPINA3, KRT14, KRT5, LYZ, KRT6B, VWF, PECAM1
```
