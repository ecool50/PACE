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
if (FALSE)  anchorGenes(fit, spe)  # \dontrun{}
```
