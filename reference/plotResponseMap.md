# Tissue maps and scatters for one exemplar image per arm

The spatial view behind
[`plotResponseCurve()`](https://ecool50.github.io/PACE/reference/plotResponseCurve.md).
For one image per condition arm: a smooth kernel-density map of the
neighbour type with the focal cells overlaid and coloured by expression
(top), and expression against neighbour density for the same cells with
a fitted line (bottom). Colour scales are shared across the panels so
the arms are directly comparable.

## Usage

``` r
plotResponseMap(
  object,
  spe,
  gene,
  focal,
  neighbour,
  images = NULL,
  kde_h = 60,
  ngrid = 200,
  scalebar_len = 100
)
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  fitted with a `condition_col`.

- spe:

  The
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  that was fitted.

- gene, focal, neighbour:

  Gene, focal cell type and neighbour cell type.

- images:

  Optional named character vector of image identifiers, one per arm. The
  default picks, for each arm, the image with the most focal cells.

- kde_h, ngrid:

  Kernel bandwidth and grid size for the density map.

- scalebar_len:

  Scale-bar length in micrometres.

## Value

A `patchwork` object.

## Examples

``` r
spe <- readRDS(system.file("extdata", "mel_cosmx_subset.rds", package = "PACE"))
# \donttest{
fit <- paceFit(spe, celltype_col = "cellType", condition_col = "Responder",
               image_col = "image", kernel_per_image = TRUE,
               image_re = "intercept", verbose = FALSE)
#>  - Computing 180 x 298 likelihood matrix.
#>  - Likelihood calculations took 0.04 seconds.
#>  - Fitting model with 298 mixture components.
#>  - Model fitting took 0.14 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] ResponderPD:B_Cell: 180 genes shrunk; sig (lfsr<0.05) = 0
#>  - Computing 180 x 375 likelihood matrix.
#>  - Likelihood calculations took 0.05 seconds.
#>  - Fitting model with 375 mixture components.
#>  - Model fitting took 0.29 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 180 x 579 likelihood matrix.
#>  - Likelihood calculations took 0.07 seconds.
#>  - Fitting model with 579 mixture components.
#>  - Model fitting took 0.54 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] ResponderPD:Endothelial: 180 genes shrunk; sig (lfsr<0.05) = 24
#>  - Computing 180 x 364 likelihood matrix.
#>  - Likelihood calculations took 0.05 seconds.
#>  - Fitting model with 364 mixture components.
#>  - Model fitting took 0.23 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] ResponderPD:Fibroblast: 180 genes shrunk; sig (lfsr<0.05) = 1
#>  - Computing 180 x 243 likelihood matrix.
#>  - Likelihood calculations took 0.03 seconds.
#>  - Fitting model with 243 mixture components.
#>  - Model fitting took 0.13 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 180 x 375 likelihood matrix.
#>  - Likelihood calculations took 0.05 seconds.
#>  - Fitting model with 375 mixture components.
#>  - Model fitting took 0.30 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] ResponderPD:Macrophage: 180 genes shrunk; sig (lfsr<0.05) = 13
#>  - Computing 180 x 364 likelihood matrix.
#>  - Likelihood calculations took 0.05 seconds.
#>  - Fitting model with 364 mixture components.
#>  - Model fitting took 0.09 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 180 x 562 likelihood matrix.
#>  - Likelihood calculations took 0.07 seconds.
#>  - Fitting model with 562 mixture components.
#>  - Model fitting took 0.13 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] ResponderPD:T_Cell: 180 genes shrunk; sig (lfsr<0.05) = 13
#>  - Computing 180 x 265 likelihood matrix.
#>  - Likelihood calculations took 0.03 seconds.
#>  - Fitting model with 265 mixture components.
#>  - Model fitting took 0.23 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 180 x 409 likelihood matrix.
#>  - Likelihood calculations took 0.05 seconds.
#>  - Fitting model with 409 mixture components.
#>  - Model fitting took 0.61 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] ResponderPD:Tumour: 180 genes shrunk; sig (lfsr<0.05) = 25
plotResponseMap(fit, spe, "SPP1", "Macrophage", "Tumour")

# }
```
