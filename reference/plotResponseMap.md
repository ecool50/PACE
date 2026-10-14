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
plotResponseMap(fit, spe, "SPP1", "Macrophage", "Tumour")

# }
```
