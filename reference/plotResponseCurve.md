# Response curve: binned means with the PACE slope overlaid

Population view of one gene in `focal` cells against `neighbour`
density, split by the fit's condition. The solid line and ribbon are
binned means and standard errors of the observed data; the dashed line
is the PACE slope read straight from the fitted model, anchored at each
arm's data centroid. It is the model's estimate drawn over the data, not
a smoother refitted to it.

## Usage

``` r
plotResponseCurve(
  object,
  spe,
  gene,
  focal,
  neighbour,
  n_bins = 10,
  min_bin_n = 10,
  colours = c("#3B6FB6", "#C0392B")
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

- n_bins:

  Number of density quantile bins (default 10).

- min_bin_n:

  Bins with fewer focal cells than this are dropped.

- colours:

  Colours for the two arms.

## Value

A `ggplot` object.

## Examples

``` r
spe <- readRDS(system.file("extdata", "mel_cosmx_subset.rds", package = "PACE"))
# \donttest{
fit <- paceFit(spe, celltype_col = "cellType", condition_col = "Responder",
               image_col = "image", kernel_per_image = TRUE,
               image_re = "intercept", verbose = FALSE)
plotResponseCurve(fit, spe, "SPP1", "Macrophage", "Tumour")

# }
```
