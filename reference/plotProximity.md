# Expression versus number of neighbours (per-bin boxplot)

Reproduces the manuscript proximity figure (verbatim style port of
`density_bin_boxplot` in the analysis helpers): for each focal cell the
number of `neighbour`-type cells within `radius` micrometres is counted
and binned, and each gene's raw counts are drawn as a boxplot per bin
with the per-bin mean overlaid as a point. Because these genes are
zero-inflated the box median sits at zero in most bins, so the mean
point carries the trend.

## Usage

``` r
plotProximity(
  object,
  spe,
  genes,
  focal,
  neighbour,
  radius = 30,
  breaks = c(0, 2, 4, 6, 8, Inf),
  box_colour = "#4F8B5E",
  assay_name = "counts"
)
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  (supplies the cell-type column used at fit time).

- spe:

  The
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  that was fitted.

- genes:

  Character vector of genes (shown as facets).

- focal, neighbour:

  Focal and neighbour cell types.

- radius:

  Neighbour search radius in micrometres (default 30).

- breaks:

  Neighbour-count bin breaks (default `c(0, 2, 4, 6, 8, Inf)`).

- box_colour:

  Box and mean-point colour.

- assay_name:

  Counts assay name (default `"counts"`).

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE)  plotProximity(fit, spe, c("MRC1", "APOC1"), "Macrophage", "Tumour")  # \dontrun{}
```
