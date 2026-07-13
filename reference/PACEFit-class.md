# PACEFit: a fitted PACE model

An S4 container for the result of a PACE fit. It holds the raw penalised
quasi-likelihood solver state together with the three interpretable
output layers: shrunken per-(gene, focal, neighbour) proximity slopes,
the per-gene variance decomposition, and the per-pair driver tables.
Construct one with
[`paceFit()`](https://ecool50.github.io/PACE/reference/paceFit.md); read
it out with
[`neighbourSlopes()`](https://ecool50.github.io/PACE/reference/neighbourSlopes.md),
[`varianceDecomposition()`](https://ecool50.github.io/PACE/reference/varianceDecomposition.md),
and
[`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md).

## Usage

``` r
# S4 method for class 'PACEFit'
show(object)
```

## Slots

- `fit`:

  The raw streaming PQL solver state (fixed and random effects, per-cell
  contamination loadings, gene-wise overdispersions, variance
  components).

- `neighbourSlopes`:

  A data frame of shrunken proximity coefficients, one row per (gene,
  focal cell type, neighbour cell type, condition term), with the raw
  estimate, shrunken estimate, and local false sign rate (lfsr).

- `varianceDecomposition`:

  A list of per-gene, per-focal variance decomposition tables (cell-type
  identity, spatial cell state, contamination, residual; and, for
  condition cohorts, a responder spatial block).

- `topDrivers`:

  A list of per-pair driver-score tables ranking the genes that mediate
  each focal-neighbour spatial relationship.

- `cellTypes`:

  The cell types modelled, in fitting order.

- `params`:

  The fitting parameters actually used (bandwidths, contamination and
  dispersion settings, condition column, and so on).

- `context`:

  The working frame and fixed-effect design retained from the fit so
  that the downstream stages
  ([`paceShrink()`](https://ecool50.github.io/PACE/reference/paceShrink.md),
  [`paceDecompose()`](https://ecool50.github.io/PACE/reference/paceDecompose.md),
  [`paceDrivers()`](https://ecool50.github.io/PACE/reference/paceDrivers.md))
  can run without refitting.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
fit
#> class: PACEFit
#> cell types (8): B_Cell, Dendritic_Cell, Endothelial, Macrophage, Myoepithelial, Stromal, T_Cell, Tumour
#> kernels: h_bio = 30 um, h_tech = 5 um | contamination: percell_hc; dispersion: nb1
#> pipeline: model -> shrink -> decompose -> drivers
#>   neighbour slopes: 17792 rows (63 at lfsr < 0.05)
cellTypes <- fit@cellTypes
```
