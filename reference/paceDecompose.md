# Variance decomposition

Partitions per-gene expression variance into cell-type identity, spatial
cell state, contamination, and residual (plus a responder block for
condition cohorts), populating
[`varianceDecomposition()`](https://ecool50.github.io/PACE/reference/varianceDecomposition.md).
Needs the fitted object plus the same `SpatialExperiment` used for
[`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md)
(to read the counts).

## Usage

``` r
paceDecompose(object, ...)

# S4 method for class 'PACEFit'
paceDecompose(object, spe, ...)
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  from
  [`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md).

- ...:

  Unused.

- spe:

  The
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  that was fitted.

## Value

The `PACEFit` with the variance decomposition added.
