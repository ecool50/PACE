# Build the neighbourhood kernels

Returns the per-cell kernel-weighted neighbour abundances: the Gaussian
biological kernel `K_bio` (used for the proximity coefficients) and the
short-range exponential technical kernel `K_tech`. Useful for inspecting
the neighbourhood before or independently of a fit.

## Usage

``` r
buildNeighbourhood(object, ...)

# S4 method for class 'SpatialExperiment'
buildNeighbourhood(
  object,
  celltype_col,
  h_bio = 30,
  h_tech = 5,
  eps = NULL,
  ...
)
```

## Arguments

- object:

  A
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- ...:

  Passed to the kernel builder.

- celltype_col:

  colData column with the cell-type annotation.

- h_bio, h_tech:

  Bandwidths in micrometres (defaults 30, 5).

- eps:

  Truncation radius; defaults to `3 * h_bio`.

## Value

A list with `K_bio` and `K_tech` (cells x cell types).

## Examples

``` r
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
nb <- buildNeighbourhood(spe, celltype_col = "cellType")
dim(nb$K_bio)
#> [1] 7898    8
```
