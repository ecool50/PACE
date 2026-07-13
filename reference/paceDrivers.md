# Per-pair driver scores

Ranks the genes mediating each focal-neighbour relationship by driver
score, populating
[`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md).
Requires
[`paceShrink()`](https://ecool50.github.io/PACE/reference/paceShrink.md)
and
[`paceDecompose()`](https://ecool50.github.io/PACE/reference/paceDecompose.md)
to have run first.

## Usage

``` r
paceDrivers(object, ...)

# S4 method for class 'PACEFit'
paceDrivers(object, pairs = NULL, ...)
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  with shrunken slopes and a decomposition.

- ...:

  Unused.

- pairs:

  Optional list of focal-neighbour pairs to score; `NULL` scores all
  pairs.

## Value

The `PACEFit` with the driver tables added.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
fit <- paceDrivers(fit)
#> Warning: topDrivers could not be computed: 'x' must be an array of at least two dimensions
names(topDrivers(fit))
#> NULL
```
