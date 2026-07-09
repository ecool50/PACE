# Shrink the neighbour slopes

Stabilises the fitted per-(gene, focal, neighbour) proximity slopes with
multivariate adaptive shrinkage, populating
[`neighbourSlopes()`](https://ecool50.github.io/PACE/reference/neighbourSlopes.md).

## Usage

``` r
paceShrink(object, ...)

# S4 method for class 'PACEFit'
paceShrink(object, ...)
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  from
  [`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md).

- ...:

  Unused.

## Value

The `PACEFit` with the shrunken neighbour slopes added.
