# Pairwise spatial variance table

Returns the tidy per-(focal, neighbour) spatial variance table that
[`plotPairHeatmap()`](https://ecool50.github.io/PACE/reference/plotPairHeatmap.md)
draws: for each focal cell type, the percentage of its total expression
variance contributed by spatial interaction with each neighbour,
obtained by splitting the focal's spatial share across neighbours by a
normalised (off-diagonal-renormalised) Pratt attribution. Off-diagonal
pairs only.

## Usage

``` r
pairVariance(object, ...)

# S4 method for class 'PACEFit'
pairVariance(object, block = c("spatial", "responder"), ...)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- ...:

  Unused.

- block:

  Either `"spatial"` (default) or `"responder"` (condition cohorts).

## Value

A data frame with columns `focal`, `neighbour`, and `val` (spatial % of
the focal type's total variance).

## Examples

``` r
if (FALSE)  pairVariance(fit)  # \dontrun{}
```
