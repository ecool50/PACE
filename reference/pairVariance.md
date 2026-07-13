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
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
head(pairVariance(fit))
#> # A tibble: 6 × 3
#>   focal  neighbour          val
#>   <chr>  <chr>            <dbl>
#> 1 B_Cell Dendritic_Cell 0      
#> 2 B_Cell Endothelial    0      
#> 3 B_Cell Macrophage     0.00619
#> 4 B_Cell Myoepithelial  0      
#> 5 B_Cell Stromal        0      
#> 6 B_Cell T_Cell         2.72   
```
