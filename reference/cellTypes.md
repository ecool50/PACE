# Cell types modelled

The cell types the fit was built on, in fitting order. The order
matters: it indexes the rows of the random-effect design, so the
readouts and the solver state agree on what a neighbour column refers
to.

## Usage

``` r
cellTypes(object, ...)

# S4 method for class 'PACEFit'
cellTypes(object, ...)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- ...:

  Unused.

## Value

A character vector of cell type names.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
cellTypes(fit)
#> [1] "B_Cell"         "Dendritic_Cell" "Endothelial"    "Macrophage"    
#> [5] "Myoepithelial"  "Stromal"        "T_Cell"         "Tumour"        
```
