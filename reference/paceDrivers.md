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
paceDrivers(object, spe = NULL, pairs = NULL, ...)
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  with shrunken slopes and a decomposition.

- ...:

  Unused.

- spe:

  The
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  that was fitted. Needed only when the fit does not store its fitted
  means.

- pairs:

  Optional list of focal-neighbour pairs to score; `NULL` scores all
  pairs.

## Value

The `PACEFit` with the driver tables added.

## Details

The driver scores read the fitted means. A fit saved without its `n x G`
matrices is rebuilt exactly from the fit and `spe`, as in
[`paceDecompose()`](https://ecool50.github.io/PACE/reference/paceDecompose.md),
so pass `spe` for such a fit.

## Examples

``` r
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
fit <- paceDrivers(fit, spe)
names(topDrivers(fit))
#>  [1] "B_Cell_Dendritic_Cell"        "B_Cell_Endothelial"          
#>  [3] "B_Cell_Macrophage"            "B_Cell_Myoepithelial"        
#>  [5] "B_Cell_Stromal"               "B_Cell_T_Cell"               
#>  [7] "B_Cell_Tumour"                "Dendritic_Cell_B_Cell"       
#>  [9] "Dendritic_Cell_Endothelial"   "Dendritic_Cell_Macrophage"   
#> [11] "Dendritic_Cell_Myoepithelial" "Dendritic_Cell_Stromal"      
#> [13] "Dendritic_Cell_T_Cell"        "Dendritic_Cell_Tumour"       
#> [15] "Endothelial_B_Cell"           "Endothelial_Dendritic_Cell"  
#> [17] "Endothelial_Macrophage"       "Endothelial_Myoepithelial"   
#> [19] "Endothelial_Stromal"          "Endothelial_T_Cell"          
#> [21] "Endothelial_Tumour"           "Macrophage_B_Cell"           
#> [23] "Macrophage_Dendritic_Cell"    "Macrophage_Endothelial"      
#> [25] "Macrophage_Myoepithelial"     "Macrophage_Stromal"          
#> [27] "Macrophage_T_Cell"            "Macrophage_Tumour"           
#> [29] "Myoepithelial_B_Cell"         "Myoepithelial_Dendritic_Cell"
#> [31] "Myoepithelial_Endothelial"    "Myoepithelial_Macrophage"    
#> [33] "Myoepithelial_Stromal"        "Myoepithelial_T_Cell"        
#> [35] "Myoepithelial_Tumour"         "Stromal_B_Cell"              
#> [37] "Stromal_Dendritic_Cell"       "Stromal_Endothelial"         
#> [39] "Stromal_Macrophage"           "Stromal_Myoepithelial"       
#> [41] "Stromal_T_Cell"               "Stromal_Tumour"              
#> [43] "T_Cell_B_Cell"                "T_Cell_Dendritic_Cell"       
#> [45] "T_Cell_Endothelial"           "T_Cell_Macrophage"           
#> [47] "T_Cell_Myoepithelial"         "T_Cell_Stromal"              
#> [49] "T_Cell_Tumour"                "Tumour_B_Cell"               
#> [51] "Tumour_Dendritic_Cell"        "Tumour_Endothelial"          
#> [53] "Tumour_Macrophage"            "Tumour_Myoepithelial"        
#> [55] "Tumour_Stromal"               "Tumour_T_Cell"               
```
