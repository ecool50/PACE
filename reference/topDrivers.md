# Per-pair driver tables

Returns the driver-score tables ranking the genes that mediate each
focal-neighbour spatial relationship.

## Usage

``` r
topDrivers(object, ...)

# S4 method for class 'PACEFit'
topDrivers(object, ...)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- ...:

  Unused.

## Value

A list of per-pair driver tables.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
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
