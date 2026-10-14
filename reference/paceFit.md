# Fit a PACE model to a SpatialExperiment

`paceFit()` runs the full pipeline in one call:
[`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md)
-\>
[`paceShrink()`](https://ecool50.github.io/PACE/reference/paceShrink.md)
-\>
[`paceDecompose()`](https://ecool50.github.io/PACE/reference/paceDecompose.md)
-\>
[`paceDrivers()`](https://ecool50.github.io/PACE/reference/paceDrivers.md).
For step-by-step control (inspecting the fitted model, re-running the
downstream without refitting), call those stages directly.

## Usage

``` r
paceFit(object, ...)

# S4 method for class 'SpatialExperiment'
paceFit(object, ..., pairs = NULL)
```

## Arguments

- object:

  A
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- ...:

  Arguments passed to
  [`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md)
  (`celltype_col`, `image_col`, `contamination`, `dispersion`, ...).

- pairs:

  Optional list of focal-neighbour pairs for the driver tables.

## Value

A fully populated
[PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

## Examples

``` r
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
# \donttest{
fit <- paceFit(spe, celltype_col = "cellType", verbose = FALSE)
head(neighbourSlopes(fit))
#>     gene  focal neighbour   term      estimate   std.error estimate_shrunk
#> 1 ABCC11 B_Cell    B_Cell B_Cell -1.968907e-05 0.001376136               0
#> 2  ACTA2 B_Cell    B_Cell B_Cell -9.976956e-05 0.004245639               0
#> 3  ACTG2 B_Cell    B_Cell B_Cell  1.418293e-03 0.005170493               0
#> 4  ADAM9 B_Cell    B_Cell B_Cell  3.277669e-05 0.003072945               0
#> 5 ADGRE5 B_Cell    B_Cell B_Cell  8.125365e-04 0.005304314               0
#> 6  ADH1B B_Cell    B_Cell B_Cell -3.142850e-05 0.002131910               0
#>   sd_shrunk lfsr
#> 1         0    1
#> 2         0    1
#> 3         0    1
#> 4         0    1
#> 5         0    1
#> 6         0    1
# }
```
