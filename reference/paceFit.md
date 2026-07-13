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
#>  - Computing 278 x 313 likelihood matrix.
#>  - Likelihood calculations took 0.08 seconds.
#>  - Fitting model with 313 mixture components.
#>  - Model fitting took 0.11 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] B_Cell: 278 genes shrunk; sig (lfsr<0.05) = 0
#>  - Computing 278 x 92 likelihood matrix.
#>  - Likelihood calculations took 0.02 seconds.
#>  - Fitting model with 92 mixture components.
#>  - Model fitting took 0.03 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Dendritic_Cell: 278 genes shrunk; sig (lfsr<0.05) = 0
#>  - Computing 278 x 404 likelihood matrix.
#>  - Likelihood calculations took 0.10 seconds.
#>  - Fitting model with 404 mixture components.
#>  - Model fitting took 0.15 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Endothelial: 278 genes shrunk; sig (lfsr<0.05) = 2
#>  - Computing 278 x 417 likelihood matrix.
#>  - Likelihood calculations took 0.10 seconds.
#>  - Fitting model with 417 mixture components.
#>  - Model fitting took 0.22 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Macrophage: 278 genes shrunk; sig (lfsr<0.05) = 5
#>  - Computing 278 x 92 likelihood matrix.
#>  - Likelihood calculations took 0.01 seconds.
#>  - Fitting model with 92 mixture components.
#>  - Model fitting took 0.03 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Myoepithelial: 278 genes shrunk; sig (lfsr<0.05) = 0
#>  - Computing 278 x 391 likelihood matrix.
#>  - Likelihood calculations took 0.09 seconds.
#>  - Fitting model with 391 mixture components.
#>  - Model fitting took 0.12 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Stromal: 278 genes shrunk; sig (lfsr<0.05) = 3
#>  - Computing 278 x 430 likelihood matrix.
#>  - Likelihood calculations took 0.10 seconds.
#>  - Fitting model with 430 mixture components.
#>  - Model fitting took 0.14 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 278 x 628 likelihood matrix.
#>  - Likelihood calculations took 0.16 seconds.
#>  - Fitting model with 628 mixture components.
#>  - Model fitting took 0.19 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] T_Cell: 278 genes shrunk; sig (lfsr<0.05) = 21
#>  - Computing 278 x 404 likelihood matrix.
#>  - Likelihood calculations took 0.10 seconds.
#>  - Fitting model with 404 mixture components.
#>  - Model fitting took 0.34 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 278 x 590 likelihood matrix.
#>  - Likelihood calculations took 0.15 seconds.
#>  - Fitting model with 590 mixture components.
#>  - Model fitting took 0.68 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>   [mashr] Tumour: 278 genes shrunk; sig (lfsr<0.05) = 29
head(neighbourSlopes(fit))
#>     gene  focal neighbour   term      estimate   std.error estimate_shrunk
#> 1 ABCC11 B_Cell    B_Cell B_Cell -1.968913e-05 0.001376136               0
#> 2  ACTA2 B_Cell    B_Cell B_Cell -9.977158e-05 0.004245640               0
#> 3  ACTG2 B_Cell    B_Cell B_Cell  1.418290e-03 0.005170495               0
#> 4  ADAM9 B_Cell    B_Cell B_Cell  3.277721e-05 0.003072945               0
#> 5 ADGRE5 B_Cell    B_Cell B_Cell  8.125332e-04 0.005304315               0
#> 6  ADH1B B_Cell    B_Cell B_Cell -3.142839e-05 0.002131910               0
#>   sd_shrunk lfsr
#> 1         0    1
#> 2         0    1
#> 3         0    1
#> 4         0    1
#> 5         0    1
#> 6         0    1
# }
```
