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

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
fit <- paceShrink(fit)
#>  - Computing 278 x 313 likelihood matrix.
#>  - Likelihood calculations took 0.07 seconds.
#>  - Fitting model with 313 mixture components.
#>  - Model fitting took 0.09 seconds.
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
#>  - Likelihood calculations took 0.09 seconds.
#>  - Fitting model with 404 mixture components.
#>  - Model fitting took 0.14 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Endothelial: 278 genes shrunk; sig (lfsr<0.05) = 2
#>  - Computing 278 x 404 likelihood matrix.
#>  - Likelihood calculations took 0.09 seconds.
#>  - Fitting model with 404 mixture components.
#>  - Model fitting took 0.21 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Macrophage: 278 genes shrunk; sig (lfsr<0.05) = 5
#>  - Computing 278 x 92 likelihood matrix.
#>  - Likelihood calculations took 0.00 seconds.
#>  - Fitting model with 92 mixture components.
#>  - Model fitting took 0.03 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Myoepithelial: 278 genes shrunk; sig (lfsr<0.05) = 0
#>  - Computing 278 x 391 likelihood matrix.
#>  - Likelihood calculations took 0.09 seconds.
#>  - Fitting model with 391 mixture components.
#>  - Model fitting took 0.11 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Stromal: 278 genes shrunk; sig (lfsr<0.05) = 3
#>  - Computing 278 x 430 likelihood matrix.
#>  - Likelihood calculations took 0.10 seconds.
#>  - Fitting model with 430 mixture components.
#>  - Model fitting took 0.13 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 278 x 628 likelihood matrix.
#>  - Likelihood calculations took 0.16 seconds.
#>  - Fitting model with 628 mixture components.
#>  - Model fitting took 0.18 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] T_Cell: 278 genes shrunk; sig (lfsr<0.05) = 21
#>  - Computing 278 x 404 likelihood matrix.
#>  - Likelihood calculations took 0.09 seconds.
#>  - Fitting model with 404 mixture components.
#>  - Model fitting took 0.36 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 278 x 590 likelihood matrix.
#>  - Likelihood calculations took 0.15 seconds.
#>  - Fitting model with 590 mixture components.
#>  - Model fitting took 0.60 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>   [mashr] Tumour: 278 genes shrunk; sig (lfsr<0.05) = 29
head(neighbourSlopes(fit))
#>     gene  focal neighbour   term      estimate   std.error estimate_shrunk
#> 1 ABCC11 B_Cell    B_Cell B_Cell -1.968614e-05 0.001376136               0
#> 2  ACTA2 B_Cell    B_Cell B_Cell -1.011682e-04 0.004245694               0
#> 3  ACTG2 B_Cell    B_Cell B_Cell  1.420043e-03 0.005170530               0
#> 4  ADAM9 B_Cell    B_Cell B_Cell  3.256417e-05 0.003072893               0
#> 5 ADGRE5 B_Cell    B_Cell B_Cell  8.103512e-04 0.005304599               0
#> 6  ADH1B B_Cell    B_Cell B_Cell -3.150005e-05 0.002131909               0
#>   sd_shrunk lfsr
#> 1         0    1
#> 2         0    1
#> 3         0    1
#> 4         0    1
#> 5         0    1
#> 6         0    1
```
