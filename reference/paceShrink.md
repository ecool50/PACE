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

  Further arguments passed to the shrinkage step, notably
  `null_correlation` (default `TRUE`: estimate the correlation between
  focal cell types under the null and pass it to mash as `V`; `FALSE`
  treats them as independent), `data_driven`, and `shrink_threads`.

  `shrink_threads` (default 4) is how many R processes shrink the
  neighbour slices at once. The slices are independent and carry about
  90% of the work, so this is where the time goes: on the full breast
  cancer cohort it takes the shrinkage from 55.4 to 21.1 seconds, and
  the called slopes agree with the serial path to 8.3e-13 with identical
  calls and no sign flips.

  It is not, however, the serial computation. With `data_driven = TRUE`
  the parallel path cannot reproduce the serial one: `cov_pca()` draws
  its starting vectors from the RNG, and a serial run consumes that
  stream slice by slice. Above 1 each slice is seeded with its own
  index, so the parallel result is reproducible run to run but differs
  from the serial one – by ~1e-12 at cohort scale, but by up to ~1e-6 on
  small slices, where there is less data to swamp the difference. Set
  `shrink_threads = 1` when you need the serial stream exactly, which is
  what the package's own fixture tests do. With `data_driven = FALSE`
  nothing draws from the stream and the two are identical.

## Value

The `PACEFit` with the shrunken neighbour slopes added.

## Details

One mash model is fitted per neighbour cell type, and every one of them
is fitted over the same genes. mash calibrates against the genes it is
given, so a common gene set keeps the lfsr comparable across neighbours;
a gene that cannot be used for one neighbour is dropped for all of them,
with a message.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
fit <- paceShrink(fit)
#>  - Computing 278 x 92 likelihood matrix.
#>  - Likelihood calculations took 0.03 seconds.
#>  - Fitting model with 92 mixture components.
#>  - Model fitting took 0.09 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 278 x 404 likelihood matrix.
#>  - Likelihood calculations took 0.18 seconds.
#>  - Fitting model with 404 mixture components.
#>  - Model fitting took 1.20 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 278 x 92 likelihood matrix.
#>  - Likelihood calculations took 0.01 seconds.
#>  - Fitting model with 92 mixture components.
#>  - Model fitting took 0.05 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.00 seconds.
#>  - Computing 278 x 313 likelihood matrix.
#>  - Likelihood calculations took 0.14 seconds.
#>  - Fitting model with 313 mixture components.
#>  - Model fitting took 3.60 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 278 x 391 likelihood matrix.
#>  - Likelihood calculations took 0.18 seconds.
#>  - Fitting model with 391 mixture components.
#>  - Model fitting took 2.20 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 278 x 417 likelihood matrix.
#>  - Likelihood calculations took 0.17 seconds.
#>  - Fitting model with 417 mixture components.
#>  - Model fitting took 5.90 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.09 seconds.
#>  - Computing 278 x 430 likelihood matrix.
#>  - Likelihood calculations took 0.20 seconds.
#>  - Fitting model with 430 mixture components.
#>  - Model fitting took 1.61 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.05 seconds.
#>  - Computing 278 x 628 likelihood matrix.
#>  - Likelihood calculations took 0.28 seconds.
#>  - Fitting model with 628 mixture components.
#>  - Model fitting took 1.28 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 278 x 404 likelihood matrix.
#>  - Likelihood calculations took 0.19 seconds.
#>  - Fitting model with 404 mixture components.
#>  - Model fitting took 5.90 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
#>  - Computing 278 x 590 likelihood matrix.
#>  - Likelihood calculations took 0.12 seconds.
#>  - Fitting model with 590 mixture components.
#>  - Model fitting took 0.55 seconds.
#>  - Computing posterior matrices.
#>  - Computation allocated took 0.01 seconds.
head(neighbourSlopes(fit))
#>     gene  focal neighbour   term      estimate   std.error estimate_shrunk
#> 1 ABCC11 B_Cell    B_Cell B_Cell -1.968907e-05 0.001376136               0
#> 2  ACTA2 B_Cell    B_Cell B_Cell -9.976965e-05 0.004245639               0
#> 3  ACTG2 B_Cell    B_Cell B_Cell  1.418293e-03 0.005170493               0
#> 4  ADAM9 B_Cell    B_Cell B_Cell  3.277669e-05 0.003072945               0
#> 5 ADGRE5 B_Cell    B_Cell B_Cell  8.125370e-04 0.005304314               0
#> 6  ADH1B B_Cell    B_Cell B_Cell -3.142864e-05 0.002131910               0
#>   sd_shrunk lfsr
#> 1         0    1
#> 2         0    1
#> 3         0    1
#> 4         0    1
#> 5         0    1
#> 6         0    1
```
