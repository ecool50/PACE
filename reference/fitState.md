# Raw solver state

The penalised quasi-likelihood solver's output: fixed and random effects
and their standard errors, gene-wise overdispersions, variance
components, the per-cell contamination loadings, and the convergence
history.

## Usage

``` r
fitState(object, ...)

# S4 method for class 'PACEFit'
fitState(object, ...)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- ...:

  Unused.

## Value

A named list of solver state.

## Details

This is the fit itself rather than a readout of it, and it is exposed
for inspection and diagnostics. Prefer
[`neighbourSlopes()`](https://ecool50.github.io/PACE/reference/neighbourSlopes.md),
[`varianceDecomposition()`](https://ecool50.github.io/PACE/reference/varianceDecomposition.md)
and
[`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md)
for interpretation: those apply the shrinkage and attribution the method
is defined by, and the raw BLUPs here do not.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
dim(fitState(fit)$U)
#> [1]  72 278
fitState(fit)$converged
#> [1] TRUE
```
