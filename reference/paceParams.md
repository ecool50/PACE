# Fitting parameters

The settings the fit actually used, rather than the ones requested:
kernel bandwidths, the contamination and dispersion models, the column
names read from the `SpatialExperiment`, and the condition term if there
is one. Worth reading before comparing two fits, since several of these
change the answer.

## Usage

``` r
paceParams(object, ...)

# S4 method for class 'PACEFit'
paceParams(object, ...)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- ...:

  Unused.

## Value

A named list of parameters.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
paceParams(fit)$h_bio
#> [1] 30
```
