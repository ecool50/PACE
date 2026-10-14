# Stash PACE readouts in a SpatialExperiment

Copies a fit's interpretable output into `metadata(spe)`, so that
results travel with the data they came from and downstream code can find
them without being handed a second object.

## Usage

``` r
stashPACE(spe, fit, name = "PACE", overwrite = FALSE)
```

## Arguments

- spe:

  A
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  – the one the fit was built on. The cell count is checked against the
  fit.

- fit:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- name:

  The `metadata()` element to write. Defaults to `"PACE"`.

- overwrite:

  Whether to replace an existing element of that name.

## Value

`spe`, with `metadata(spe)[[name]]` populated.

## Details

The solver state is deliberately not copied. A fit is large, it is often
saved and reloaded independently of the data, and nothing downstream
reads it; see
[`fitState()`](https://ecool50.github.io/PACE/reference/fitState.md) if
you need it. What is stashed is the three reporting layers plus enough
provenance to tell which fit produced them: the cell types, the fitting
parameters, and the package version.

## See also

[`neighbourSlopes()`](https://ecool50.github.io/PACE/reference/neighbourSlopes.md),
[`varianceDecomposition()`](https://ecool50.github.io/PACE/reference/varianceDecomposition.md),
[`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md)

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
spe <- stashPACE(spe, fit)
names(S4Vectors::metadata(spe)$PACE)
#> [1] "neighbourSlopes"       "varianceDecomposition" "topDrivers"           
#> [4] "cellTypes"             "params"                "packageVersion"       
#> [7] "stashedAt"            
head(S4Vectors::metadata(spe)$PACE$neighbourSlopes, 3)
#>     gene  focal neighbour   term      estimate   std.error estimate_shrunk
#> 1 ABCC11 B_Cell    B_Cell B_Cell -1.968907e-05 0.001376136               0
#> 2  ACTA2 B_Cell    B_Cell B_Cell -9.976965e-05 0.004245639               0
#> 3  ACTG2 B_Cell    B_Cell B_Cell  1.418293e-03 0.005170493               0
#>   sd_shrunk lfsr
#> 1         0    1
#> 2         0    1
#> 3         0    1
```
