# Shrunken neighbour proximity slopes

Returns the shrunken per-(gene, focal, neighbour) proximity coefficients
with their local false sign rate (lfsr). For condition cohorts the
responder interaction terms are included.

## Usage

``` r
neighbourSlopes(object, ...)

# S4 method for class 'PACEFit'
neighbourSlopes(object, ...)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- ...:

  Unused.

## Value

A data frame of shrunken slopes.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
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
