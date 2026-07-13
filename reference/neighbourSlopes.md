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
