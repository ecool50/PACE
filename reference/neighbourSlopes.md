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
if (FALSE)  neighbourSlopes(fit)  # \dontrun{}
```
