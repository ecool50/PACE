# Per-gene variance decomposition

Returns the per-gene, per-focal variance decomposition. By default the
observed single-frame decomposition (cell-type identity, spatial cell
state, contamination, residual; and a responder spatial block for
condition cohorts) is returned; set `which = "blocks"` for the
underlying link-scale block tables.

## Usage

``` r
varianceDecomposition(object, ...)

# S4 method for class 'PACEFit'
varianceDecomposition(object, which = c("perGene", "blocks"), ...)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- ...:

  Unused.

- which:

  Either `"perGene"` (default, the observed single-frame table) or
  `"blocks"` (the raw link-scale decomposition).

## Value

A data frame (for `"perGene"`) or a list (for `"blocks"`).

## Examples

``` r
if (FALSE)  varianceDecomposition(fit)  # \dontrun{}
```
