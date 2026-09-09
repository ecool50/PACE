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
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
head(varianceDecomposition(fit))
#>          focal    gene Cell type %    Spatial % Spillover % Residual %
#> SEC11C  B_Cell  SEC11C  18.7147181 14.494678992   0.6373641   66.15324
#> DAPK3   B_Cell   DAPK3  14.5989691  0.002760858   1.0526045   84.34567
#> TCIM    B_Cell    TCIM  85.9791365  0.002251981   0.5945541   13.42406
#> NKG7    B_Cell    NKG7   0.3900736  0.002239134   8.2668799   91.34081
#> RAPGEF3 B_Cell RAPGEF3  14.2060602  0.008126525   3.2152159   82.57060
#> PPARG   B_Cell   PPARG   9.2327061  0.013567628   4.4774151   86.27631
#>           SS_lineage  SS_within      denom n_focal
#> SEC11C   638.5745510 2773.57704 3412.15159     330
#> DAPK3     14.4759924   84.68164   99.15763     330
#> TCIM    2153.2706156  351.13999 2504.41061     330
#> NKG7       0.3704028   94.58676   94.95716     330
#> RAPGEF3   37.1808640  224.54451  261.72537     330
#> PPARG     31.1183156  305.92605  337.04436     330
```
