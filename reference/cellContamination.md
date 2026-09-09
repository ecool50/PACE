# Per-cell contamination

Returns the per-cell contamination loading `rho_i` and the fraction of a
cell's expected counts attributed to its local ambient field. Under the
per-cell contamination model the expected count is
`mu_ig = mu_bio_ig + rho_i * a_ig`, where `a_ig` is the cross-cell-type
ambient field at cell `i` and `rho_i` is a single empirical-Bayes
shrunken loading shared across genes. The contamination fraction
summarises that second term over genes,
`contamFraction_i = sum_g mu_spill_ig / sum_g mu_ig`, and is the
quantity reported as `contam_frac` in the solver's fitting trace.

## Usage

``` r
cellContamination(object, ...)

# S4 method for class 'PACEFit'
cellContamination(object, ...)
```

## Arguments

- object:

  A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
  fitted with `contamination = "percell_hc"`.

- ...:

  Unused.

## Value

A data frame with one row per cell, in fitting order: `cell` (the
working-frame cell identifier, or the row index if the frame is
unnamed), `celltype` (the label the cell was fitted under), `rho` (the
contamination loading) and `contamFraction`.

## Details

A high contamination fraction marks a cell whose profile is
substantially explained by its neighbours rather than by its own cell
type, which is the expected signature of a segmentation or
transcript-assignment error. It is not a complete cell-quality score:
because `rho_i` carries a fixed gene-direction (the local ambient), it
will under-report contamination whose direction departs from the
neighbourhood average, and a cell that is wholly a segmentation artefact
may fit some incorrect cell type with a low `rho_i`. Pair it with
ordinary per-cell quality control.

Read `contamFraction` rather than `rho` on its own. The ambient field is
cross-cell-type, so a cell with no differently-typed neighbour inside
the technical kernel has `a_ig = 0` for every gene and therefore no
contamination at all. Its `rho_i` is then unidentified and shrinks to
the empirical Bayes prior mean, giving every such cell the same
apparently middling loading. On the shipped example fit this is 4,589 of
7,898 cells, 84% of the tumour cells in a tumour-dominated crop.
`contamFraction` reports them as zero, which is correct; `rho` alone
does not.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
cc <- cellContamination(fit)
head(cc)
#>   cell   celltype         rho contamFraction
#> 1  442     Tumour 0.115044201    0.000000000
#> 2  444     B_Cell 0.002827062    0.001557514
#> 3  446 Macrophage 0.587792638    0.123800294
#> 4  449     Tumour 0.115044201    0.000000000
#> 5  451 Macrophage 0.089269111    0.102356337
#> 6  455     Tumour 0.115044201    0.000000000
# Contamination is highest where cells sit against a different type.
tapply(cc$contamFraction, cc$celltype, median)
#>         B_Cell Dendritic_Cell    Endothelial     Macrophage  Myoepithelial 
#>     0.09101894     0.08295640     0.01808932     0.07264364     0.03763919 
#>        Stromal         T_Cell         Tumour 
#>     0.03767385     0.07218389     0.00000000 
```
