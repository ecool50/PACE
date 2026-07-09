# Per-focal variance decomposition plot

Reproduces the manuscript per-focal decomposition figure: a full stacked
bar of the cell-type identity, spatial cell state, and contamination
(spillover) components, with a zoomed grouped bar of the two small
spatial components.

## Usage

``` r
plotDecomposition(object, title = "Per-focal decomposition")
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- title:

  Plot title.

## Value

A `patchwork` / `ggplot` object.

## Examples

``` r
if (FALSE)  plotDecomposition(fit)  # \dontrun{}
```
