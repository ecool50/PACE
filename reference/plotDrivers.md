# Per-pair driver composite

Reproduces the manuscript per-pair driver figure: the top genes ranked
by driver score (MCSD) for a focal-neighbour pair, alongside their
per-gene single-frame variance decomposition.

## Usage

``` r
plotDrivers(
  object,
  focal,
  neighbour,
  n_top = 5,
  panels = c("both", "mcsd", "gene")
)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- focal, neighbour:

  The focal and neighbour cell types.

- n_top:

  Number of top genes to show.

- panels:

  One of `"both"`, `"mcsd"`, or `"gene"`.

## Value

A `patchwork` / `ggplot` object.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
plotDrivers(fit, "Macrophage", "Tumour")
```
