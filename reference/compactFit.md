# Drop the per-cell matrices from a fit

A fit stores the fitted means and the contamination offsets as dense
cells by genes matrices. On a real cohort those dominate its size, and
nothing downstream needs them kept:
[`paceDecompose()`](https://ecool50.github.io/PACE/reference/paceDecompose.md)
and
[`paceDrivers()`](https://ecool50.github.io/PACE/reference/paceDrivers.md)
read per-cell-type summaries the fit also holds, and rebuild the
matrices from the counts when they are absent.

## Usage

``` r
compactFit(
  object,
  which = c("mu", "mu_bio", "mu_spill", "technical_offset_mat", "bleed_offset_mat")
)
```

## Arguments

- object:

  A
  [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

- which:

  The matrices to drop. The default is all of them.

## Value

`object`, without the named matrices.

## Details

Use this before saving a fit. Every readout still works afterwards, and
re-decomposing a compacted fit reproduces the tables it already carries.

## See also

[`stashPACE()`](https://ecool50.github.io/PACE/reference/stashPACE.md)
for keeping the readouts with the data.

## Examples

``` r
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
small <- compactFit(fit)
c(before = object.size(fit), after = object.size(small)) |> format(units = "MB")
#>    before     after 
#> "7086784" "7086784" 
```
