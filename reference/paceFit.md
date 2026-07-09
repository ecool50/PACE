# Fit a PACE model to a SpatialExperiment

`paceFit()` runs the full pipeline in one call:
[`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md)
-\>
[`paceShrink()`](https://ecool50.github.io/PACE/reference/paceShrink.md)
-\>
[`paceDecompose()`](https://ecool50.github.io/PACE/reference/paceDecompose.md)
-\>
[`paceDrivers()`](https://ecool50.github.io/PACE/reference/paceDrivers.md).
For step-by-step control (inspecting the fitted model, re-running the
downstream without refitting), call those stages directly.

## Usage

``` r
paceFit(object, ...)

# S4 method for class 'SpatialExperiment'
paceFit(object, ..., pairs = NULL)
```

## Arguments

- object:

  A
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- ...:

  Arguments passed to
  [`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md)
  (`celltype_col`, `image_col`, `contamination`, `dispersion`, ...).

- pairs:

  Optional list of focal-neighbour pairs for the driver tables.

## Value

A fully populated
[PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- paceFit(spe, celltype_col = "cellType")
neighbourSlopes(fit)
} # }
```
