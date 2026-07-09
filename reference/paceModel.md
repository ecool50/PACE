# Fit the PACE model

`paceModel()` builds the biological and technical neighbourhood kernels
from the spatial coordinates and cell-type labels and fits the
hierarchical negative binomial mixed model by streaming penalised
quasi-likelihood (with the per-cell contamination correction). It
returns a
[PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
carrying only the fitted model; the reporting layers are added by
[`paceShrink()`](https://ecool50.github.io/PACE/reference/paceShrink.md),
[`paceDecompose()`](https://ecool50.github.io/PACE/reference/paceDecompose.md),
and
[`paceDrivers()`](https://ecool50.github.io/PACE/reference/paceDrivers.md)
(or all at once by
[`paceFit()`](https://ecool50.github.io/PACE/reference/paceFit.md)).

## Usage

``` r
paceModel(object, ...)

# S4 method for class 'SpatialExperiment'
paceModel(
  object,
  celltype_col,
  image_col = NULL,
  condition_col = NULL,
  assay_name = "counts",
  h_bio = 30,
  h_tech = 5,
  contamination = c("percell_hc", "none"),
  dispersion = c("nb1", "nb2"),
  kernel_per_image = FALSE,
  image_re = c("none", "intercept", "slopes", "condition_slopes"),
  types = NULL,
  resp_term = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  with a counts assay and two-dimensional spatial coordinates.

- ...:

  Further arguments passed to the streaming fitter (`n_iter`, `threads`,
  `drop_sparse_neff`, `within_image`, `edge_correct`,
  `data_informed_tau`, `tau_shrinkage`, ...).

- celltype_col:

  colData column with the discrete cell-type annotation.

- image_col:

  Optional colData column grouping cells into images/samples. `NULL`
  (default) treats the object as one section.

- condition_col:

  Optional colData column with a binary condition; enables the responder
  spatial block for multi-sample cohorts.

- assay_name:

  Counts assay name (default `"counts"`).

- h_bio, h_tech:

  Biological and technical kernel bandwidths in micrometres (defaults 30
  and 5).

- contamination:

  `"percell_hc"` (default) or `"none"`.

- dispersion:

  `"nb1"` (default) or `"nb2"`.

- kernel_per_image:

  If `TRUE`, neighbourhoods are built within each image.

- image_re:

  Second image-level random-effect block: `"none"`, `"intercept"`,
  `"slopes"`, or `"condition_slopes"`.

- types:

  Cell types in the order they should index the random-effect blocks.
  `NULL` (default) uses every observed type in alphabetical order; pass
  an explicit vector to reproduce a locked fit, whose block order the
  caller chose.

- resp_term:

  Responder interaction term prefix for condition cohorts; if `NULL` it
  is derived from `condition_col`.

- verbose:

  Whether to print progress.

## Value

A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
with the fitted model (reporting layers empty).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- paceModel(spe, celltype_col = "cellType") |>
         paceShrink() |> paceDecompose(spe) |> paceDrivers()
} # }
```
