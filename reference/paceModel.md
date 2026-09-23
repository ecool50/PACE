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
  `chunk_size`, `drop_sparse_neff`, `within_image`, `edge_correct`,
  `data_informed_tau`, `tau_shrinkage`, and the memory and approximation
  settings described below).

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

## Memory and approximation settings

Three pass-through arguments change how much memory a fit needs, or
trade exactness for time. Two of them are approximations that are ON by
default.

- `ambient_mode`:

  `"cache"` (default) or `"stream"`. The contamination model needs an
  ambient field per cell and gene; `"cache"` materialises that n by G
  product, `"stream"` recomputes each chunk's columns from the weights
  and the counts. The numbers are identical. On a 1.2M cell, 5,001 gene
  panel the cached product alone is over 5 GB, so `"stream"` is what
  makes a full transcriptome panel fit in memory at all.

- `alpha_warmup`:

  Default `6`. The per-gene dispersion MLE is re-fitted only on the
  first `alpha_warmup` iterations and on the last one; in between, alpha
  is frozen at its warmed-up value. Alpha typically settles within five
  iterations while the MLE is a large share of each iteration's cost, so
  this is on by default. **It is an approximation**: pass `Inf` to
  re-fit the dispersion on every iteration. Lowering it below the
  default does move results – on the breast cancer cohort `4` changes
  the number of calls at `lfsr < 0.05`.

- `alpha_max_n`:

  Default `Inf`, meaning the dispersion MLE sees every cell. A finite
  value caps it at an even, deterministic subsample; the estimator, the
  Brent search and its tolerance are unchanged. The dispersion is one
  scalar per gene and its standard error falls as `1/sqrt(n)`, so most
  cells add little. **Validate any cap on your own cohort before
  trusting it.** A cap of 50,000 reproduced two cohorts of roughly 10^5
  cells exactly, and on a 1.2M cell panel – where the same cap is a far
  smaller fraction of the data – it left the number of calls unchanged
  while swapping the identity of 72 of them. The binding quantity is not
  the absolute subsample.

## Examples

``` r
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
# \donttest{
fit <- paceModel(spe, celltype_col = "cellType", verbose = FALSE)
fit
#> class: PACEFit
#> cell types (8): B_Cell, Dendritic_Cell, Endothelial, Macrophage, Myoepithelial, Stromal, T_Cell, Tumour
#> kernels: h_bio = 30 um, h_tech = 5 um | contamination: percell_hc; dispersion: nb1
#> pipeline: model
# }
```
