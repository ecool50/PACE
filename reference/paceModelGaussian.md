# Fit the PACE model to continuous intensities

The Gaussian, identity-link form of
[`paceModel()`](https://ecool50.github.io/PACE/reference/paceModel.md),
for continuous measurements such as imaging mass cytometry protein
intensities rather than counts.

## Usage

``` r
paceModelGaussian(object, ...)

# S4 method for class 'SpatialExperiment'
paceModelGaussian(
  object,
  celltype_col,
  image_col = NULL,
  condition_col = NULL,
  assay_name = "intensity",
  h_bio = 30,
  h_tech = 5,
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
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- ...:

  Passed to
  [`pace_fit_streaming_gaussian()`](https://ecool50.github.io/PACE/reference/pace_fit_streaming_gaussian.md).

- celltype_col:

  Column of `colData` holding the cell type labels.

- image_col:

  Column of `colData` identifying the image. Defaults to the single
  image when there is only one.

- condition_col:

  Optional column of `colData` giving a condition contrast.

- assay_name:

  The assay holding the intensities.

- h_bio, h_tech:

  Bandwidths of the biological and technical kernels.

- kernel_per_image:

  Build the kernels within each image separately.

- image_re:

  An optional second random-effect block over images.

- types:

  Cell types, in the order the fit should use them.

- resp_term:

  Optional response term passed to the model frame.

- verbose:

  Print the fitting trace.

## Value

A [PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md).

## Details

The working response on the identity link is the intensity itself and
the working weight is the per-gene residual variance, so the mean
converges in a single inner solve and the outer iteration only moves the
variance components.

This path does not model technical contamination. Correct the
intensities first, against an ambient field built with
[`ambientField()`](https://ecool50.github.io/PACE/reference/ambientField.md),
and pass the corrected matrix here.
