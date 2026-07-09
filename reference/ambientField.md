# Build the ambient contamination field

Returns the sparse cross-cell-type ambient weight matrix that defines
the per-cell contamination term (the short-range technical field a cell
receives from its heterotypic neighbours).

## Usage

``` r
ambientField(object, ...)

# S4 method for class 'SpatialExperiment'
ambientField(
  object,
  celltype_col,
  image_col = NULL,
  h_tech = 5,
  assay_name = "counts",
  ...
)
```

## Arguments

- object:

  A
  [SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- ...:

  Passed to the ambient-field builder.

- celltype_col:

  colData column with the cell-type annotation.

- image_col:

  Optional colData column grouping cells into images.

- h_tech:

  Technical bandwidth in micrometres (default 5).

- assay_name:

  Counts assay name (default `"counts"`).

## Value

The ambient-field object (sparse weight matrix and image index).
