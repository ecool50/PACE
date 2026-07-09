# PACE: Proximity-Associated Changes in Expression

PACE quantifies how a cell's gene expression changes with proximity to
specific neighbouring cell types in imaging-based spatial
transcriptomics. It fits hierarchical negative binomial mixed models
with partial pooling across cell types, corrects for transcript
contamination between adjacent cells with a per-cell ambient term, and
decomposes expression variance into cell-type identity, spatial cell
state, contamination, and residual components.

## Details

The user-facing entry point is
[`paceFit()`](https://ecool50.github.io/PACE/reference/paceFit.md),
which takes a
[SpatialExperiment::SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
and returns a
[PACEFit](https://ecool50.github.io/PACE/reference/PACEFit-class.md)
object. The results are read out with
[`neighbourSlopes()`](https://ecool50.github.io/PACE/reference/neighbourSlopes.md),
[`varianceDecomposition()`](https://ecool50.github.io/PACE/reference/varianceDecomposition.md),
and
[`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md).

## See also

Useful links:

- <https://github.com/ecool50/PACE>

- <https://ecool50.github.io/PACE>

- Report bugs at <https://github.com/ecool50/PACE/issues>

## Author

**Maintainer**: Elijah Willie <ewillie0004@gmail.com>

Authors:

- Shreya Rajesh Rao

- John Ormerod

- Ellis Patrick
