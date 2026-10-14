# PACE for Python

Bindings for the PACE compiled core. Nothing here reimplements the method: the
extension compiles `../src/core/*.cpp` directly, the same translation units the
R package uses, so the two cannot drift apart.

## Install

```sh
pip install -e python/
```

### macOS: pin the toolchain

A conda or Homebrew `clang` usually brings an old `ld64` along with it, and an
old `ld64` links a shared object that loads nowhere:

```
ImportError: dlopen(...): mis-aligned LINKEDIT string pool
```

The build succeeds, so this only shows up at import. Use Apple's toolchain:

```sh
CC=/usr/bin/clang CXX=/usr/bin/clang++ \
  LDSHARED="/usr/bin/clang -bundle -undefined dynamic_lookup" \
  pip install -e python/
```

`ld -v` should report `PROJECT:ld-27037` or newer. `ld64-530` is the one that
fails. The R package has the same requirement for the same reason.

### Pin Eigen to 3.4

The core's solver uses Eigen, which the R package gets from RcppEigen. A Python
build needs it separately, and the version is not a free choice: against Eigen
3.5 the random effects move by about 3e-4, so the two languages no longer agree.
The build warns when it picks up anything but 3.4.

Wherever R is installed, RcppEigen is already a 3.4 tree:

```sh
EIGEN_INCLUDE_DIR=$(Rscript -e 'cat(system.file("include", package="RcppEigen"))') \
  pip install -e python/
```

## What is bound so far

| function | does |
|---|---|
| `mash` | canonical multivariate adaptive shrinkage for one neighbour slice |
| `driver_scores` | per-gene driver scores for one focal-neighbour pair |
| `pair_variance_pratt` | Pratt attribution of spatial variance across neighbours |
| `neighbour_kernels` | Gaussian and technical neighbour kernels from coordinates |
| `ambient_field` | sparse cross-type contamination field, compressed-column |
| `random_design_block` | one block of the random-effect design, plus the grouping |
| `fit_irls` | the penalised quasi-likelihood loop: the fit itself |
| `sparse_product` | sparse-sparse product, in the core's summation order |
| `pace_fit` | every stage above, in order: counts and coordinates to a fit |
| `fit_anndata` | the same, reading an `AnnData` |
| `pace_shrink` | multivariate adaptive shrinkage of the fitted slopes |
| `pace_decompose` | the per-gene, per-cell-type variance decomposition |
| `single_frame_decomposition` | the observed decomposition R reports |
| `pace_drivers` | per-pair driver tables, ranked |
| `fit_statistics` | the per-cell-type summaries the readouts are built on |

`pace_fit()` runs the whole pipeline from counts, coordinates and cell type
labels, and reproduces `paceModel()` bit for bit: every one of the 20,016
random effects on the packaged breast cancer subset, and every dispersion.
The caveat is the trigamma function, which R takes from its own math library
and the core computes itself; the two agree to 1.8e-15, which thirteen IRLS
iterations amplify to 8e-07 in the fitted effects. Against an R build using the
core's trigamma, nothing differs at all.

```python
import pace

fit = pace.fit_anndata(adata, "cell_type")             # from an AnnData
fit = pace.pace_fit(counts, x, y, cell_type)           # or from arrays
fit.U[fit.row_names.index("Macrophage::Tumour")]
```

`fit_anndata()` reads the coordinates from `obsm["spatial"]`, the labels from
`obs`, the gene names from `var_names`, and the counts from `X` or a `layer`.
It never modifies the object, and `anndata` is not a dependency: it uses only
the attributes an `AnnData` exposes, so anything shaped like one works.

It checks that the counts are counts, which is the one thing the R package
never has to: `counts(spe)` names the assay it wants, whereas `adata.X` holds
whatever the last preprocessing step left there. A negative binomial fit of
log-normalised values converges and reports slopes that mean nothing, so
non-integer values are refused rather than fitted. Point `layer` at the raw
counts, usually `layer="counts"`.

It also normalises a sparse matrix that is not in canonical form. scipy will
hand over a CSC holding duplicate entries for one `(cell, gene)`, and the
solver reads stored entries one at a time, so those count twice: such a matrix
can be elementwise identical to a good one, hold none but non-negative
integers, give the same library sizes, and still fit 30% wrong. R is safe from
this only because `Matrix` refuses the object outright. Prefer a sparse `X`:
a dense one is converted here, at about six times its own size at peak.

Two things to know before comparing a Python fit to an R one row by row:

**Cell type order.** `types=None` sorts the observed labels, which is Python's
byte order. R sorts under `LC_COLLATE`, so a label set mixing cases orders
differently: R puts `epithelial` before `T-NK`, Python after. The fit is
equivariant under a permutation of types, so the biology is unchanged, but the
row order is not. Pass `types=` to match a given R fit, as you would in R.

**Neighbour scope.** `kernel_per_image` defaults to `False`, as R's does: even
with several images the neighbour search is global, which is what a cohort of
fields of view from one section wants. Set it `True` when the images are
separate samples. The contamination field is per-image either way.

**Condition cohorts.** Name the reference arm. Unlike the cell type order,
which only relabels the output, the condition's first level is what every other
arm is contrasted against: it sets what the baseline neighbour slopes mean and
which way the interaction points. R chooses it by locale collation and Python
would choose it by byte order, and the two disagree on labels as ordinary as
`{"Responder", "non-responder"}` -- where the fitted slopes then differ by 2
rather than being reordered. Leaving `condition_levels` unset warns for exactly
that reason.

A numeric or boolean condition is refused: R would model such a column
continuously, where treatment contrasts give a dummy per level, and the two are
different fits.

Pass `condition_key` (or `condition=` and
`condition_name=` for arrays) and the model gains a condition main effect in
the fixed design and a condition-by-neighbour interaction in the random-effect
block, which is the responder layer: `~ 1 + condition * (types)` over
`~ 1 + condition`, as in R. The terms take their prefix from the column name,
so a column `Responder` with levels PD and PR gives `ResponderPR` and
`ResponderPR:Tumour`, and a Categorical's category order is honoured as R
honours a factor's. Checked against R on a two-level and a three-level
condition; the three-level case is the one that pins the column order, since
`model.matrix()` varies the condition level fastest within each type.

The ambient field is always materialised, with no equivalent of R's
`ambient_mode = "stream"`, so a very large cohort needs the memory for the
whole cells-by-genes product.

## Shrinking the slopes

A fitted slope on its own is mostly noise. `pace_shrink()` pools them the way
R's `paceShrink()` does: one mash fit per neighbour, with the focal cell types
as its conditions, over a single gene set shared by every slice so that one
lfsr cut means the same thing across them.

```python
slopes = pace.pace_shrink(fit, null_correlation=False)
slopes.calls()                      # how many effects, and how many wrong signs
frame = slopes.to_dataframe()       # pandas, imported only if you ask for this
frame.nsmallest(10, "lfsr")[["gene", "focal", "neighbour", "estimate_shrunk"]]
```

On the packaged subset this calls the **same gene-focal pairs** as R, not merely
the same number of them, under both `null_correlation` settings. Three things
to know:

- The mixture is the compiled core's, which reproduces mashr's canonical
  covariances but not the data-driven ones, so this is R at
  `data_driven = FALSE`. On the two canonical cohorts that costs 11 of 1,637
  breast cancer calls and none of 46 melanoma ones.
- `null_correlation` defaults to `True`, as R's does. The published breast
  cancer counts were produced with `False`, so pass `False` to match them.
- On a condition fit the slopes shrunk are `<resp_term>:<neighbour>`, so what
  comes back is the **responder layer and not the baseline proximity effect**.
  That is R's behaviour, and `resp_term` defaults to the condition's last level
  either way.

R's table arrives in alphabetical gene order, an artifact of its splitting the
long frame by gene; Python keeps the fit's gene order. Sort if you are comparing
row by row.

## Decomposing the variance

`pace_decompose()` splits each gene's fitted variance, for each focal cell
type, into five blocks: the cell type's own offset, the baseline spatial state
(the neighbour slopes), the responder spatial state (their interaction with a
condition), the contamination, and what the dispersion accounts for.

```python
decomposition = pace.pace_decompose(fit, counts)
decomposition.aggregates("specw")        # per focal type, specificity-weighted
decomposition.to_dataframe()             # per (focal, gene)
```

It matches R's five-block table and all three of its aggregate weightings. The
fitted means it needs are not stored on a fit -- on a real cohort they are a
dense cells-by-genes matrix -- so they are rebuilt a gene block at a time, which
is what R does too; `fit_statistics()` exposes that step if you want it
separately.

`single_frame_decomposition()` gives the other table, the one R's
`varianceDecomposition()` actually returns:

```python
frame = pace.single_frame_decomposition(fit, counts)
frame.by_focal()            # one row per focal type, sum-of-squares weighted
```

Here the magnitudes are **observed** rather than fitted -- the cell-type share
is the between-type sum of squares of `log1p(CP10k)` -- and only the split of
what remains inside a cell type comes from the model. That is why its
percentages differ from the five-block ones, and why its observed columns match
R to 1e-11 while the fitted split matches to the trigamma floor. The five-block
table above is what sits underneath it, reached in R as
`fit@varianceDecomposition$blocks$gene_focal_5block`.


`ambient_field()` returns compressed-column arrays rather than a matrix, so
nothing here depends on scipy:

```python
from scipy.sparse import csc_matrix
f = pace.ambient_field(x, y, cell_type, image, n_images=1)
W = csc_matrix((f["data"], f["indices"], f["indptr"]), shape=f["shape"])
```
The other three functions take fitted output and summarise it.

## Tests

```sh
python -m pytest python/tests
```

The mash tests hold the bindings to the same target as the R package: every
neighbour slice of both canonical cohorts, compared gene by gene against mashr
rather than by count, giving 1634 breast cancer and 426 melanoma calls.

## Ranking the drivers of a pair

`pace_drivers()` ranks, for each focal-neighbour pair, the genes mediating it:

```python
slopes = pace.pace_shrink(fit, null_correlation=False)
stats = pace.fit_statistics(fit, counts)
decomposition = pace.pace_decompose(fit, counts, stats=stats)
drivers = pace.pace_drivers(fit, slopes, decomposition, stats)

drivers["Macrophage", "Tumour"].to_dataframe().head()
```

The score is `b^2 * spec^2 * focal_mean`, with `b` the shrunken slope. Squaring
the slope makes it direction-blind, so a gene pushed down near the neighbour
ranks alongside one pushed up; the specificity and level keep a gene that is
barely expressed in the focal type off the top of the table on a noisy slope.
The filter is lfsr < 0.05 alone -- there is no specificity cut at this layer,
since `spec` is already in the score.

Against R this gives the same pairs, the same genes in the same rank order, and
all fourteen columns to 2e-7. Each pair also carries its status and how many of
its calls are expected to have the wrong sign.

## The whole pipeline

```python
fit    = pace.fit_anndata(adata, "cell_type")
slopes = pace.pace_shrink(fit, null_correlation=False)
stats  = pace.fit_statistics(fit, counts)
blocks = pace.pace_decompose(fit, counts, stats=stats)
frame  = pace.single_frame_decomposition(fit, counts, blocks)
drive  = pace.pace_drivers(fit, slopes, blocks, stats)
```
