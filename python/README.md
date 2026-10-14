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

`pace_fit()` runs the whole pipeline from counts, coordinates and cell type
labels, and reproduces the R package's fit bit for bit: every one of the 20,016
random effects on the packaged breast cancer subset, and every dispersion.
The one caveat is the trigamma function, which R takes from its own math
library and the core computes itself; the two agree to 1.8e-15, which thirteen
IRLS iterations amplify to 8e-07 in the fitted effects.

A condition term is not implemented. It needs R's `model.matrix` semantics,
and approximating those would make the two disagree quietly, so `condition`
raises instead.

```python
import pace
fit = pace.pace_fit(counts, x, y, cell_type)   # counts is cells by genes
fit.U[fit.row_names.index("Macrophage::Tumour")]
```

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
