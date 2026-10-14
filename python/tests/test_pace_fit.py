"""`pace_fit()` end to end, against the fit R produces from the same raw data.

test_fit.py checks the solver when it is handed R's own captured inputs. This
checks the layer above: starting from nothing but counts, coordinates and cell
type labels, Python builds every design, kernel, offset, anchor mask and prior
weight itself. Agreement here means the orchestration matches R's as well as
the numerics do.

The reference is an R build using the core's trigamma, because that is the one
remaining difference between the two languages. Against R's stock trigamma the
fit still agrees, but only to the tolerance test_fit.py documents.
"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

from pace.fit import pace_fit

FIX = Path(__file__).parent / "fixtures"


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def lines(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


@pytest.fixture(scope="module")
def fit():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    return pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                    lines("raw_celltype.txt"), n_threads=1)


def test_the_random_effects_are_bit_identical_to_r(fit):
    reference = vec("irls_ab_U.txt").reshape(fit.U.shape, order="F")
    assert np.array_equal(fit.U, reference)


def test_the_dispersions_are_bit_identical_to_r(fit):
    assert np.array_equal(fit.alpha, vec("irls_ab_alpha.txt"))


def test_the_fit_agrees_with_rs_stock_trigamma_to_the_documented_tolerance(fit):
    """R's released build uses Rf_trigamma, the core uses its own.

    Thirteen IRLS iterations amplify a 1.8e-15 difference in the trigamma into
    this, so the tolerance is the one test_fit.py justifies, not a guess.
    """
    np.testing.assert_allclose(fit.U, vec("raw_expect_U.txt").reshape(fit.U.shape, order="F"),
                               rtol=0, atol=1e-06)
    np.testing.assert_allclose(fit.alpha, vec("raw_expect_alpha.txt"), rtol=0, atol=1e-05)


def test_the_random_effect_rows_are_named_the_way_r_names_them(fit):
    assert fit.row_names == lines("raw_expect_rownames.txt")


def test_the_fit_converges_the_way_r_did(fit):
    n_iter, converged = vec("raw_expect_meta.txt", int)
    assert fit.n_iter == int(n_iter)
    assert fit.converged is bool(converged)


def test_the_contamination_loadings_match_rs(fit):
    """One loading per cell, against R's `percell_bleed_rho`.

    R's reference here is its stock build, so this carries the same trigamma
    difference as the fit itself: 1.1e-06 absolute on loadings that run from
    1e-04 to 3.7.
    """
    reference = vec("raw_expect_rho.txt")
    assert fit.rho.size == reference.size
    np.testing.assert_allclose(fit.rho, reference, rtol=1e-05, atol=0)


def test_a_condition_term_refuses_rather_than_approximating():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    with pytest.raises(NotImplementedError, match="model.matrix"):
        pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"), lines("raw_celltype.txt"),
                 condition=np.zeros(int(dim[0])))


# --- several images -------------------------------------------------------
# The packaged subset is one image, so it cannot tell a global neighbour search
# from a per-image one. These split it on the median x and hold both settings to
# R. The two R fits differ by 0.167 in U, so the pair discriminates: taking the
# wrong branch is not a rounding error. Tolerances are the trigamma floor.

@pytest.fixture(scope="module")
def two_image_inputs():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    return counts, vec("raw_x.txt"), vec("raw_y.txt"), lines("raw_celltype.txt"), \
        lines("two_image_label.txt")


@pytest.mark.parametrize("kernel_per_image,tag", [(False, "global"), (True, "perimage")])
def test_both_neighbour_scopes_match_r_on_two_images(two_image_inputs, kernel_per_image, tag):
    counts, x, y, cell_type, image = two_image_inputs
    fit = pace_fit(counts, x, y, cell_type, image=image,
                   kernel_per_image=kernel_per_image, n_threads=1)
    reference = vec(f"two_image_{tag}_U.txt").reshape(fit.U.shape, order="F")
    np.testing.assert_allclose(fit.U, reference, rtol=0, atol=1e-06)
    np.testing.assert_allclose(fit.alpha, vec(f"two_image_{tag}_alpha.txt"), rtol=0, atol=1e-05)


def test_the_default_neighbour_scope_is_global_like_rs(two_image_inputs):
    """R's `kernel_per_image` defaults to FALSE even with several images."""
    counts, x, y, cell_type, image = two_image_inputs
    fit = pace_fit(counts, x, y, cell_type, image=image, n_threads=1)
    reference = vec("two_image_global_U.txt").reshape(fit.U.shape, order="F")
    np.testing.assert_allclose(fit.U, reference, rtol=0, atol=1e-06)
    assert fit.params["kernel_per_image"] is False


# --- what must be refused rather than guessed ----------------------------

def test_a_missing_cell_type_is_refused_not_turned_into_a_type(two_image_inputs):
    """str(nan) is "nan", which would otherwise be a cell type like any other."""
    counts, x, y, cell_type, _ = two_image_inputs
    labels = list(cell_type)
    labels[7] = np.nan
    with pytest.raises(ValueError, match="every cell needs a cell type"):
        pace_fit(counts, x, y, labels)


def test_an_explicit_type_order_is_honoured(two_image_inputs):
    """R's `types` argument, which is how a locked row order is reproduced."""
    counts, x, y, cell_type, _ = two_image_inputs
    reversed_types = sorted(set(cell_type))[::-1]
    fit = pace_fit(counts, x, y, cell_type, types=reversed_types, n_threads=1)
    assert fit.cell_types == reversed_types
    assert fit.row_names[0] == f"{reversed_types[0]}::(Intercept)"


def test_an_unknown_label_is_refused(two_image_inputs):
    counts, x, y, cell_type, _ = two_image_inputs
    with pytest.raises(ValueError, match="not in `types`"):
        pace_fit(counts, x, y, cell_type, types=sorted(set(cell_type))[1:])
