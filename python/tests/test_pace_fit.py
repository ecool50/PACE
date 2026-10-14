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
