"""The fitted-mean rebuild and the per-cell-type statistics, against R.

A fit keeps no fitted means, so every readout starts by rebuilding them. R does
the same, and the two have to agree before anything built on top of them can.
Tolerances are the fit's own trigamma difference carried through.
"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

import pace

FIX = Path(__file__).parent / "fixtures"


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def lines(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


@pytest.fixture(scope="module")
def counts():
    dim = vec("raw_y_dim.txt", int)
    return csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))


@pytest.fixture(scope="module")
def fit(counts):
    return pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                         lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                         n_threads=1)


def test_the_rebuilt_biological_means_match_rs(fit, counts):
    mu_bio, _ = pace.fitted_means(fit, counts, [0, 1, 139])
    reference = vec("mu_bio_g3.txt").reshape(mu_bio.shape, order="F")
    np.testing.assert_allclose(mu_bio, reference, rtol=1e-07, atol=0)


def test_the_rebuilt_contamination_means_match_rs(fit, counts):
    _, mu_spill = pace.fitted_means(fit, counts, [0, 1, 139])
    reference = vec("mu_spill_g3.txt").reshape(mu_spill.shape, order="F")
    np.testing.assert_allclose(mu_spill, reference, rtol=1e-04, atol=1e-06)


def test_the_statistics_match_rs(fit, counts):
    stats = pace.fit_statistics(fit, counts)
    dim = vec("stats_mu_mean_dim.txt", int)
    assert stats.cell_types == lines("stats_types.txt")
    assert np.array_equal(stats.n_cells, vec("stats_n_cells.txt", int))
    assert stats.has_contamination
    np.testing.assert_allclose(
        stats.mu_mean, vec("stats_mu_mean.txt").reshape(int(dim[0]), int(dim[1]), order="F"),
        rtol=1e-06, atol=1e-06)
    np.testing.assert_allclose(
        stats.toff_var, vec("stats_toff_var.txt").reshape(int(dim[0]), int(dim[1]), order="F"),
        rtol=1e-06, atol=1e-06)


def test_blocking_the_genes_does_not_change_the_answer(fit, counts):
    """The block size is a memory knob, so it must not touch the numbers."""
    whole = pace.fit_statistics(fit, counts, max_block_entries=10 ** 9)
    blocked = pace.fit_statistics(fit, counts, max_block_entries=7898 * 3)
    assert np.array_equal(whole.mu_mean, blocked.mu_mean)
    assert np.array_equal(whole.toff_var, blocked.toff_var)


def test_the_means_are_floored_the_way_the_solver_floors_them(fit, counts):
    mu_bio, mu_spill = pace.fitted_means(fit, counts, [0, 1, 139])
    assert mu_bio.min() >= 1e-6
    assert mu_spill.min() >= 0.0


def test_a_fit_without_the_readout_inputs_says_so(fit, counts):
    import dataclasses
    stripped = dataclasses.replace(fit, design=None)
    with pytest.raises(ValueError, match="does not carry `design`"):
        pace.fitted_means(stripped, counts)
