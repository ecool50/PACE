"""`pace_decompose()`, against R's five-block table.

Both languages call the same C++, so a disagreement here is a marshalling
error: a wrong index vector, a transposed matrix, a column-major array handed
over as row-major. That is not hypothetical -- the index matrices really were
being read transposed, and because they are square nothing raised, the table
was just built from the wrong rows of U.

Tolerances are the fit's trigamma difference carried through. The percentages
are on a 0-100 scale, so theirs is correspondingly larger.
"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

import pace

FIX = Path(__file__).parent / "fixtures"

BLOCK_COLUMNS = [
    ("v_state_baseline", "V_state_baseline"), ("v_state_responder", "V_state_responder"),
    ("v_spill", "V_spill"), ("v_disp", "V_disp"),
    ("celltype_offset_sq", "celltype_offset_sq"), ("focal_mean", "focal_mean"),
    ("max_other_mean", "max_other_mean"), ("spec", "spec"),
    ("focal_other_ratio", "focal_other_ratio"),
]
PERCENT_KEYS = ["Cell_type", "Spatial_state", "Responder_spatial_state",
                "Spillover", "Residual"]


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


def build(counts, condition):
    extra = {} if not condition else dict(
        condition=lines("cond2_label.txt"), condition_name="Arm",
        condition_levels=["PD", "PR"])
    return pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                         lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                         n_threads=1, **extra)


@pytest.fixture(scope="module", params=[("dec", False), ("decc", True)],
                ids=["no-condition", "condition"])
def case(request, counts):
    tag, condition = request.param
    fit = build(counts, condition)
    return tag, pace.pace_decompose(fit, counts)


def test_the_rows_are_the_same_as_rs(case):
    tag, decomposition = case
    assert list(decomposition.gene) == lines(f"{tag}_gene.txt")
    assert list(decomposition.focal) == lines(f"{tag}_focal.txt")


def test_every_variance_block_matches_rs(case):
    tag, decomposition = case
    for field, fixture in BLOCK_COLUMNS:
        np.testing.assert_allclose(
            getattr(decomposition, field), vec(f"{tag}_{fixture}.txt"),
            rtol=0, atol=1e-05, err_msg=f"{field} disagrees with R")


def test_the_percentages_match_rs(case):
    tag, decomposition = case
    reference = np.column_stack([vec(f"{tag}_pct_{key}.txt") for key in PERCENT_KEYS])
    np.testing.assert_allclose(decomposition.pct, reference, rtol=0, atol=1e-04)


def test_the_contamination_flag_matches_rs(case):
    tag, decomposition = case
    assert np.array_equal(decomposition.is_contaminated,
                          vec(f"{tag}_is_contaminated.txt", int) == 1)


@pytest.mark.parametrize("weighting", ["mean", "specw", "pooled"])
def test_the_per_focal_aggregates_match_rs(case, weighting):
    tag, decomposition = case
    table = {"mean": decomposition.mean_table, "specw": decomposition.specw_table,
             "pooled": decomposition.pooled_table}[weighting]
    reference = vec(f"{tag}_agg_{weighting}.txt").reshape(table.shape, order="F")
    np.testing.assert_allclose(table, reference, rtol=0, atol=1e-05)


def test_the_index_matrices_are_not_read_transposed(counts):
    """The bug this file exists for.

    `slope_rows` is n_kernel x n_focals column-major. Read transposed it picks
    another focal's slopes, and since it is square nothing complains. Checked
    by recomputing one block independently of the core.
    """
    fit = build(counts, False)
    types = fit.cell_types
    n_types = len(types)
    covariance = np.asarray(pace._core.group_covariances(
        np.asfortranarray(fit.kernel), fit.type_code, n_types, 1))
    slab = n_types * n_types
    by_hand = np.zeros((n_types, len(fit.genes)))
    for focal in range(n_types):
        sigma = covariance[focal * slab:(focal + 1) * slab].reshape(
            n_types, n_types, order="F")
        rows = [(k + 1) * n_types + focal for k in range(n_types)]
        u, se = fit.U[rows, :], fit.se_U[rows, :]
        by_hand[focal, :] = (np.einsum("kg,kl,lg->g", u, sigma, u)
                             + np.diag(sigma) @ (se ** 2))
    decomposition = pace.pace_decompose(fit, counts)
    np.testing.assert_allclose(decomposition.v_state_baseline, by_hand.ravel(),
                               rtol=0, atol=1e-08)


def test_a_condition_fit_carries_a_responder_block(counts):
    with_condition = pace.pace_decompose(build(counts, True), counts)
    without = pace.pace_decompose(build(counts, False), counts)
    assert np.any(with_condition.v_state_responder > 0)
    assert np.all(without.v_state_responder == 0)


def test_it_can_become_dataframes(case):
    pytest.importorskip("pandas")
    _, decomposition = case
    frame = decomposition.to_dataframe()
    assert "Spatial state %" in frame.columns
    assert len(frame) == len(decomposition)
    aggregate = decomposition.aggregates("specw")
    assert list(aggregate["focal"]) == decomposition.focal_levels


def test_a_fit_without_the_kernel_says_so(counts):
    import dataclasses
    stripped = dataclasses.replace(build(counts, False), kernel=None)
    with pytest.raises(ValueError, match="does not carry `kernel`"):
        pace.pace_decompose(stripped, counts)
