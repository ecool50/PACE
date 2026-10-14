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


# --- the observed single frame -------------------------------------------
# This is the table R's varianceDecomposition() actually returns, and it is
# not the five-block one: the magnitudes are observed (between-type sums of
# squares of log1p(CP10k)) and only the within-type split comes from the fit.
# So its observed columns should agree with R far more tightly than the
# fitted ones, which carry the trigamma difference.

@pytest.fixture(scope="module", params=[("sf", False), ("sfc", True)],
                ids=["no-condition", "condition"])
def single_frame(request, counts):
    tag, condition = request.param
    fit = build(counts, condition)
    return tag, pace.single_frame_decomposition(fit, counts)


def test_the_single_frame_rows_are_the_same_as_rs(single_frame):
    tag, frame = single_frame
    assert list(frame.gene) == lines(f"{tag}_gene.txt")
    assert list(frame.focal) == lines(f"{tag}_focal.txt")
    assert np.array_equal(frame.n_focal, vec(f"{tag}_n_focal.txt", int))


def test_the_observed_sums_of_squares_match_rs_almost_exactly(single_frame):
    """These come from the counts, not the fit, so nothing should drift."""
    tag, frame = single_frame
    for field, fixture in [("ss_lineage", "SS_lineage"), ("ss_within", "SS_within"),
                           ("denom", "denom"), ("pct_celltype", "Cell_type")]:
        np.testing.assert_allclose(getattr(frame, field), vec(f"{tag}_{fixture}.txt"),
                                   rtol=1e-09, atol=1e-09, err_msg=field)


def test_the_within_type_split_matches_rs(single_frame):
    tag, frame = single_frame
    for field, fixture in [("pct_spatial", "Spatial"), ("pct_spill", "Spillover"),
                           ("pct_residual", "Residual")]:
        np.testing.assert_allclose(getattr(frame, field), vec(f"{tag}_{fixture}.txt"),
                                   rtol=0, atol=1e-04, err_msg=field)


def test_only_a_condition_fit_has_a_responder_share(single_frame):
    tag, frame = single_frame
    if tag == "sfc":
        np.testing.assert_allclose(frame.pct_responder,
                                   vec("sfc_Responder_spatial.txt"), rtol=0, atol=1e-04)
        assert "Responder spatial %" in frame.block_names
    else:
        assert frame.pct_responder is None
        assert "Responder spatial %" not in frame.block_names


def test_the_shares_sum_to_one_hundred(single_frame):
    _, frame = single_frame
    total = frame.pct_celltype + frame.pct_spatial + frame.pct_spill + frame.pct_residual
    if frame.pct_responder is not None:
        total = total + frame.pct_responder
    finite = np.isfinite(total)
    np.testing.assert_allclose(total[finite], 100.0, rtol=0, atol=1e-06)


def test_the_per_focal_table_is_sum_of_squares_weighted(single_frame):
    """block = 100 * sum(pct/100 * denom) / sum(denom), per focal."""
    _, frame = single_frame
    for row, focal in enumerate(frame.focal_levels):
        rows = frame.focal == focal
        weight = frame.denom[rows]
        for column, name in enumerate(frame.block_names):
            share = {"Cell type %": frame.pct_celltype, "Spatial %": frame.pct_spatial,
                     "Responder spatial %": frame.pct_responder,
                     "Spillover %": frame.pct_spill,
                     "Residual %": frame.pct_residual}[name][rows]
            expected = 100.0 * np.sum(share / 100.0 * weight) / np.sum(weight)
            assert abs(frame.focal_table[row, column] - expected) < 1e-08


def test_it_can_become_dataframes_too(single_frame):
    pytest.importorskip("pandas")
    _, frame = single_frame
    assert "SS_lineage" in frame.to_dataframe().columns
    assert list(frame.by_focal()["focal"]) == frame.focal_levels


# --- pairwise spatial variance -------------------------------------------
# The definition is particular and easy to get wrong: the focal's spatial share
# comes from the OBSERVED single frame, and that one number is split across the
# focal's neighbours by a Pratt attribution renormalised over the OFF-DIAGONAL
# pairs alone. So a row sums to the focal's observed spatial share, and the
# own-type pair is excluded rather than shown as a large diagonal.

def _condition_fit(counts):
    return pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                         lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                         condition=lines("cond2_label.txt"), condition_name="Arm",
                         condition_levels=["PD", "PR"], n_threads=1)


@pytest.mark.parametrize("tag,condition,block", [
    ("pv", False, "spatial"), ("pvcs", True, "spatial"), ("pvcr", True, "responder")])
def test_the_pair_table_matches_rs(counts, tag, condition, block):
    fit = _condition_fit(counts) if condition else build(counts, False)
    table = pace.pair_variance(fit, counts, block=block)
    position = {key: i for i, key in enumerate(zip(table["focal"], table["neighbour"]))}
    keys = list(zip(lines(f"{tag}_focal.txt"), lines(f"{tag}_neighbour.txt")))
    assert len(keys) == table["val"].size
    order = np.array([position[key] for key in keys])
    np.testing.assert_allclose(table["val"][order], vec(f"{tag}_val.txt"),
                               rtol=0, atol=1e-07)


def test_each_focal_row_sums_to_its_observed_spatial_share(counts):
    """The definition, pinned: observed share split across neighbours."""
    fit = build(counts, False)
    frame = pace.single_frame_decomposition(fit, counts)
    table = pace.pair_variance(fit, counts, frame)
    column = frame.block_names.index("Spatial %")
    for row, focal in enumerate(frame.focal_levels):
        rows = table["focal"] == focal
        if not rows.any():
            continue
        assert abs(table["val"][rows].sum() - frame.focal_table[row, column]) < 1e-09


def test_the_own_type_pair_is_excluded(counts):
    fit = build(counts, False)
    table = pace.pair_variance(fit, counts)
    assert not np.any(table["focal"] == table["neighbour"])


def test_the_responder_block_needs_a_condition(counts):
    fit = build(counts, False)
    with pytest.raises(ValueError, match="no responder block"):
        pace.pair_variance(fit, counts, block="responder")


def test_an_unknown_block_is_refused(counts):
    fit = build(counts, False)
    with pytest.raises(ValueError, match='block must be'):
        pace.pair_variance(fit, counts, block="spillover")
