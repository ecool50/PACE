"""`pace_shrink()`, against R's `pace_shrink()` on the same fit.

The target is the called set, not the call count: two runs can agree on how
many effects pass lfsr < 0.05 while disagreeing about which. Both of R's
`null_correlation` settings are checked, because R defaults to TRUE and the
manuscript's counts were produced with FALSE.

R's table comes out in alphabetical gene order, because `mvpql_to_results_multi`
splits the long frame by gene and R's `split()` sorts. Python keeps the fit's
gene order, so the comparison aligns on (gene, focal, neighbour) rather than
assuming a shared row order.
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
        return [line.rstrip("\n") for line in handle if line.strip()]


@pytest.fixture(scope="module")
def fit():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    return pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                         lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                         n_threads=1)


def aligned(slopes, tag):
    """The index that puts Python's rows in R's row order."""
    position = {key: i for i, key in
                enumerate(zip(slopes.gene, slopes.focal, slopes.neighbour))}
    keys = list(zip(lines(f"{tag}_gene.txt"), lines(f"{tag}_focal.txt"),
                    lines(f"{tag}_neighbour.txt")))
    assert len(keys) == len(slopes), "R and Python disagree on how many rows there are"
    return np.array([position[key] for key in keys])


@pytest.mark.parametrize("tag,null_correlation", [("shrinkI", False), ("shrinkV", True)])
def test_the_called_set_is_the_same_as_rs(fit, tag, null_correlation):
    """The load-bearing assertion: the same gene-focal pairs, not just as many."""
    slopes = pace.pace_shrink(fit, null_correlation=null_correlation)
    order = aligned(slopes, tag)
    assert np.array_equal(slopes.lfsr[order] < 0.05, vec(f"{tag}_lfsr.txt") < 0.05)


@pytest.mark.parametrize("tag,null_correlation", [("shrinkI", False), ("shrinkV", True)])
def test_the_shrunken_estimates_match_rs(fit, tag, null_correlation):
    slopes = pace.pace_shrink(fit, null_correlation=null_correlation)
    order = aligned(slopes, tag)
    for field, fixture, tolerance in [("estimate", "estimate", 1e-06),
                                      ("std_error", "std_error", 1e-06),
                                      ("estimate_shrunk", "estimate_shrunk", 1e-05),
                                      ("sd_shrunk", "sd_shrunk", 1e-05),
                                      ("lfsr", "lfsr", 1e-04)]:
        np.testing.assert_allclose(getattr(slopes, field)[order],
                                   vec(f"{tag}_{fixture}.txt"), rtol=0, atol=tolerance,
                                   err_msg=f"{field} disagrees with R")


def test_every_focal_neighbour_gene_is_present(fit):
    slopes = pace.pace_shrink(fit, null_correlation=False)
    assert len(slopes) == 8 * 8 * 278
    assert set(slopes.focal) == set(fit.cell_types)
    assert set(slopes.neighbour) == set(fit.cell_types)


def test_the_term_names_are_the_plain_neighbour_without_a_condition(fit):
    slopes = pace.pace_shrink(fit, null_correlation=False)
    assert set(slopes.term) == set(fit.cell_types)


def test_the_reported_standard_error_is_floored_like_rs(fit):
    slopes = pace.pace_shrink(fit, null_correlation=False)
    assert slopes.std_error.min() >= 1e-4


def test_the_call_summary_counts_what_it_says(fit):
    slopes = pace.pace_shrink(fit, null_correlation=False)
    summary = slopes.calls()
    assert summary["n_calls"] == int(np.sum(slopes.lfsr < 0.05))
    assert 0 < summary["expected_false_sign"] < summary["n_calls"]


def test_it_can_become_a_dataframe_with_rs_column_names(fit):
    pytest.importorskip("pandas")
    frame = pace.pace_shrink(fit, null_correlation=False).to_dataframe()
    assert list(frame.columns) == ["gene", "focal", "neighbour", "term", "estimate",
                                   "std.error", "estimate_shrunk", "sd_shrunk", "lfsr"]
    assert len(frame) == 8 * 8 * 278


# --- a condition fit shrinks the responder layer, not the baseline -------

def test_a_condition_fit_shrinks_the_interaction_term():
    """R looks up `<resp_term>:<neighbour>`, so what comes back is responder-only."""
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    conditioned = pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                                lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                                condition=lines("cond2_label.txt"), condition_name="Arm",
                                n_threads=1)
    slopes = pace.pace_shrink(conditioned, null_correlation=False)
    assert all(term.startswith("ArmPR:") for term in set(slopes.term))
    assert set(slopes.term) == {f"ArmPR:{t}" for t in conditioned.cell_types}


def test_the_responder_term_defaults_to_the_last_level():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    conditioned = pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                                lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                                condition=lines("cond3_label.txt"), condition_name="Arm3",
                                n_threads=1)
    slopes = pace.pace_shrink(conditioned, null_correlation=False)
    assert all(term.startswith("Arm3PR:") for term in set(slopes.term))
