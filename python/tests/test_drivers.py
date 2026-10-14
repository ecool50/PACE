"""`pace_drivers()`, against R's per-pair driver tables.

The ranking is the load-bearing part: the same genes in the same order, not
merely the same scores. The filter is lfsr < 0.05 alone, with no specificity
cut at this layer, since spec already enters the score.
"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

import pace

FIX = Path(__file__).parent / "fixtures"

COLUMNS = {"driver_score": "driver_score", "driver_score_spec4": "driver_score_spec4",
           "b_clean": "b_clean", "spec": "spec", "focal_mean": "focal_mean",
           "R2_S": "r2_s", "mu_bar": "mu_bar", "alpha": "alpha", "V_S": "v_s",
           "V_resid": "v_resid", "V_total": "v_total", "lfsr": "lfsr",
           "sd_shrunk": "sd_shrunk"}


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def lines(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


@pytest.fixture(scope="module")
def pipeline():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    fit = pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                        lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                        n_threads=1)
    slopes = pace.pace_shrink(fit, null_correlation=False)
    stats = pace.fit_statistics(fit, counts)
    decomposition = pace.pace_decompose(fit, counts, stats=stats)
    drivers = pace.pace_drivers(fit, slopes, decomposition, stats)
    return counts, fit, slopes, stats, decomposition, drivers


def flatten(drivers, keys):
    """Every pair's rows in R's key order, as R stacks them."""
    pair, gene, columns = [], [], {name: [] for name in COLUMNS}
    for key in keys:
        table = drivers.pairs.get(key)
        if table is None or not len(table):
            continue
        pair += [key] * len(table)
        gene += list(table.gene)
        for name, attribute in COLUMNS.items():
            columns[name] += list(getattr(table, attribute))
    return pair, gene, columns


def test_the_same_pairs_and_genes_in_the_same_rank_order(pipeline):
    *_, drivers = pipeline
    pair, gene, _ = flatten(drivers, lines("drv_keys.txt"))
    assert pair == lines("drv_pair.txt")
    assert gene == lines("drv_gene.txt")


def test_every_driver_column_matches_rs(pipeline):
    *_, drivers = pipeline
    _, _, columns = flatten(drivers, lines("drv_keys.txt"))
    for name in COLUMNS:
        np.testing.assert_allclose(
            np.asarray(columns[name], dtype=float), vec(f"drv_{name}.txt"),
            rtol=0, atol=1e-05, err_msg=f"{name} disagrees with R")


def test_the_pair_statuses_match_rs(pipeline):
    *_, drivers = pipeline
    keys = lines("drv_keys.txt")
    mine = [drivers.pairs[k].status if k in drivers.pairs else "dropped (n_eff)"
            for k in keys]
    assert mine == lines("drv_status.txt")


def test_the_score_is_the_locked_definition(pipeline):
    """driver_score = b^2 * spec^2 * max(focal_mean, 0), and nothing else."""
    *_, drivers = pipeline
    for table in drivers.pairs.values():
        if not len(table):
            continue
        expected = (table.b_clean ** 2) * (table.spec ** 2) * np.maximum(table.focal_mean, 0)
        np.testing.assert_allclose(table.driver_score, expected, rtol=1e-09, atol=0)


def test_the_ranking_is_by_descending_score(pipeline):
    *_, drivers = pipeline
    for table in drivers.pairs.values():
        if len(table) > 1:
            assert np.all(np.diff(table.driver_score) <= 1e-12)


def test_only_lfsr_filters_the_table(pipeline):
    """No specificity cut at the reporting layer, by design."""
    *_, drivers = pipeline
    specs = np.concatenate([t.spec for t in drivers.pairs.values() if len(t)])
    assert np.all(np.concatenate([t.lfsr for t in drivers.pairs.values() if len(t)]) < 0.05)
    assert specs.min() < 0.5, "a spec filter would have removed these rows"


def test_a_pair_is_never_its_own_neighbour(pipeline):
    *_, drivers = pipeline
    for table in drivers.pairs.values():
        assert table.focal != table.neighbour


def test_a_requested_subset_of_pairs_is_honoured(pipeline):
    counts, fit, slopes, stats, decomposition, _ = pipeline
    one = pace.pace_drivers(fit, slopes, decomposition, stats,
                            pairs=[("Macrophage", "Tumour")])
    assert list(one.pairs) == ["Macrophage_Tumour"]
    assert one["Macrophage", "Tumour"] is one.pairs["Macrophage_Tumour"]


def test_the_expected_false_signs_are_the_summed_lfsr(pipeline):
    *_, drivers = pipeline
    keys = lines("drv_keys.txt")
    reference = vec("drv_efs.txt")
    for key, expected in zip(keys, reference):
        table = drivers.pairs.get(key)
        got = table.expected_false_sign if table is not None else 0.0
        assert abs(got - expected) < 1e-06


def test_it_can_become_a_dataframe(pipeline):
    pytest.importorskip("pandas")
    *_, drivers = pipeline
    frame = drivers.to_dataframe()
    assert list(frame.columns[:4]) == ["focal", "neighbour", "rank", "gene"]
    assert len(frame) == len(lines("drv_gene.txt"))
