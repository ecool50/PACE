"""The figures.

These are not pixel comparisons against R. What they check is that each panel
draws the numbers the ported table holds, so a plot cannot quietly show
something the readout does not say.
"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

import pace

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FIX = Path(__file__).parent / "fixtures"


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def lines(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


@pytest.fixture(scope="module")
def readouts():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    fit = pace.pace_fit(counts, vec("raw_x.txt"), vec("raw_y.txt"),
                        lines("raw_celltype.txt"), genes=lines("raw_genes.txt"),
                        n_threads=1)
    slopes = pace.pace_shrink(fit, null_correlation=False)
    stats = pace.fit_statistics(fit, counts)
    blocks = pace.pace_decompose(fit, counts, stats=stats)
    frame = pace.single_frame_decomposition(fit, counts, blocks)
    drivers = pace.pace_drivers(fit, slopes, blocks, stats)
    pairs = pace.pair_variance(fit, counts, frame)
    return counts, fit, slopes, blocks, frame, drivers, pairs


def test_the_decomposition_leaves_the_residual_out(readouts):
    """R's figure stops near 40%, not 100: the residual block is not drawn.

    It is the convention, not an oversight -- what the figure compares is how
    much of a cell's variance the model attributes, and a residual bar towering
    over the rest says nothing.
    """
    *_, frame, _, _ = readouts
    _, ax = plt.subplots()
    pace.plot_decomposition(frame, ax=ax, zoom=False)
    drawn = {text.get_text() for text in ax.get_legend().get_texts()}
    assert "Residual" not in drawn
    assert drawn == {"Cell type", "Spatial cell state", "Spillover"}
    total = sum(patch.get_height() for container in ax.containers for patch in container)
    assert total < 100 * len(frame.focal_levels)


def test_the_decomposition_is_ordered_by_spatial_share(readouts):
    """R orders focals by descending spatial block, not alphabetically."""
    *_, frame, _, _ = readouts
    _, ax = plt.subplots()
    pace.plot_decomposition(frame, ax=ax, zoom=False)
    drawn = [label.get_text() for label in ax.get_xticklabels()]
    spatial = frame.block_names.index("Spatial %")
    expected = [frame.focal_levels[i].replace("_", " ")
                for i in np.argsort(-frame.focal_table[:, spatial], kind="stable")]
    assert drawn == expected


def test_the_decomposition_bars_are_the_table(readouts):
    *_, frame, _, _ = readouts
    _, ax = plt.subplots()
    pace.plot_decomposition(frame, ax=ax, zoom=False)
    order = np.argsort(-frame.focal_table[:, frame.block_names.index("Spatial %")],
                       kind="stable")
    for name, container in zip(["Cell type %", "Spatial %", "Spillover %"], ax.containers):
        heights = np.array([patch.get_height() for patch in container])
        np.testing.assert_allclose(
            heights, frame.focal_table[order, frame.block_names.index(name)], rtol=1e-12)


def test_the_zoom_panel_leaves_out_the_two_big_blocks(readouts):
    """It carries no legend of its own, as R's does not: it shares the main one."""
    *_, frame, _, _ = readouts
    _, zoom = pace.plot_decomposition(frame)
    drawn = {container.get_label() for container in zoom.containers}
    assert "Cell type" not in drawn and "Residual" not in drawn
    assert "Spatial cell state" in drawn
    assert zoom.get_legend() is None
    plt.close("all")


def test_the_heatmap_is_oriented_as_rs_is(readouts):
    """Focal on x, neighbour on y, so a COLUMN is one focal type."""
    *_, fit, _, _, _, _, pairs = readouts
    _, ax = plt.subplots()
    pace.plot_pair_heatmap(pairs, cell_types=fit.cell_types, ax=ax)
    assert ax.get_xlabel() == "Focal"
    assert ax.get_ylabel() == "Neighbour"
    grid = ax.images[0].get_array()
    types = list(fit.cell_types)
    for focal, neighbour, value in zip(pairs["focal"], pairs["neighbour"], pairs["val"]):
        assert abs(grid[types.index(neighbour), types.index(focal)] - value) < 1e-12
    for i in range(len(types)):
        assert np.ma.is_masked(grid[i, i]) or np.isnan(grid[i, i])
    plt.close("all")


def test_the_driver_bars_are_ranked_scores(readouts):
    *_, drivers, _ = readouts
    table = max(drivers.pairs.values(), key=len)
    _, ax = plt.subplots()
    pace.plot_drivers(drivers, table.focal, table.neighbour, ax=ax, n=5)
    widths = [patch.get_width() for patch in ax.containers[0]]
    np.testing.assert_allclose(widths, table.driver_score[:5], rtol=1e-12)
    plt.close("all")


def test_the_driver_panel_can_show_what_each_gene_is_made_of(readouts):
    """R's companion panel: the check that a top gene is not just abundant."""
    *_, frame, drivers, _ = readouts
    bars, blocks = pace.plot_drivers(drivers, "Macrophage", "Tumour", single_frame=frame, n=5)
    assert "Per-gene decomposition" in blocks.get_title()
    genes = [label.get_text() for label in blocks.get_xticklabels()]
    assert genes == list(drivers["Macrophage", "Tumour"].gene[:5])
    plt.close("all")


def test_an_empty_pair_says_so_rather_than_drawing_nothing(readouts):
    *_, drivers, _ = readouts
    empty = [k for k, t in drivers.pairs.items() if not len(t)]
    if not empty:
        pytest.skip("every pair has calls on this crop")
    table = drivers.pairs[empty[0]]
    _, ax = plt.subplots()
    pace.plot_drivers(drivers, table.focal, table.neighbour, ax=ax)
    assert any("no calls" in text.get_text() for text in ax.texts)
    plt.close("all")


def test_proximity_uses_rs_literal_breaks(readouts):
    """R bins at 0,2,4,6,8,Inf, with the last label naming the observed max."""
    counts, fit, *_ = readouts
    panels = pace.plot_proximity(fit, counts, ["MRC1", "APOC1"], "Macrophage", "Tumour")
    ticks = [label.get_text() for label in panels[0].get_xticklabels()]
    assert ticks[:4] == ["[0,2]", "(2,4]", "(4,6]", "(6,8]"]
    assert ticks[4].startswith("(8,") and not ticks[4].endswith("inf]")
    assert "Number of Tumour neighbours" in panels[0].get_xlabel()
    plt.close("all")


def test_the_shrinkage_panel_counts_the_calls(readouts):
    _, _, slopes, *_ = readouts
    _, ax = plt.subplots()
    pace.plot_shrinkage(slopes, ax=ax)
    assert str(int(np.sum(slopes.lfsr < 0.05))) in ax.get_title(loc="left")
    plt.close("all")


def test_a_fit_without_coordinates_cannot_draw_proximity(readouts):
    import dataclasses
    counts, fit, *_ = readouts
    stripped = dataclasses.replace(fit, x=None)
    with pytest.raises(ValueError, match="coordinates"):
        pace.plot_proximity(stripped, counts, ["MRC1"], "Macrophage", "Tumour")
