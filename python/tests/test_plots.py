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


def test_the_decomposition_bars_are_the_table(readouts):
    *_, frame, _, _ = readouts
    _, ax = plt.subplots()
    pace.plot_decomposition(frame, ax=ax, zoom=False)
    # one container per block, each bar the table's value
    for column, container in enumerate(ax.containers):
        heights = np.array([patch.get_height() for patch in container])
        np.testing.assert_allclose(heights, frame.focal_table[:, column], rtol=1e-12)
    plt.close("all")


def test_the_zoom_panel_leaves_out_the_two_big_blocks(readouts):
    """Its whole purpose: the spatial blocks are unreadable on the full bar."""
    *_, frame, _, _ = readouts
    full, zoom = pace.plot_decomposition(frame)
    drawn = {text.get_text() for text in zoom.get_legend().get_texts()}
    assert "Cell type" not in drawn and "Residual" not in drawn
    assert "Spatial" in drawn
    plt.close("all")


def test_the_heatmap_holds_the_pair_values(readouts):
    *_, fit, _, _, _, _, pairs = readouts
    _, ax = plt.subplots()
    pace.plot_pair_heatmap(pairs, cell_types=fit.cell_types, ax=ax)
    grid = ax.images[0].get_array()
    types = list(fit.cell_types)
    for focal, neighbour, value in zip(pairs["focal"], pairs["neighbour"], pairs["val"]):
        assert abs(grid[types.index(focal), types.index(neighbour)] - value) < 1e-12
    # the own-type diagonal is absent, not zero
    for i in range(len(types)):
        assert np.ma.is_masked(grid[i, i]) or np.isnan(grid[i, i])
    plt.close("all")


def test_the_driver_bars_are_ranked_scores(readouts):
    *_, drivers, _ = readouts
    table = max(drivers.pairs.values(), key=len)
    _, ax = plt.subplots()
    pace.plot_drivers(drivers, table.focal, table.neighbour, ax=ax, n=5)
    widths = [patch.get_width() for patch in ax.containers[0]]
    np.testing.assert_allclose(widths, table.driver_score[:5][::-1], rtol=1e-12)
    plt.close("all")


def test_the_driver_colour_shows_the_direction(readouts):
    """The score is direction-blind, so colour is the only sign there is."""
    *_, drivers, _ = readouts
    table = max(drivers.pairs.values(), key=len)
    _, ax = plt.subplots()
    pace.plot_drivers(drivers, table.focal, table.neighbour, ax=ax, n=5)
    colours = [patch.get_facecolor() for patch in ax.containers[0]]
    up = table.b_clean[:5][::-1] > 0
    distinct = {tuple(c) for c in colours}
    assert len(distinct) == (2 if (up.any() and (~up).any()) else 1)
    plt.close("all")


def test_an_empty_pair_says_so_rather_than_drawing_nothing(readouts):
    *_, fit, _, _, _, drivers, _ = readouts
    empty = [k for k, t in drivers.pairs.items() if not len(t)]
    if not empty:
        pytest.skip("every pair has calls on this crop")
    focal, neighbour = drivers.pairs[empty[0]].focal, drivers.pairs[empty[0]].neighbour
    _, ax = plt.subplots()
    pace.plot_drivers(drivers, focal, neighbour, ax=ax)
    assert any("no calls" in text.get_text() for text in ax.texts)
    plt.close("all")


def test_proximity_bins_by_neighbour_count_not_the_centred_kernel(readouts):
    """The kernel is centred within image, so its axis would run negative."""
    counts, fit, *_ = readouts
    _, ax = plt.subplots()
    pace.plot_proximity(fit, counts, ["MRC1", "APOC1"], "Macrophage", "Tumour", ax=ax)
    ticks = [label.get_text() for label in ax.get_xticklabels()]
    assert ticks[0] == "0"
    assert all(not tick.startswith("-") for tick in ticks)
    assert "um" in ax.get_xlabel()
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
