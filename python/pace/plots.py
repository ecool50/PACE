"""Figures, drawn from the ported tables.

Each function here draws one of the R package's panels and takes the same
readout R's version takes, so the conventions that matter -- which blocks
appear, how a pair is normalised, what a bar is a fraction of -- come from the
tables rather than being decided again here. Where R's panel makes a choice
that is not implied by its data, this says so at the call.

matplotlib is imported lazily, so it is needed only if you draw.
"""
from __future__ import annotations

import numpy as np

# Tableau 10, which is what the R figures use.
TABLEAU = ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
           "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"]
BLOCK_COLOURS = {
    "Cell type %": "#4E79A7",
    "Spatial %": "#F28E2B",
    "Spatial state %": "#F28E2B",
    "Responder spatial %": "#E15759",
    "Responder spatial state %": "#E15759",
    "Spillover %": "#76B7B2",
    "Residual %": "#BAB0AC",
}


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:                       # pragma: no cover
        raise ImportError(
            "drawing needs matplotlib, which is not a hard dependency of pace; "
            "install it with `pip install pace-spatial[plot]`") from error
    return plt


def _bare(axis):
    axis.spines[["top", "right"]].set_visible(False)
    return axis


def plot_decomposition(single_frame, ax=None, title="Per-focal decomposition",
                       zoom=True):
    """Stacked bars of each focal type's variance blocks.

    Takes a `single_frame_decomposition()` result, because that is the table R
    draws this from: the shares are observed, pooled over genes by sum of
    squares. Every block the table carries is drawn, residual included, which
    is what R's `plotDecomposition()` does, so the bars sum to 100%.

    `zoom` adds R's companion panel. It is on by default for a reason worth
    stating: cell type and residual take almost all of the height, so the
    spatial and spillover blocks -- the ones PACE exists to estimate -- are a
    sliver on the full bar and unreadable without it.
    """
    plt = _pyplot()
    if ax is None:
        if zoom:
            _, axes = plt.subplots(1, 2, figsize=(12, 4.2),
                                   gridspec_kw={"width_ratios": [1.35, 1]})
            ax, zoom_ax = axes
        else:
            _, ax = plt.subplots(figsize=(7, 4))
            zoom_ax = None
    else:
        zoom_ax = None
    focals = single_frame.focal_levels
    bottom = np.zeros(len(focals))
    for column, name in enumerate(single_frame.block_names):
        values = single_frame.focal_table[:, column]
        ax.bar(focals, values, bottom=bottom, label=name.replace(" %", ""),
               color=BLOCK_COLOURS.get(name, TABLEAU[column % len(TABLEAU)]),
               edgecolor="white", linewidth=0.5)
        bottom = bottom + values
    ax.set_ylabel("% of observed variance")
    ax.set_title(title, loc="left")
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    _bare(ax)
    if zoom_ax is None:
        return ax

    # the small blocks on their own scale: everything but cell type and residual
    small = [(i, name) for i, name in enumerate(single_frame.block_names)
             if name not in ("Cell type %", "Residual %")]
    width = 0.8 / max(len(small), 1)
    positions = np.arange(len(focals))
    for slot, (column, name) in enumerate(small):
        zoom_ax.bar(positions + slot * width - 0.4 + width / 2,
                    single_frame.focal_table[:, column], width,
                    label=name.replace(" %", ""),
                    color=BLOCK_COLOURS.get(name, TABLEAU[column % len(TABLEAU)]))
    zoom_ax.set_xticks(positions, focals, rotation=45, ha="right")
    zoom_ax.set_ylabel("% of observed variance")
    zoom_ax.set_title("the spatial and spillover blocks, rescaled", loc="left")
    zoom_ax.legend(frameon=False, fontsize=8)
    _bare(zoom_ax)
    return ax, zoom_ax


def plot_pair_heatmap(pairs, cell_types=None, ax=None, cmap="inferno",
                      title="Spatial % by neighbour"):
    """The focal-by-neighbour heatmap of `pair_variance()`.

    Rows are focal types, columns neighbours, and a row sums to that focal's
    observed spatial share. The own-type diagonal is absent from the table by
    definition and is left blank rather than drawn as zero, since zero would
    read as "measured and nil".
    """
    plt = _pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 5))
    focals = list(dict.fromkeys(pairs["focal"])) if cell_types is None else list(cell_types)
    neighbours = list(dict.fromkeys(pairs["neighbour"])) if cell_types is None else list(cell_types)
    grid = np.full((len(focals), len(neighbours)), np.nan)
    for focal, neighbour, value in zip(pairs["focal"], pairs["neighbour"], pairs["val"]):
        grid[focals.index(focal), neighbours.index(neighbour)] = value

    image = ax.imshow(grid, cmap=cmap, aspect="auto")
    ax.set_xticks(range(len(neighbours)), neighbours, rotation=45, ha="right")
    ax.set_yticks(range(len(focals)), focals)
    ax.set_xlabel("neighbour")
    ax.set_ylabel("focal")
    ax.set_title(title, loc="left")
    bar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    bar.set_label("% of focal's variance")
    return ax


def plot_drivers(drivers, focal, neighbour, n=12, ax=None):
    """The top genes mediating one pair, by driver score.

    The bars are the score, which is direction-blind by construction, so the
    sign of the slope is shown by colour rather than by bar direction: a gene
    suppressed near the neighbour would otherwise look like a weak driver.
    """
    plt = _pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))
    table = drivers[focal, neighbour]
    if not len(table):
        ax.text(0.5, 0.5, f"no calls for {focal} <- {neighbour}", ha="center",
                va="center", transform=ax.transAxes)
        return _bare(ax)
    take = slice(0, min(n, len(table)))
    genes = table.gene[take][::-1]
    scores = table.driver_score[take][::-1]
    up = table.b_clean[take][::-1] > 0
    ax.barh(range(genes.size), scores,
            color=np.where(up, TABLEAU[1], TABLEAU[0]))
    ax.set_yticks(range(genes.size), genes)
    ax.set_xlabel("driver score  ($b^2 \\cdot spec^2 \\cdot focal\\ mean$)")
    ax.set_title(f"{focal} <- {neighbour}  ({table.status})", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=TABLEAU[1]),
               plt.Rectangle((0, 0), 1, 1, color=TABLEAU[0])]
    ax.legend(handles, ["up near neighbour", "down near neighbour"],
              frameon=False, fontsize=8, loc="lower right")
    return _bare(ax)


def plot_proximity(fit, counts, genes, focal, neighbour, ax=None, radius=30.0,
                   max_neighbours=6):
    """Raw counts of a gene against how many of a neighbour type a cell has.

    Binned by the NUMBER of neighbouring cells of that type within `radius`,
    which is what R's panel shows and what a reader can act on. The fitted
    kernel would be the other candidate, but it is centred within image, so its
    axis runs negative and means nothing on its own.

    The counts are zero-inflated, so a per-bin median is almost always zero and
    says nothing; the mean is drawn instead, as R does.
    """
    from scipy.spatial import cKDTree

    plt = _pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    if isinstance(genes, str):
        genes = [genes]
    if fit.x is None:
        raise ValueError("this fit does not carry coordinates; re-fit to draw this")

    types = list(fit.cell_types)
    focal_cells = np.flatnonzero(fit.type_code == types.index(focal))
    neighbour_cells = np.flatnonzero(fit.type_code == types.index(neighbour))
    tree = cKDTree(np.column_stack([fit.x[neighbour_cells], fit.y[neighbour_cells]]))
    found = tree.query_ball_point(
        np.column_stack([fit.x[focal_cells], fit.y[focal_cells]]), r=radius)
    n_near = np.array([len(hit) for hit in found])
    binned = np.clip(n_near, 0, max_neighbours)

    gene_row = {gene: i for i, gene in enumerate(fit.genes)}
    csr = counts.tocsr()
    present = [b for b in range(max_neighbours + 1) if np.any(binned == b)]
    for index, gene in enumerate(genes):
        column = np.asarray(csr[:, gene_row[gene]].todense()).ravel()[focal_cells]
        means = [column[binned == b].mean() for b in present]
        ax.plot(present, means, marker="o", color=TABLEAU[index % len(TABLEAU)], label=gene)
    ax.set_xticks(present,
                  [str(b) if b < max_neighbours else f"{max_neighbours}+" for b in present])
    ax.set_xlabel(f"{neighbour} cells within {radius:g} um")
    ax.set_ylabel("mean raw count")
    ax.set_title(f"{focal} cells, by {neighbour} proximity", loc="left")
    ax.legend(frameon=False, fontsize=8)
    return _bare(ax)


def plot_shrinkage(slopes, ax=None):
    """What the shrinkage did: the raw slope against the shrunken one."""
    plt = _pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(4.5, 4.5))
    raw, shrunk = slopes.estimate, slopes.estimate_shrunk
    # A real cohort has millions of these; drawing them all costs minutes and
    # shows nothing more than a sample does.
    if raw.size > 25_000:
        pick = np.random.default_rng(0).choice(raw.size, 25_000, replace=False)
        raw, shrunk = raw[pick], shrunk[pick]
    ax.scatter(raw, shrunk, s=5, alpha=0.3, color=TABLEAU[0], edgecolors="none")
    limit = float(np.nanpercentile(np.abs(raw), 99.5))
    ax.plot([-limit, limit], [-limit, limit], ls="--", lw=0.8, color="grey")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("raw slope")
    ax.set_ylabel("shrunken slope")
    called = int(np.sum(slopes.lfsr < 0.05))
    ax.set_title(f"{called} calls at lfsr < 0.05", loc="left")
    return _bare(ax)
