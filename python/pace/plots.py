"""Figures, drawn from the ported tables.

Each function draws one of the R package's panels from the readout R's version
takes, so the conventions come from the data rather than being decided again
here, and the palettes are R's own hex values rather than approximations.

matplotlib is imported lazily, so it is needed only if you draw.
"""
from __future__ import annotations

import numpy as np

# R/plots.R's palettes, verbatim.
BLOCK_COLOURS = {
    "Cell type": "#003049",
    "Spatial cell state": "#780000",
    "Responder spatial state": "#c1121f",
    "Spillover": "#E76F51",
}
GENE_BLOCK_COLOURS = {
    "Cell type": "#1F3B57",
    "Spatial cell state": "#8B1A1A",
    "Responder spatial state": "#E07B39",
    "Spillover": "#F6C9C9",
}
DRIVER_BAR = "#F5A623"
BOX_COLOUR = "#4F8B5E"
HEAT_COLOURS = ["#000000", "#4D0000", "#8B0000", "#E41A1C", "#FF7F00",
                "#FFD92F", "#FFFFB2"]


def _pretty(name):
    """R's figures show "B Cell" where the data says "B_Cell"."""
    return str(name).replace("_", " ")


def _block_label(column_name):
    """"Spatial %" in the table is "Spatial cell state" on the figure."""
    stripped = column_name.replace(" %", "")
    return {"Spatial": "Spatial cell state",
            "Responder spatial": "Responder spatial state"}.get(stripped, stripped)


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
    """Stacked bars of each focal type's variance blocks, as R draws them.

    Two things here are R's conventions rather than the data's, and both are
    deliberate. The RESIDUAL block is left out, so the bars stop around 40%
    instead of reaching 100: what is being compared is how much of a cell's
    variance the model attributes, and a residual bar towering over the rest
    tells the reader nothing. And the focal types are ordered by descending
    spatial share rather than alphabetically, so the types PACE has something
    to say about come first.

    `zoom` adds R's companion panel on the two small blocks, which are a
    sliver on the main bar and unreadable without it.
    """
    plt = _pyplot()
    drawn = [(column, _block_label(name))
             for column, name in enumerate(single_frame.block_names)
             if name not in ("Cell type %", "Residual %")]
    spatial = next(column for column, label in drawn if label == "Spatial cell state")
    order = np.argsort(-single_frame.focal_table[:, spatial], kind="stable")
    focals = [_pretty(single_frame.focal_levels[i]) for i in order]
    table = single_frame.focal_table[order, :]

    if ax is None:
        if zoom:
            figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            ax, zoom_ax = axes
        else:
            figure, ax = plt.subplots(figsize=(6, 4.5))
            zoom_ax = None
    else:
        zoom_ax = None

    cell_type = single_frame.block_names.index("Cell type %")
    bottom = table[:, cell_type].copy()
    ax.bar(focals, bottom, label="Cell type", color=BLOCK_COLOURS["Cell type"], width=0.7)
    for column, label in drawn:
        ax.bar(focals, table[:, column], bottom=bottom, label=label,
               color=BLOCK_COLOURS[label], width=0.7)
        bottom = bottom + table[:, column]
    ax.set_ylabel("Proportion of total variance")
    ax.set_xlabel("Focal cell type")
    ax.set_title(title, loc="center")
    ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0f}%")
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    # Below the axis label, not on top of it: the tick labels are rotated, so
    # the axis needs more room underneath than the default anchor leaves.
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.52), title="Block")
    _bare(ax)
    if zoom_ax is None:
        return ax

    width = 0.8 / len(drawn)
    positions = np.arange(len(focals))
    for slot, (column, label) in enumerate(drawn):
        offsets = positions + slot * width - 0.4 + width / 2
        values = table[:, column]
        zoom_ax.bar(offsets, values, width, label=label, color=BLOCK_COLOURS[label])
        for position, value in zip(offsets, values):
            zoom_ax.annotate(f"{value:.2f}%", (position, value), ha="center",
                             va="bottom", fontsize=6, rotation=90)
    zoom_ax.set_xticks(positions, focals, rotation=45, ha="right")
    zoom_ax.set_ylabel("Proportion of total variance")
    zoom_ax.set_xlabel("Focal cell type")
    zoom_ax.set_title(" & ".join(label for _, label in drawn), loc="center")
    zoom_ax.yaxis.set_major_formatter(lambda value, _: f"{value:.1f}%")
    zoom_ax.set_ylim(0, zoom_ax.get_ylim()[1] * 1.25)
    _bare(zoom_ax)
    return ax, zoom_ax


def plot_pair_heatmap(pairs, cell_types=None, ax=None,
                      title="Pairwise spatial % (variance, normalised Pratt)"):
    """The pair heatmap, oriented as R draws it: FOCAL on x, NEIGHBOUR on y.

    Each cell is the percentage of the focal's total variance contributed by
    spatial interaction with that neighbour, so a COLUMN sums to the focal's
    observed spatial share. A pair the table does not carry is drawn as zero
    rather than left blank, which is R's choice; the own-type diagonal is left
    blank, since it is excluded by definition rather than measured and nil.
    """
    from matplotlib.colors import LinearSegmentedColormap

    plt = _pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5.5))
    types = list(cell_types) if cell_types is not None else sorted(
        set(pairs["focal"]) | set(pairs["neighbour"]))
    grid = np.zeros((len(types), len(types)))
    for focal, neighbour, value in zip(pairs["focal"], pairs["neighbour"], pairs["val"]):
        grid[types.index(neighbour), types.index(focal)] = value
    np.fill_diagonal(grid, np.nan)

    cmap = LinearSegmentedColormap.from_list("pace_heat", HEAT_COLOURS)
    cmap.set_bad("white")
    image = ax.imshow(grid, cmap=cmap, aspect="auto")
    for row in range(len(types)):
        for column in range(len(types)):
            if row == column:
                continue
            value = grid[row, column]
            ax.text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if value < np.nanmax(grid) * 0.6 else "black")
    ax.set_xticks(range(len(types)), types, rotation=45, ha="right")
    ax.set_yticks(range(len(types)), types)
    ax.set_xlabel("Focal")
    ax.set_ylabel("Neighbour")
    ax.set_title(title, loc="center", fontsize=10)
    bar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    bar.set_label("spatial %\nof total")
    return ax


def plot_drivers(drivers, focal, neighbour, single_frame=None, n=5, ax=None):
    """The top genes mediating one pair, as R draws them.

    Two panels. The left ranks the genes by driver score; the right shows what
    each one's variance is made of, from the observed single frame, which is
    the check on the ranking: a top gene whose variance is nearly all cell-type
    identity is being driven by how much the focal expresses it rather than by
    the neighbour. Pass `single_frame` to get that panel.

    The bars are one colour, as R's are. The score is direction-blind by
    construction, so the sign lives in the table rather than the figure.
    """
    plt = _pyplot()
    table = drivers[focal, neighbour]
    if ax is None:
        if single_frame is not None:
            _, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            ax, block_ax = axes
        else:
            _, ax = plt.subplots(figsize=(5.5, 4.5))
            block_ax = None
    else:
        block_ax = None

    if not len(table):
        ax.text(0.5, 0.5, f"no calls for {focal} <- {neighbour}", ha="center",
                va="center", transform=ax.transAxes)
        return _bare(ax)

    take = slice(0, min(n, len(table)))
    genes = list(table.gene[take])
    scores = table.driver_score[take]
    ax.barh(np.arange(len(genes))[::-1], scores, color=DRIVER_BAR, height=0.7)
    for position, value in zip(np.arange(len(genes))[::-1], scores):
        ax.annotate(f"{value:.3f}", (value, position), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=7)
    ax.set_yticks(np.arange(len(genes))[::-1], genes)
    ax.set_xlabel("Driver score (spatial cell state)")
    ax.set_title(f"{_pretty(focal)} <- {_pretty(neighbour)}", loc="left")
    ax.set_xlim(0, scores.max() * 1.2)
    _bare(ax)
    if block_ax is None:
        return ax

    # what each top gene's variance is made of, for this focal type
    rows = single_frame.focal == focal
    position = {gene: i for i, gene in enumerate(single_frame.gene[rows])}
    shares = {"Cell type": single_frame.pct_celltype[rows],
              "Spatial cell state": single_frame.pct_spatial[rows],
              "Spillover": single_frame.pct_spill[rows]}
    if single_frame.pct_responder is not None:
        shares["Responder spatial state"] = single_frame.pct_responder[rows]
    order = ["Cell type", "Spatial cell state"] + (
        ["Responder spatial state"] if single_frame.pct_responder is not None else []) + ["Spillover"]
    index = [position[gene] for gene in genes if gene in position]
    drawn_genes = [gene for gene in genes if gene in position]
    bottom = np.zeros(len(drawn_genes))
    for label in order:
        values = shares[label][index]
        block_ax.bar(drawn_genes, values, bottom=bottom, label=label,
                     color=GENE_BLOCK_COLOURS[label], width=0.7)
        bottom = bottom + values
    block_ax.set_ylabel("Proportion of total variance")
    block_ax.set_xlabel("Gene")
    block_ax.set_title("Per-gene decomposition (single-frame)", loc="center")
    block_ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0f}%")
    block_ax.tick_params(axis="x", rotation=45)
    for label in block_ax.get_xticklabels():
        label.set_horizontalalignment("right")
    block_ax.legend(frameon=False, fontsize=8, title="Block",
                    bbox_to_anchor=(1.01, 1), loc="upper left")
    _bare(block_ax)
    return ax, block_ax


def plot_proximity(fit, counts, genes, focal, neighbour, ax=None, radius=30.0,
                   breaks=(0, 2, 4, 6, 8, np.inf)):
    """A gene's raw counts against how many of a neighbour type a cell has.

    One panel per gene, drawn as R draws it: boxplots over bins of the NEIGHBOUR
    COUNT within `radius`, with the bin mean as a filled point. The mean is
    what carries the trend here -- these counts are zero-inflated, so the box
    median sits at zero in most bins and says nothing.

    The fitted kernel is the other thing one could bin on, but it is centred
    within image, so its axis runs negative and means nothing on its own.
    """
    from scipy.spatial import cKDTree

    plt = _pyplot()
    if isinstance(genes, str):
        genes = [genes]
    if fit.x is None:
        raise ValueError("this fit does not carry coordinates; re-fit to draw this")
    if ax is None:
        _, axes = plt.subplots(1, len(genes), figsize=(5.5 * len(genes), 4.2), squeeze=False)
        axes = list(axes[0])
    else:
        axes = [ax] * len(genes) if len(genes) == 1 else list(np.atleast_1d(ax))

    types = list(fit.cell_types)
    focal_cells = np.flatnonzero(fit.type_code == types.index(focal))
    neighbour_cells = np.flatnonzero(fit.type_code == types.index(neighbour))
    tree = cKDTree(np.column_stack([fit.x[neighbour_cells], fit.y[neighbour_cells]]))
    n_near = np.array([len(hit) for hit in tree.query_ball_point(
        np.column_stack([fit.x[focal_cells], fit.y[focal_cells]]), r=radius)])

    # R's default breaks are literal, not derived from the data: c(0,2,4,6,8,Inf),
    # closed on the left of the first bin. Its last label names the observed
    # maximum rather than printing Inf, so "(8,19]" rather than "(8,Inf]".
    edges = np.asarray(breaks, dtype=float)
    which = np.clip(np.digitize(n_near, edges[1:-1], right=True), 0, len(edges) - 2)
    labels = [f"[{edges[0]:g},{edges[1]:g}]"] + [
        f"({edges[b]:g},{edges[b + 1]:g}]" for b in range(1, len(edges) - 1)]
    labels[-1] = f"({edges[-2]:g},{n_near.max():d}]"

    gene_row = {gene: i for i, gene in enumerate(fit.genes)}
    csr = counts.tocsr()
    for axis, gene in zip(axes, genes):
        column = np.asarray(csr[:, gene_row[gene]].todense()).ravel()[focal_cells]
        groups = [column[which == b] for b in range(len(labels))]
        axis.boxplot(groups, tick_labels=labels, showfliers=True, widths=0.6,
                     boxprops=dict(color=BOX_COLOUR), whiskerprops=dict(color=BOX_COLOUR),
                     capprops=dict(color=BOX_COLOUR), medianprops=dict(color=BOX_COLOUR),
                     flierprops=dict(marker=".", markersize=3,
                                     markerfacecolor=BOX_COLOUR,
                                     markeredgecolor="none", alpha=0.35))
        axis.plot(np.arange(1, len(groups) + 1),
                  [group.mean() if group.size else np.nan for group in groups],
                  "o", color=BOX_COLOUR, markersize=8, linestyle="none")
        axis.set_title(gene, loc="center")
        axis.set_xlabel(f"Number of {_pretty(neighbour)} neighbours")
        axis.set_ylabel(f"{_pretty(focal)} RAW counts")
        _bare(axis)
    return axes[0] if len(axes) == 1 else axes


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
    ax.scatter(raw, shrunk, s=5, alpha=0.3, color=BLOCK_COLOURS["Cell type"],
               edgecolors="none")
    limit = float(np.nanpercentile(np.abs(raw), 99.5))
    ax.plot([-limit, limit], [-limit, limit], ls="--", lw=0.8, color="grey")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("raw slope")
    ax.set_ylabel("shrunken slope")
    called = int(np.sum(slopes.lfsr < 0.05))
    ax.set_title(f"{called} calls at lfsr < 0.05", loc="left")
    return _bare(ax)
