"""Sequencing the stages of a PACE fit.

Every numerical stage lives in the compiled core and is shared with the R
package. What this module does is order them, build the designs, and hand the
solver what it expects -- the job `pace_fit_streaming()` does in R.

`fit_anndata()` is a thin reader on top of `pace_fit()` for callers holding
an AnnData rather than arrays.

A condition cohort needs the part of R's model.matrix that PACE actually uses:
treatment contrasts, the expansion of `cond * (types)`, and the order the
resulting columns come out in. That is reproduced here for the two formulas
PACE builds rather than in general, and checked against R on a two-level and a
three-level condition, because only the three-level case can show that the
condition level varies fastest within each type.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.sparse import csc_matrix, issparse

from . import _core


@dataclass
class PaceFit:
    """A fitted model, with the pieces the readouts need."""

    B: np.ndarray            # fixed effects, p x genes
    U: np.ndarray            # random effects, q x genes
    se_B: np.ndarray
    se_U: np.ndarray
    alpha: np.ndarray        # per-gene dispersion
    tau: np.ndarray
    rho: np.ndarray          # per-cell contamination loading
    cell_types: list[str]
    genes: list[str]
    row_names: list[str]     # "<type>::<term>", one per random-effect row
    n_iter: int
    converged: bool
    params: dict[str, Any] = field(default_factory=dict)
    # What the readouts need to rebuild the fitted means, which the fit does
    # not store: R keeps the same pieces on `re_meta` and rebuilds mu on
    # demand rather than carrying an n x G matrix around. All of these are
    # small -- the designs are sparse and the rest are one value per cell.
    design: dict[str, Any] | None = None          # Z, the n x q random-effect design
    ambient_weights: dict[str, Any] | None = None  # W, the cross-type field
    fixed: np.ndarray | None = None               # the n x p fixed design
    offset: np.ndarray | None = None              # log library size
    type_code: np.ndarray | None = None           # each cell's 0-based type
    kernel: np.ndarray | None = None              # n x types, the Gaussian kernel
    x: np.ndarray | None = None                   # cell coordinates, for the plots
    y: np.ndarray | None = None
    condition_code: np.ndarray | None = None      # each cell's level, or None

    def __repr__(self) -> str:
        return (f"PaceFit({len(self.cell_types)} cell types, {len(self.genes)} genes, "
                f"{self.n_iter} iterations, "
                f"{'converged' if self.converged else 'NOT converged'})")


def _is_missing(value) -> bool:
    """Whether a label is absent: None, NaN, or a pandas NA.

    NaN is not equal to itself, which catches float NaN and NaT. `pandas.NA`
    goes one step further: comparing it to itself yields NA again, and taking
    the truth of that raises. An undefined self-comparison is exactly what a
    missing-value sentinel looks like, so that counts as missing too.
    """
    if value is None:
        return True
    try:
        return bool(value != value)
    except (TypeError, ValueError):
        return True


def _is_numeric_condition(condition):
    """Whether a condition arrived as numbers rather than labels."""
    values = np.asarray(condition)
    return values.dtype.kind in "biufc"


def _treatment_contrasts(labels, levels, name):
    """A factor as R's `model.matrix` codes it, under the default contrasts.

    `contr.treatment` drops the first level as the reference and gives one 0/1
    column per remaining level, named `paste0(column, level)`. The column name
    is the name of the variable, which is why the prefix has to be passed in:
    in R it comes from the colData column, and the terms it produces are what
    `resp_term` and every downstream filter match against.
    """
    if len(levels) < 2:
        raise ValueError(
            f"`{name}` has only the level {levels!r}; a condition needs at least two")
    dummies = np.zeros((len(labels), len(levels) - 1), order="F")
    for column, level in enumerate(levels[1:]):
        dummies[:, column] = (labels == level).astype(np.float64)
    return dummies, [f"{name}{level}" for level in levels[1:]]


def _neighbour_of_term(term, types):
    """The neighbour type a random-effect term carries a slope on, or None.

    R decides this with `grepl("(^|:)<type>$", term)`, so "Tumour" and
    "ResponderPR:Tumour" both slope on Tumour while "(Intercept)" and
    "ResponderPR" slope on nothing. Matching the string outright rather than by
    regex is the same answer for an ordinary name and a safer one for a type
    whose name holds a regex metacharacter, which in R can take another type's
    variance.

    Be aware that this is a deliberate disagreement, not only a safer one: for
    types named like "CD4.T" and "CD4_T", R's pattern for "CD4.T" also matches
    the term "CD4_T", so R gives that row the wrong type's kernel variance and
    the two languages then build different prior weights. The fitted slopes
    barely moved on the data this was measured on, but the designs are not the
    same.
    """
    for candidate in types:
        if term == candidate or term.endswith(":" + candidate):
            return candidate
    return None


def _csc(matrix) -> dict:
    """A scipy CSC matrix as the plain arrays the core takes."""
    m = matrix.tocsc()
    # tocsc() on something already CSC hands back the same object, so a matrix
    # carrying duplicate or unsorted entries arrives exactly as it is. The core
    # walks the stored entries one at a time, so a duplicated (cell, gene) is
    # read as two observations: a matrix elementwise identical to a good one,
    # with integer values and the right library sizes, fits 30% wrong. R cannot
    # reach this, because Matrix refuses such a dgCMatrix outright. Summing the
    # duplicates is what the format means, so normalise rather than refuse, and
    # on a copy because sum_duplicates() works in place.
    if not m.has_canonical_format:
        m = m.copy()
        m.sum_duplicates()          # sorts the row indices as part of its work
    # The core indexes sparse matrices with int, and numpy would narrow an
    # int64 index array to int32 silently: a wrapped column pointer reads past
    # the end rather than failing. scipy only widens the indices when it has
    # to, so this is reachable, and R cannot reach it because a dgCMatrix is
    # int-indexed from the start.
    limit = np.iinfo(np.int32).max
    if m.nnz > limit or max(m.shape) > limit:
        raise ValueError(
            f"this matrix needs 64-bit indices ({m.nnz} stored values, shape {m.shape}) "
            f"and the core takes 32-bit; the limit is {limit}")
    return {
        "indptr": np.ascontiguousarray(m.indptr, dtype=np.int32),
        "indices": np.ascontiguousarray(m.indices, dtype=np.int32),
        "data": np.asfortranarray(m.data, dtype=np.float64),
        "shape": (int(m.shape[0]), int(m.shape[1])),
    }


def _anchors(x, y, counts_csc, type_code, image_code, n_types, n_cells,
             homotypic_fraction=0.5, owner_threshold=0.1, core_threshold=0.1,
             n_threads=1):
    """Genes that anchor each cell type, from its homotypic-core cells.

    A cell surrounded by its own kind carries little borrowed signal, so the
    mean profile over those cells is the cleanest estimate of what the type
    really expresses. A type with too few such cells falls back to its overall
    mean rather than to a noisy one.
    """
    type_means = _core.group_column_means(counts_csc, type_code, n_types, False, n_threads)
    same = _core.same_type_fraction(x, y, type_code, image_code,
                                    radius=30.0, min_image_cells=50, n_threads=n_threads)
    core = same >= homotypic_fraction

    core_code = np.full(n_cells, -1, dtype=np.int32)
    core_code[core] = type_code[core]
    core_means = _core.group_column_means(counts_csc, core_code, n_types, False, n_threads)

    sizes = np.bincount(core_code[core_code >= 0], minlength=n_types)
    means = np.array(type_means, order="F")
    enough = sizes >= 20
    means[enough, :] = core_means[enough, :]

    out = _core.anchor_mask(np.asfortranarray(means), owner_threshold, core_threshold)
    return out["mask"], int(core.sum())


def pace_fit(counts, x, y, cell_type, image=None, *, genes=None, types=None,
             h_bio=30.0, h_tech=5.0, eps=None, condition=None,
             condition_name="condition", condition_levels=None, min_effective=30.0,
             within_image=True, kernel_per_image=False, homotypic_fraction=0.5,
             ambient_mode="cache", n_iter=32, min_iter=12, chunk_size=128,
             n_threads=1, verbose=False):
    """Fit PACE to a cells-by-genes count matrix.

    Parameters
    ----------
    counts
        Sparse cells by genes counts. Note the orientation: `AnnData.X` is
        already cells by genes, but a `SpatialExperiment` assay is the other
        way round.
    x, y
        Cell coordinates, in the same units as the bandwidths.
    cell_type
        One label per cell.
    types
        The cell types, in the order they should index the random-effect rows.
        Note that R builds its formula by pasting the names together, without
        backticks, so a name that is not a syntactic R name ("T-NK", "CD8+",
        "arm 2") fails there outright. Such a panel fits here but has no R
        counterpart to be compared against.
        None sorts the observed labels, which is Python's byte order and NOT
        R's: R sorts under `LC_COLLATE`, so a label set mixing cases orders
        differently ("epithelial" sorts before "T-NK" in R, after it here).
        The fit is equivariant under a permutation of types, so the science is
        unchanged, but the row order is not, and the two stop being
        bit-identical. Pass the types explicitly to match a given R fit.
    image
        One label per cell, or None for a single image.
    kernel_per_image
        Whether to build neighbourhoods within each image rather than over all
        cells at once. False, matching R: a cohort whose images are fields of
        view of one section wants the global search. Set it True when the
        images are separate samples. The contamination field is always
        per-image regardless, as it is in R.
    ambient_mode
        `"cache"` (default) materialises the whole cells-by-genes ambient
        product, as R's default does. `"stream"` keeps only the weight matrix
        and recomputes each chunk's columns from it, which is what lets a full
        transcriptome panel fit in memory: on a million cells by five thousand
        genes the cached product alone is several GB. The numbers are the same
        either way, since cache mode forms its product with the same routine
        and the two differ only in chunking.
    condition
        One label per cell, or None. Given one, the model gains a condition
        main effect in the fixed design and a condition-by-neighbour
        interaction in the random-effect block, which is the responder layer:
        `~ 1 + condition * (types)`, as in R.
    condition_name
        The name the condition's terms are built from, so a two-level
        condition called "Responder" with levels PD and PR gives the term
        `ResponderPR`. In R this is the colData column name; here it has to be
        passed, because an array has no name. `fit_anndata()` uses the obs key.
    condition_levels
        The levels in order. The first is the reference and is dropped, as R's
        default treatment contrasts drop it.

        This is NOT the same caveat as `types`, and the difference matters. A
        permutation of `types` relabels the output and leaves the science
        alone. The condition's first level is the arm everything else is
        contrasted against, so it decides what the baseline neighbour slopes
        mean and which way the interaction points. R picks it with
        `levels(factor(...))` under `LC_COLLATE`; leaving this None picks it by
        byte order, and the two disagree on labels as ordinary as
        {"Responder", "non-responder"}, where the fitted slopes then differ by
        2 rather than being reordered. Leaving it None warns for that reason.

    Returns
    -------
    PaceFit
    """
    if hasattr(counts, "obsm"):
        raise TypeError(
            "pace_fit() takes a count matrix, not an AnnData; use "
            "fit_anndata(adata, celltype_key) instead.")

    x = np.asfortranarray(np.asarray(x, dtype=np.float64).ravel())
    y = np.asfortranarray(np.asarray(y, dtype=np.float64).ravel())
    n_cells = x.size
    if n_cells == 0:
        raise ValueError("there are no cells to fit")

    labels = np.asarray(cell_type, dtype=object).ravel()
    # A missing label must not become a cell type. str(nan) is "nan", and if
    # the types came from the same str() it would be a type like any other,
    # with its own kernel column and random-effect block. R refuses here, and
    # so does this.
    missing = [i for i, v in enumerate(labels) if _is_missing(v)]
    if missing:
        raise ValueError(
            f"every cell needs a cell type; {len(missing)} are missing "
            f"(first at index {missing[0]})")
    labels = np.asarray([str(v) for v in labels], dtype=object)

    if types is None:
        types = sorted(set(labels))
    else:
        types = [str(t) for t in types]
        if len(set(types)) != len(types):
            raise ValueError("types must not repeat")
        unknown = sorted(set(labels) - set(types))
        if unknown:
            raise ValueError(f"these cell types are not in `types`: {unknown}")
    n_types = len(types)
    index = {t: i for i, t in enumerate(types)}
    type_code = np.ascontiguousarray([index[v] for v in labels], dtype=np.int32)

    condition_labels = None
    if condition is not None:
        if _is_numeric_condition(condition):
            raise ValueError(
                "a numeric or boolean condition is ambiguous: R would model such a "
                "column CONTINUOUSLY, giving one fixed column of the raw values, "
                "where this applies treatment contrasts and gives one dummy per "
                "level. The two fits differ. Pass the labels as strings to get the "
                "contrasts, and note that R's own responder term is derived wrongly "
                "for a continuous condition, so there is no reference to match.")
        condition_labels = np.asarray(condition, dtype=object).ravel()
        if condition_labels.size != n_cells:
            raise ValueError(
                f"condition has {condition_labels.size} entries for {n_cells} cells")
        if any(_is_missing(v) for v in condition_labels):
            raise ValueError("every cell needs a condition; none may be missing")
        condition_labels = np.asarray([str(v) for v in condition_labels], dtype=object)
        if condition_levels is None:
            condition_levels = sorted(set(condition_labels))
            # Unlike `types`, this is not a relabelling. The first level is the
            # reference the other arms are contrasted against, so it sets what
            # the baseline slopes mean and which way the interaction points, and
            # R picks it under LC_COLLATE where this picks it by byte order. On
            # labels like {"Responder", "non-responder"} the two disagree and
            # the fits are genuinely different, not permuted.
            warnings.warn(
                f"condition_levels was not given, so {condition_levels[0]!r} is the "
                f"reference arm (byte order of {condition_levels}). R chooses it by "
                "locale collation and may choose differently, which changes the fit "
                "rather than reordering it. Pass condition_levels to be certain.",
                stacklevel=2)
        else:
            condition_levels = [str(level) for level in condition_levels]
            if len(set(condition_levels)) != len(condition_levels):
                raise ValueError("condition_levels must not repeat")
            unknown = sorted(set(condition_labels) - set(condition_levels))
            if unknown:
                raise ValueError(f"these condition labels are not a level: {unknown}")
            # And the other direction. A level no cell carries gives an all-zero
            # dummy, which makes the fixed design singular and every coefficient
            # NaN, with nothing to say why.
            unused = [level for level in condition_levels
                      if level not in set(condition_labels)]
            if unused:
                raise ValueError(
                    f"no cell carries the condition level(s) {unused}, so they would "
                    "add an all-zero column and make the fit singular; drop them")

    if image is None:
        image_code = np.zeros(n_cells, dtype=np.int32)
        n_images = 1
    else:
        img = np.asarray(image, dtype=object).ravel()
        if any(_is_missing(v) for v in img):
            raise ValueError("every cell needs an image; none may be missing")
        img = np.asarray([str(v) for v in img], dtype=object)
        levels = sorted(set(img))
        lookup = {v: i for i, v in enumerate(levels)}
        image_code = np.ascontiguousarray([lookup[v] for v in img], dtype=np.int32)
        n_images = len(levels)

    counts_csc = _csc(counts)
    if counts_csc["shape"][0] != n_cells:
        raise ValueError(
            f"counts has {counts_csc['shape'][0]} rows but there are {n_cells} cells; "
            "counts must be cells by genes")
    n_genes = counts_csc["shape"][1]
    if genes is None:
        genes = [f"gene{i}" for i in range(n_genes)]
    elif len(genes) != n_genes:
        raise ValueError(f"genes has {len(genes)} names for {n_genes} columns of counts")

    # Library size, the offset the model carries at coefficient one. Summed
    # straight off the sparse matrix rather than through tocsr(), which would
    # copy the whole thing; the counts are integers, so the sum is exact
    # whatever order it accumulates in.
    library = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    if np.any(library <= 0):
        raise ValueError("every cell needs at least one count; drop empty cells first")
    offset = np.asfortranarray(np.log(library))

    if eps is None:
        eps = 3.0 * h_bio

    # --- kernels, then the two corrections that make the slopes interpretable
    # `image_code` is carried whether or not it is used: the core ignores the
    # grouping entirely when per_group is false, and searches one global grid,
    # which is what R does by default even with many images.
    kernels = _core.neighbour_kernels(x, y, type_code, n_types, image_code,
                                      per_group=kernel_per_image, h_bio=h_bio, h_tech=h_tech,
                                      eps=eps, n_threads=n_threads)
    k_bio = kernels["k_bio"]
    n_dropped = 0
    if min_effective > 0:
        dropped = _core.drop_sparse_kernel(k_bio, type_code, min_effective)
        k_bio, n_dropped = dropped["kernel"], dropped["n_dropped"]
    if within_image:
        k_bio = _core.centre_within_groups(k_bio, image_code, type_code, n_images, n_types)

    # --- the contamination field and what pins it
    weights_matrix = _core.ambient_field(x, y, type_code, image_code, n_images=n_images,
                                         h_tech=h_tech, edge_correct=True, n_threads=n_threads)
    # The solver reads the ambient PRODUCT, cells by genes, not the weight
    # matrix itself. Formed with the core's product because the summation order
    # is part of the answer; scipy's would differ in the last bits.
    # Cached, or streamed: see `ambient_mode`. Streaming hands the solver the
    # weight matrix itself and lets it form each chunk's columns, so the n x G
    # product is never materialised.
    if ambient_mode not in ("cache", "stream"):
        raise ValueError('ambient_mode must be "cache" or "stream"')
    field = (weights_matrix if ambient_mode == "stream"
             else _core.sparse_product(weights_matrix, counts_csc, n_threads))
    mask, n_core = _anchors(x, y, counts_csc, type_code, image_code, n_types, n_cells,
                            homotypic_fraction=homotypic_fraction, n_threads=n_threads)

    # --- the fixed design, and the random-effect term matrix
    # Without a condition this is R's `~ 1 + <types>` over an intercept-only
    # fixed design. With one it is `~ 1 + cond * (<types>)` over `~ 1 + cond`,
    # and the column order is R's: main effects in formula order, then the
    # interactions with the CONDITION LEVEL VARYING FASTEST within each type,
    # which is what model.matrix() emits and is not the obvious order.
    columns = [np.ones(n_cells)]
    term_names = ["(Intercept)"]
    if condition_labels is None:
        fixed = np.ones((n_cells, 1))
        fixed_names = ["(Intercept)"]
    else:
        dummies, dummy_names = _treatment_contrasts(
            condition_labels, condition_levels, condition_name)
        fixed = np.column_stack([np.ones(n_cells), dummies])
        fixed_names = ["(Intercept)"] + dummy_names
        columns.extend(dummies[:, j] for j in range(dummies.shape[1]))
        term_names.extend(dummy_names)
    columns.extend(k_bio[:, j] for j in range(n_types))
    term_names.extend(types)
    if condition_labels is not None:
        for j, neighbour in enumerate(types):
            for d, dummy_name in enumerate(dummy_names):
                columns.append(dummies[:, d] * k_bio[:, j])
                term_names.append(f"{dummy_name}:{neighbour}")

    terms = np.asfortranarray(np.column_stack(columns))
    fixed = np.asfortranarray(fixed)
    p_fixed = fixed.shape[1]
    n_terms = terms.shape[1]
    design = _core.random_design_block(terms, type_code, n_types)
    q = n_types * n_terms
    # How many random-effect rows have no data behind them at all: their design
    # column is identically zero, so the row is identified only by the prior and
    # comes back as whatever the prior says. A cell type absent from one
    # condition arm produces a lot of these, and nothing else in the output
    # distinguishes them from a genuine estimate of zero.
    column_pointer = np.asarray(design["indptr"])
    n_structurally_zero = int(np.sum(np.diff(column_pointer) == 0))

    # Row (type, term) of U, in the order the design lays them out.
    row_names = [f"{types[g]}::{term_names[t]}" for t in range(n_terms) for g in range(n_types)]
    # ONE-based, like R's match(): zero means "this row belongs to no cell-type
    # block" and the weight stays at one. Passing zero-based codes silently
    # gives every row weight one, which looks like a fit with no prior at all.
    focal_of_row = np.ascontiguousarray(
        np.tile(np.arange(1, n_types + 1, dtype=np.int32), n_terms))

    # --- prior weights: how well a gene is seen in the focal type, scaled by
    # how much the neighbour kernel actually varies over that focal's cells. A
    # slope estimated where the neighbour density barely moves gets little
    # prior weight; an intercept term takes the detection unscaled.
    sizes = np.bincount(type_code, minlength=n_types)
    detection = np.array(
        _core.group_column_means(counts_csc, type_code, n_types, True, n_threads), order="F")
    # A cell type with no cells has no detection rate, and the core says so with
    # NaN. R starts the matrix at zero and fills only the types that have cells,
    # so an empty type takes zero and the prior weights stay finite. Letting the
    # NaN through poisons the whole weight matrix, and from there every
    # coefficient in the fit: an empty type is reachable whenever `types` is
    # passed explicitly, which is exactly what matching an R row order asks for.
    detection[sizes == 0, :] = 0.0
    covariance = _core.group_covariances(np.asfortranarray(k_bio), type_code, n_types, n_threads)
    k_var = np.full((n_types, n_types), np.nan)
    slab = n_types * n_types
    for f in range(n_types):
        if sizes[f] < 5:        # too few cells to carry a usable variance
            continue
        block_cov = covariance[f * slab:(f + 1) * slab].reshape(n_types, n_types, order="F")
        k_var[f, :] = np.diag(block_cov)
    k_var_norm = _core.normalise_rows_to_max(np.asfortranarray(k_var))

    # Zero, as R does: a row outside every cell-type block takes no scale. With
    # the single block here every row is assigned below, so this only matters
    # once an image-level block is added.
    scale_of_row = np.zeros(q)
    for t, term in enumerate(term_names):
        # Which neighbour this term slopes on, by name rather than by position:
        # with a condition the terms are no longer intercept-then-types, and
        # "ResponderPR:Tumour" has to take Tumour's variance. A term that
        # slopes on no neighbour, the intercept or a condition main effect,
        # takes the detection profile unscaled.
        neighbour = _neighbour_of_term(term, types)
        scale = 1.0 if neighbour is None else None
        column = None if neighbour is None else types.index(neighbour)
        for f in range(n_types):
            scale_of_row[t * n_types + f] = scale if column is None else k_var_norm[f, column]
    weights = _core.data_informed_weights(detection, focal_of_row,
                                          np.asfortranarray(scale_of_row))

    if verbose:
        print(f"pace: {n_cells} cells x {n_genes} genes, {n_types} types, {n_images} image(s); "
              f"{n_dropped} sparse pairs dropped; {n_core} homotypic-core cells")
        if n_structurally_zero:
            print(f"pace: {n_structurally_zero} of {q} random-effect rows have an "
                  "all-zero design column and are identified only by the prior")

    # --- the fit
    # Mind the indexing: the solver reads `cell_group` ONE-based, while
    # `cells_by_group` and `group_start` are zero-based offsets, and
    # random_design_block() above took the zero-based codes. R does the same,
    # and the solver refuses an out-of-range index rather than reading past the
    # end, which is how this showed up.
    block = {
        "col_offset": 0, "n_terms": n_terms, "n_groups": n_types, "n_cols": q,
        "terms": terms,
        "cell_group": np.ascontiguousarray(type_code + 1, dtype=np.int32),
        "group_start": design["group_start"],
        "cells_by_group": design["cells_by_group"],
    }
    fitted = _core.fit_irls(
        x_fixed=fixed, p=p_fixed,
        z_design={"indptr": design["indptr"], "indices": design["indices"],
                  "data": design["data"], "shape": design["shape"]},
        blocks=[block], counts=counts_csc,
        ambient={"indptr": field["indptr"], "indices": field["indices"],
                 "data": field["data"], "shape": field["shape"]},
        offset=offset,
        alpha_init=np.asfortranarray(np.ones(n_genes)),
        mask=mask, mask_index=np.ascontiguousarray(type_code + 1, dtype=np.int32),
        n_mask_rows=n_types, data_informed_weights=weights,
        n_cells=n_cells, q_total=q, n_genes=n_genes,
        options=dict(n_iter=n_iter, min_iter=min_iter, chunk_size=chunk_size,
                     n_threads=n_threads, verbose=verbose,
                     ambient_streamed=ambient_mode == "stream"),
    )

    # R raises every one of these as a warning. Dropping them is what made a
    # NaN fit silent.
    for text in fitted.get("warnings", []):
        warnings.warn(f"pace: {text}", stacklevel=2)

    # R's final pass refuses a fit whose linear predictor is not finite
    # ("eta must be finite", src/core/statistics.cpp), so a degenerate design
    # errors there rather than returning anything. Python has no final pass to
    # catch it, so check the coefficients directly: without this a fit that
    # diverged comes back as a PaceFit full of NaN, which reads like a result.
    # The standard errors are exempt, since they are deliberately NaN until the
    # last iteration writes them.
    bad_genes = ~(np.all(np.isfinite(fitted["B"]), axis=0)
                  & np.all(np.isfinite(fitted["U"]), axis=0))
    if np.any(bad_genes):
        first = int(np.argmax(bad_genes))
        raise FloatingPointError(
            f"the fit did not converge to finite coefficients for "
            f"{int(bad_genes.sum())} of {n_genes} genes (first: {genes[first]!r}). "
            "This is a degenerate design rather than a numerical accident: the "
            "usual causes are a cell type or condition level that no cell "
            "carries, a condition level perfectly aligned with a cell type, or "
            "counts that are not counts. R refuses the same fit.")

    return PaceFit(
        B=fitted["B"], U=fitted["U"], se_B=fitted["se_B"], se_U=fitted["se_U"],
        alpha=fitted["alpha"], tau=fitted["tau_flat"], rho=fitted["rho"],
        cell_types=types, genes=list(genes), row_names=row_names,
        n_iter=int(fitted["n_iter"]), converged=bool(fitted["converged"]),
        design={"indptr": design["indptr"], "indices": design["indices"],
                "data": design["data"], "shape": design["shape"]},
        ambient_weights=weights_matrix, fixed=fixed, offset=offset, type_code=type_code,
        kernel=k_bio, x=x, y=y,
        condition_code=None if condition_labels is None else np.ascontiguousarray(
            [condition_levels.index(v) for v in condition_labels], dtype=np.int32),
        params=dict(h_bio=h_bio, h_tech=h_tech, eps=eps, min_effective=min_effective,
                    within_image=within_image, kernel_per_image=kernel_per_image,
                    ambient_mode=ambient_mode,
                    homotypic_fraction=homotypic_fraction,
                    n_images=n_images, n_dropped=n_dropped, n_core_cells=n_core,
                    n_structurally_zero=n_structurally_zero,
                    condition_name=condition_name if condition_labels is not None else None,
                    condition_levels=condition_levels, fixed_names=fixed_names,
                    term_names=term_names),
    )


# ---------------------------------------------------------------------------
# AnnData
# ---------------------------------------------------------------------------
# `anndata` is not imported, and is not a dependency: everything below needs
# only the attributes an AnnData exposes, so anything shaped like one works.


def _anndata_counts(adata, layer):
    """The count matrix, as compressed-column sparse, having checked it is counts.

    This is the one check R never needs. `counts(spe)` names the assay it
    wants, whereas `adata.X` is whatever the last preprocessing step left
    there, very often normalised or log1p-transformed. PACE is a negative
    binomial model on raw counts, so fitting those would converge, report
    slopes, and mean nothing at all.
    """
    if layer is None:
        matrix, source = adata.X, "adata.X"
    else:
        if layer not in adata.layers:
            raise KeyError(
                f'adata.layers has no "{layer}"; it has {list(adata.layers)}')
        matrix, source = adata.layers[layer], f'adata.layers["{layer}"]'
    if matrix is None:
        raise ValueError(f"{source} is empty")
    if getattr(adata, "isbacked", False):
        raise ValueError(
            "this AnnData is backed, so the counts are on disk rather than in "
            "memory. A backed sparse X is not a scipy matrix and would be "
            "misreported as the wrong shape, and a backed dense X would be read "
            "in whole without saying so. Call adata.to_memory() first, or pass "
            "the counts to pace_fit() yourself.")
    if hasattr(matrix, "compute") and not issparse(matrix):
        raise ValueError(
            f"{source} is a lazy array. Materialising it here would hide how "
            "much memory the fit needs, so call .compute() (or .persist()) "
            "first and assign the result back.")

    if issparse(matrix):
        counts = matrix.tocsc()
        values = counts.data
    else:
        values = np.asarray(matrix)
        if values.ndim != 2:
            raise ValueError(f"{source} must be two-dimensional, cells by genes")
        counts = csc_matrix(values)
        values = counts.data

    if values.size:
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{source} holds non-finite values")
        if values.min() < 0:
            raise ValueError(
                f"{source} holds negative values, so it is not counts. PACE fits "
                "counts; point `layer` at the raw ones.")
        if not np.array_equal(values, np.rint(values)):
            raise ValueError(
                f"{source} is not integer-valued, so it looks normalised or "
                "log-transformed rather than raw counts, and a negative binomial "
                "fit of it would be meaningless. Pass the raw counts instead, "
                'usually layer="counts". For counts kept in adata.raw, hand '
                "pace_fit() adata.raw.X with genes=adata.raw.var_names. If the "
                "values really are counts and merely non-integer, pace_fit() "
                "takes the arrays without this check.")
    return counts


def _anndata_coordinates(adata, spatial_key):
    """The two spatial coordinates, from `obsm`."""
    if spatial_key not in adata.obsm:
        raise KeyError(
            f'adata.obsm has no "{spatial_key}"; it has {list(adata.obsm)}. '
            "Spatial coordinates are expected there, as cells by two columns, "
            "which is where squidpy and the readers put them.")
    coordinates = np.asarray(adata.obsm[spatial_key], dtype=np.float64)
    if coordinates.ndim != 2:
        raise ValueError(f'adata.obsm["{spatial_key}"] must be two-dimensional')
    if coordinates.shape[1] != 2:
        raise ValueError(
            f'adata.obsm["{spatial_key}"] has {coordinates.shape[1]} columns and '
            "PACE's kernels are two-dimensional. Rather than drop a dimension "
            "quietly, add a two-column key, for instance "
            f'adata.obsm["spatial_2d"] = adata.obsm["{spatial_key}"][:, :2], and '
            'pass spatial_key="spatial_2d".')
    if not np.all(np.isfinite(coordinates)):
        bad = int(np.argmax(~np.isfinite(coordinates).all(axis=1)))
        raise ValueError(
            f'adata.obsm["{spatial_key}"] holds non-finite coordinates, first at '
            f"cell {bad}. The core refuses them too, but without naming the key.")
    return coordinates[:, 0], coordinates[:, 1]


def _anndata_column(adata, key, what):
    """One `obs` column as objects, keeping a missing value missing.

    Deliberately not `.astype(str)`: that renders a missing label as the string
    "nan", which would then be a cell type like any other.
    """
    if key not in adata.obs:
        raise KeyError(f'adata.obs has no "{key}" for the {what}')
    return np.asarray(adata.obs[key], dtype=object)


def _categorical_levels(column):
    """The levels R would use for a condition, or None to let `pace_fit` sort.

    R takes them from `levels(factor(df[[condition_col]]))`, and `factor()` on
    something already a factor keeps its order while dropping unused levels. A
    pandas Categorical is the same idea, so its category order is honoured;
    a plain column of strings has no order to honour and gets sorted.

    An unused category is dropped here, which follows `.pace_prepare()` but not
    `pace_fit_streaming()`: R hands the untouched factor to `model.matrix`,
    which keeps the unused level as an all-zero column and makes the design
    singular. R is inconsistent with itself on this, so the two cannot both be
    matched; dropping is the one that fits.
    """
    if not hasattr(column, "cat"):
        return None
    observed = set(column.dropna().astype(str))
    return [str(level) for level in column.cat.categories if str(level) in observed]


def fit_anndata(adata, celltype_key, *, image_key=None, condition_key=None,
                spatial_key="spatial", layer=None, **kwargs):
    """Fit PACE to an `AnnData`.

    A reader over `pace_fit()`, which does the work. Coordinates come from
    `obsm`, labels from `obs`, gene names from `var_names`, and the counts from
    `X` or a layer.

    Parameters
    ----------
    adata
        An `AnnData`, or anything exposing `X`, `obs`, `obsm` and `var_names`.
        It is never modified.
    celltype_key
        The `obs` column holding the cell type labels.
    image_key
        The `obs` column grouping cells into images or samples, or None for a
        single image.
    condition_key
        The `obs` column holding the condition, or None. Its name becomes the
        prefix of the condition's terms, exactly as the colData column name
        does in R, so a column "Responder" with levels PD and PR gives the term
        `ResponderPR`. A Categorical column's category order is used as the
        level order, as R uses a factor's; a column of plain strings is sorted.
    spatial_key
        The `obsm` key holding the coordinates, cells by two columns. Defaults
        to `"spatial"`, the squidpy convention.
    layer
        The layer holding raw counts, or None for `X`. Either way the values
        are checked for being counts, because `X` is often normalised.

        A dense matrix is converted to sparse here, which costs about six times
        its own size at peak, scipy building a coordinate form on the way. A
        large cohort is better stored sparse, which costs nothing to read.
    **kwargs
        Passed to `pace_fit()`: `types`, `kernel_per_image`, `h_bio`, `h_tech`
        and the rest. Read its docstring for the two places Python's defaults
        are not R's.

    Returns
    -------
    PaceFit
    """
    for attribute in ("X", "obs", "obsm", "var_names"):
        if not hasattr(adata, attribute):
            raise TypeError(
                f"this does not look like an AnnData: no `{attribute}`. To fit "
                "plain arrays, call pace_fit() instead.")

    # The cheap lookups first: a missing obsm key should not be reported only
    # after the whole count matrix has been converted.
    x, y = _anndata_coordinates(adata, spatial_key)
    cell_type = _anndata_column(adata, celltype_key, "cell types")
    image = _anndata_column(adata, image_key, "images") if image_key is not None else None
    if condition_key is not None:
        condition = _anndata_column(adata, condition_key, "condition")
        kwargs.setdefault("condition_name", condition_key)
        levels = _categorical_levels(adata.obs[condition_key])
        if levels is not None:
            kwargs.setdefault("condition_levels", levels)
    else:
        condition = None

    # Gene names, as strings. Duplicates are refused for the reason R gives
    # when it refuses them: two different genes would share one row of the
    # statistics and one row of the output, and nothing downstream could tell
    # them apart.
    genes = [str(name) for name in adata.var_names]
    if len(set(genes)) != len(genes):
        repeated = sorted({name for name in genes if genes.count(name) > 1})
        raise ValueError(
            f"var_names must be unique; {len(repeated)} are repeated "
            f"(first {repeated[:3]}). adata.var_names_make_unique() fixes it.")

    counts = _anndata_counts(adata, layer)
    if counts.shape[1] == 0:
        raise ValueError("there are no genes to fit")
    if counts.shape[0] != x.size:
        raise ValueError(
            f"the counts have {counts.shape[0]} rows and there are {x.size} sets of "
            "coordinates; an AnnData is cells by genes, so these should agree")

    return pace_fit(counts, x, y, cell_type, image=image, condition=condition,
                    genes=genes, **kwargs)
