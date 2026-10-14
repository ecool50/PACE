"""Sequencing the stages of a PACE fit.

Every numerical stage lives in the compiled core and is shared with the R
package. What this module does is order them, build the designs, and hand the
solver what it expects -- the job `pace_fit_streaming()` does in R.

`fit_anndata()` is a thin reader on top of `pace_fit()` for callers holding
an AnnData rather than arrays.

Only the no-condition path is implemented. A condition term brings in R's
model.matrix semantics (contrasts, interaction expansion, level ordering), and
guessing at those is how a translation diverges quietly, so `condition` raises
rather than approximating.
"""
from __future__ import annotations

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


def _csc(matrix) -> dict:
    """A scipy CSC matrix as the plain arrays the core takes."""
    m = matrix.tocsc()
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
             h_bio=30.0, h_tech=5.0, eps=None, condition=None, min_effective=30.0,
             within_image=True, kernel_per_image=False, homotypic_fraction=0.5,
             n_iter=32, min_iter=12, chunk_size=128, n_threads=1, verbose=False):
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
    condition
        Not implemented. Raises, rather than approximating R's model.matrix.

    Returns
    -------
    PaceFit
    """
    if condition is not None:
        raise NotImplementedError(
            "a condition term is not implemented yet. It needs R's model.matrix "
            "semantics -- contrasts, interaction expansion, level ordering -- and "
            "approximating those would make the fit disagree with the R package "
            "in ways that are hard to see.")

    if hasattr(counts, "obsm"):
        raise TypeError(
            "pace_fit() takes a count matrix, not an AnnData; use "
            "fit_anndata(adata, celltype_key) instead.")

    x = np.asfortranarray(np.asarray(x, dtype=np.float64).ravel())
    y = np.asfortranarray(np.asarray(y, dtype=np.float64).ravel())
    n_cells = x.size

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
    # NOTE: this materialises the whole n x G product, as R's default
    # `ambient_mode = "cache"` does. R can also stream it, which is what makes a
    # full transcriptome panel fit in memory; that mode is not exposed here, so
    # a cohort of ~1M cells by ~5,000 genes needs several GB for this one array
    # on top of the counts.
    field = _core.sparse_product(weights_matrix, counts_csc, n_threads)
    mask, n_core = _anchors(x, y, counts_csc, type_code, image_code, n_types, n_cells,
                            homotypic_fraction=homotypic_fraction, n_threads=n_threads)

    # --- the random-effect design: an intercept plus one slope per neighbour type
    terms = np.asfortranarray(np.column_stack([np.ones(n_cells), k_bio]))
    n_terms = terms.shape[1]
    design = _core.random_design_block(terms, type_code, n_types)
    q = n_types * n_terms

    # Row (type, term) of U, in the order the design lays them out.
    term_names = ["(Intercept)"] + list(types)
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
    detection = _core.group_column_means(counts_csc, type_code, n_types, True, n_threads)
    covariance = _core.group_covariances(np.asfortranarray(k_bio), type_code, n_types, n_threads)
    k_var = np.full((n_types, n_types), np.nan)
    sizes = np.bincount(type_code, minlength=n_types)
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
    for t in range(n_terms):
        for f in range(n_types):
            row = t * n_types + f
            # term 0 is the intercept and carries no neighbour; the rest are
            # the slope on neighbour type (t - 1)
            scale_of_row[row] = 1.0 if t == 0 else k_var_norm[f, t - 1]
    weights = _core.data_informed_weights(detection, focal_of_row,
                                          np.asfortranarray(scale_of_row))

    if verbose:
        print(f"pace: {n_cells} cells x {n_genes} genes, {n_types} types, {n_images} image(s); "
              f"{n_dropped} sparse pairs dropped; {n_core} homotypic-core cells")

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
        x_fixed=np.asfortranarray(np.ones((n_cells, 1))), p=1,
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
                     n_threads=n_threads, verbose=verbose),
    )

    return PaceFit(
        B=fitted["B"], U=fitted["U"], se_B=fitted["se_B"], se_U=fitted["se_U"],
        alpha=fitted["alpha"], tau=fitted["tau_flat"], rho=fitted["rho"],
        cell_types=types, genes=list(genes), row_names=row_names,
        n_iter=int(fitted["n_iter"]), converged=bool(fitted["converged"]),
        params=dict(h_bio=h_bio, h_tech=h_tech, eps=eps, min_effective=min_effective,
                    within_image=within_image, kernel_per_image=kernel_per_image,
                    homotypic_fraction=homotypic_fraction,
                    n_images=n_images, n_dropped=n_dropped, n_core_cells=n_core),
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
                'usually layer="counts" or adata.raw. If the values really are '
                "counts and merely non-integer, call pace_fit() with the arrays "
                "directly, which does not check.")
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
    return coordinates[:, 0], coordinates[:, 1]


def _anndata_column(adata, key, what):
    """One `obs` column as objects, keeping a missing value missing.

    Deliberately not `.astype(str)`: that renders a missing label as the string
    "nan", which would then be a cell type like any other.
    """
    if key not in adata.obs:
        raise KeyError(f'adata.obs has no "{key}" for the {what}')
    return np.asarray(adata.obs[key], dtype=object)


def fit_anndata(adata, celltype_key, *, image_key=None, spatial_key="spatial",
                layer=None, **kwargs):
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
    spatial_key
        The `obsm` key holding the coordinates, cells by two columns. Defaults
        to `"spatial"`, the squidpy convention.
    layer
        The layer holding raw counts, or None for `X`. Either way the values
        are checked for being counts, because `X` is often normalised.
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

    counts = _anndata_counts(adata, layer)
    x, y = _anndata_coordinates(adata, spatial_key)
    cell_type = _anndata_column(adata, celltype_key, "cell types")
    image = _anndata_column(adata, image_key, "images") if image_key is not None else None

    if counts.shape[0] != x.size:
        raise ValueError(
            f"the counts have {counts.shape[0]} rows and there are {x.size} sets of "
            "coordinates; an AnnData is cells by genes, so these should agree")

    return pace_fit(counts, x, y, cell_type, image=image,
                    genes=list(adata.var_names), **kwargs)
