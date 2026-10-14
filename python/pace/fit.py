"""Sequencing the stages of a PACE fit.

Every numerical stage lives in the compiled core and is shared with the R
package. What this module does is order them, build the designs, and hand the
solver what it expects -- the job `pace_fit_streaming()` does in R.

Only the no-condition path is implemented. A condition term brings in R's
model.matrix semantics (contrasts, interaction expansion, level ordering), and
guessing at those is how a translation diverges quietly, so `condition` raises
rather than approximating.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

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


def _csc(matrix) -> dict:
    """A scipy CSC matrix as the plain arrays the core takes."""
    m = matrix.tocsc()
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


def pace_fit(counts, x, y, cell_type, image=None, *, genes=None, h_bio=30.0, h_tech=5.0,
             eps=None, condition=None, min_effective=30.0, within_image=True,
             homotypic_fraction=0.5, n_iter=32, min_iter=12, chunk_size=128,
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
        One label per cell. The fitting order is the sorted unique set, which
        is what R uses when `types` is not given.
    image
        One label per cell, or None for a single image.
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

    x = np.asfortranarray(np.asarray(x, dtype=np.float64).ravel())
    y = np.asfortranarray(np.asarray(y, dtype=np.float64).ravel())
    n_cells = x.size

    labels = np.asarray(cell_type, dtype=object).ravel()
    types = sorted({str(v) for v in labels})
    n_types = len(types)
    index = {t: i for i, t in enumerate(types)}
    if any(v is None or str(v) not in index for v in labels):
        raise ValueError("every cell needs a cell type; none may be missing")
    type_code = np.ascontiguousarray([index[str(v)] for v in labels], dtype=np.int32)

    if image is None:
        image_code = np.zeros(n_cells, dtype=np.int32)
        n_images = 1
    else:
        img = np.asarray(image, dtype=object).ravel()
        levels = sorted({str(v) for v in img})
        lookup = {v: i for i, v in enumerate(levels)}
        image_code = np.ascontiguousarray([lookup[str(v)] for v in img], dtype=np.int32)
        n_images = len(levels)

    counts_csc = _csc(counts)
    if counts_csc["shape"][0] != n_cells:
        raise ValueError(
            f"counts has {counts_csc['shape'][0]} rows but there are {n_cells} cells; "
            "counts must be cells by genes")
    n_genes = counts_csc["shape"][1]
    if genes is None:
        genes = [f"gene{i}" for i in range(n_genes)]

    # Library size, the offset the model carries at coefficient one.
    library = np.asarray(counts.tocsr().sum(axis=1)).ravel().astype(np.float64)
    if np.any(library <= 0):
        raise ValueError("every cell needs at least one count; drop empty cells first")
    offset = np.asfortranarray(np.log(library))

    if eps is None:
        eps = 3.0 * h_bio

    # --- kernels, then the two corrections that make the slopes interpretable
    kernels = _core.neighbour_kernels(x, y, type_code, n_types, image_code,
                                      per_group=n_images > 1, h_bio=h_bio, h_tech=h_tech,
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

    scale_of_row = np.ones(q)
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
                    within_image=within_image, homotypic_fraction=homotypic_fraction,
                    n_images=n_images, n_dropped=n_dropped, n_core_cells=n_core),
    )
