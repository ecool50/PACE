"""The fitted means, and the per-cell-type statistics the readouts are built on.

A fit does not keep its fitted means. On a real cohort `mu` is a dense cells by
genes matrix, which is the largest thing in sight and is recoverable from what
the fit does keep, so R rebuilds it on demand and so does this. What the
reporting layer actually reads is not `mu` itself but two small summaries of
it per cell type, so the rebuild runs a block of genes at a time and keeps only
those: memory is cells by block, never cells by genes.

The arithmetic is the solver's own final pass, and deliberately the same
ambient product: taking it from scipy instead would put a second
implementation on the one path the fit's bit-identity rests on.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import csc_matrix

from . import _core

MU_FLOOR = 1e-6          # the solver's floor on the fitted mean


@dataclass
class FitStatistics:
    """Per-(cell type, gene) summaries of the fitted means."""

    cell_types: list[str]
    genes: list[str]
    n_cells: np.ndarray      # cells per type
    mu_mean: np.ndarray      # types x genes, mean fitted mean
    toff_var: np.ndarray     # types x genes, variance of the contamination offset
    has_contamination: bool

    def __repr__(self) -> str:
        return (f"FitStatistics({len(self.cell_types)} cell types, {len(self.genes)} genes, "
                f"contamination {'present' if self.has_contamination else 'absent'})")


def _as_csc(arrays):
    return csc_matrix((arrays["data"], arrays["indices"], arrays["indptr"]),
                      shape=arrays["shape"])


def fitted_means(fit, counts, genes=None):
    """Rebuild the fitted means for a block of genes.

    Returns `mu_bio` and `mu_spill`, cells by the genes asked for, exactly as
    the solver's final pass formed them. `mu_spill` is None for a fit with no
    contamination term.
    """
    for name in ("design", "fixed", "offset"):
        if getattr(fit, name) is None:
            raise ValueError(
                "this fit was not made by pace_fit() in this version, so it does "
                f"not carry `{name}` and the fitted means cannot be rebuilt")
    columns = np.arange(len(fit.genes)) if genes is None else np.asarray(genes, dtype=int)

    # mu_bio = exp(eta + offset), eta = X B + Z U, as pace_fit_streaming forms it
    eta = fit.fixed @ fit.B[:, columns] + _as_csc(fit.design) @ fit.U[:, columns]
    mu_bio = np.maximum(np.exp(eta + fit.offset[:, None]), MU_FLOOR)
    if fit.ambient_weights is None or fit.rho is None:
        return mu_bio, None

    # The ambient product through the core, not scipy: the summation order is
    # part of the answer, and the fit's own path goes through here.
    block = counts.tocsc()[:, columns]
    ambient = _as_csc(_core.sparse_product(
        fit.ambient_weights,
        {"indptr": np.ascontiguousarray(block.indptr, dtype=np.int32),
         "indices": np.ascontiguousarray(block.indices, dtype=np.int32),
         "data": np.asfortranarray(block.data, dtype=np.float64),
         "shape": (int(block.shape[0]), int(block.shape[1]))}, 1)).toarray()
    mu_spill = np.maximum(ambient * fit.rho[:, None], 0.0)
    return mu_bio, mu_spill


def fit_statistics(fit, counts, *, max_block_entries=5_000_000):
    """The per-cell-type statistics, rebuilding the means a gene block at a time.

    `max_block_entries` caps cells times genes per block, so the peak is that
    many doubles however many genes the panel has.
    """
    n_cells = fit.offset.size
    n_genes = len(fit.genes)
    n_types = len(fit.cell_types)
    per_block = max(1, int(max_block_entries // max(n_cells, 1)))

    sizes = np.bincount(fit.type_code, minlength=n_types)
    mu_mean = np.zeros((n_types, n_genes))
    toff_var = np.full((n_types, n_genes), np.nan)
    any_contamination = False

    for start in range(0, n_genes, per_block):
        columns = np.arange(start, min(start + per_block, n_genes))
        mu_bio, mu_spill = fitted_means(fit, counts, columns)
        if mu_spill is None:
            total, offset = mu_bio, None
        else:
            total = np.maximum(mu_bio + mu_spill, MU_FLOOR)
            # the contamination's share, on the scale the model adds it
            offset = np.log1p(mu_spill / mu_bio)
            any_contamination = any_contamination or bool(np.any(mu_spill > 0))
        for t in range(n_types):
            cells = fit.type_code == t
            if not cells.any():
                continue
            mu_mean[t, columns] = total[cells, :].mean(axis=0)
            if offset is not None:
                toff_var[t, columns] = offset[cells, :].var(axis=0, ddof=1)

    return FitStatistics(cell_types=list(fit.cell_types), genes=list(fit.genes),
                         n_cells=sizes, mu_mean=mu_mean, toff_var=toff_var,
                         has_contamination=any_contamination)
