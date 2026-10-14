"""The per-gene variance decomposition.

For one focal cell type and one gene, the fitted variance splits into five
blocks: the cell type's own offset, the baseline spatial state (the neighbour
slopes), the responder spatial state (their interaction with a condition), the
contamination, and what the negative binomial dispersion accounts for. The
arithmetic is `pace::variance_decomposition` in the shared core, so this module
only assembles its twenty-odd inputs, which is the whole difficulty: the index
vectors are one-based and point into rows of `U`, and getting one wrong gives a
plausible table built from the wrong coefficients.

Mirrors `mvpql_variance_decomposition_stats()` in R, which is where the
assembly comes from.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import _core
from .stats import fit_statistics

BLOCKS = ("Cell type %", "Spatial state %", "Responder spatial state %",
          "Spillover %", "Residual %")
COMPONENTS = ("celltype_offset_sq", "V_state_baseline", "V_state_responder",
              "V_spill", "V_disp")


@dataclass
class Decomposition:
    """The per-(focal, gene) table, and its per-focal aggregates."""

    gene: np.ndarray
    focal: np.ndarray
    n_focal: np.ndarray
    v_state_baseline: np.ndarray
    v_state_responder: np.ndarray
    v_spill: np.ndarray
    v_disp: np.ndarray
    celltype_offset_sq: np.ndarray
    focal_mean: np.ndarray
    max_other_mean: np.ndarray
    spec: np.ndarray
    focal_other_ratio: np.ndarray
    is_contaminated: np.ndarray
    pct: np.ndarray             # rows x 5, in BLOCKS order
    focal_levels: list[str]
    mean_table: np.ndarray      # focals x 5
    specw_table: np.ndarray
    pooled_table: np.ndarray
    total_ss: np.ndarray
    n_rows: np.ndarray
    n_specific: np.ndarray

    def __len__(self) -> int:
        return self.gene.size

    def __repr__(self) -> str:
        return (f"Decomposition({len(self)} gene-focal rows, "
                f"{len(self.focal_levels)} focal types)")

    def to_dataframe(self):
        """The per-(focal, gene) table, with R's column names."""
        import pandas as pd
        frame = pd.DataFrame({
            "gene": self.gene, "focal": self.focal, "n_focal": self.n_focal,
            "V_state_baseline": self.v_state_baseline,
            "V_state_responder": self.v_state_responder,
            "V_spill": self.v_spill, "V_disp": self.v_disp,
            "celltype_offset_sq": self.celltype_offset_sq,
            "focal_mean": self.focal_mean, "max_other_mean": self.max_other_mean,
            "spec": self.spec, "focal_other_ratio": self.focal_other_ratio,
            "is_contaminated": self.is_contaminated})
        for i, name in enumerate(BLOCKS):
            frame[name] = self.pct[:, i]
        return frame

    def aggregates(self, weighting="specw"):
        """The per-focal table, as a DataFrame.

        `weighting` picks the plain mean of each block's percentage, the
        specificity-weighted mean, or the pooled share.
        """
        import pandas as pd
        table = {"mean": self.mean_table, "specw": self.specw_table,
                 "pooled": self.pooled_table}[weighting]
        frame = pd.DataFrame(table, columns=list(BLOCKS))
        frame.insert(0, "focal", self.focal_levels)
        frame["n_genes"] = self.n_rows
        frame["n_specific"] = self.n_specific
        return frame


def _term_index(fit):
    """Each random-effect term's one-based index, from the fit's row names."""
    seen = {}
    for name in fit.row_names:
        _, _, term = name.partition("::")
        if term not in seen:
            seen[term] = len(seen) + 1
    return seen


def pace_decompose(fit, counts, *, resp_term=None, stats=None, nb1=True, n_threads=1):
    """Decompose each gene's fitted variance, per focal cell type.

    Parameters
    ----------
    fit
        A `PaceFit`.
    counts
        The same cells-by-genes counts the fit was made from.
    resp_term
        Which condition term's interaction counts as the responder block. None
        takes R's default, the condition's last level. Ignored without a
        condition.
    stats
        Pre-computed `fit_statistics()`, if you already have them.

    Returns
    -------
    Decomposition
    """
    for name in ("kernel", "type_code"):
        if getattr(fit, name) is None:
            raise ValueError(
                f"this fit does not carry `{name}`, so it cannot be decomposed; "
                "re-fit with this version of pace_fit()")

    types = list(fit.cell_types)
    n_types = len(types)
    n_genes = len(fit.genes)
    if stats is None:
        stats = fit_statistics(fit, counts, )
    terms = _term_index(fit)

    # Rows of U are laid out (term - 1) * n_types + focal, one-based, which is
    # R's ct_col(). Every index handed to the core is one-based for that reason.
    def row_of(term_name, focal_index):      # focal_index is one-based
        return (terms[term_name] - 1) * n_types + focal_index

    focal_group = np.arange(1, n_types + 1, dtype=np.int32)
    group_size = np.bincount(fit.type_code, minlength=n_types).astype(np.int32)

    slope_rows = np.asfortranarray(
        [[row_of(neighbour, f) for f in focal_group] for neighbour in types],
        dtype=np.int32)
    intercept_rows = np.ascontiguousarray(
        [row_of("(Intercept)", f) for f in focal_group], dtype=np.int32)

    ct_means = np.asfortranarray(
        _core.group_column_means(
            {"indptr": np.ascontiguousarray(counts.tocsc().indptr, dtype=np.int32),
             "indices": np.ascontiguousarray(counts.tocsc().indices, dtype=np.int32),
             "data": np.asfortranarray(counts.tocsc().data, dtype=np.float64),
             "shape": (int(counts.shape[0]), int(counts.shape[1]))},
            fit.type_code, n_types, False, n_threads))
    kernel_cov = np.asfortranarray(
        _core.group_covariances(np.asfortranarray(fit.kernel), fit.type_code,
                                n_types, n_threads))

    # The responder block, when the fit has a condition and its interaction
    # terms are present. R looks for terms named "<resp_term>:<neighbour>".
    responder_rows = np.asfortranarray(np.zeros((0, 0), dtype=np.int32))
    responder_keep = np.ascontiguousarray([], dtype=np.int32)
    responder_cov = np.asfortranarray(np.zeros(0))
    n_responder = 0
    if fit.params.get("condition_name") is not None:
        if resp_term is None:
            resp_term = f"{fit.params['condition_name']}{fit.params['condition_levels'][-1]}"
        wanted = [f"{resp_term}:{neighbour}" for neighbour in types]
        keep = [i for i, term in enumerate(wanted) if term in terms]
        if keep:
            n_responder = len(keep)
            # one-based, and into the KERNEL columns, which is what the core
            # uses them to index
            responder_keep = np.ascontiguousarray([i + 1 for i in keep], dtype=np.int32)
            responder_rows = np.asfortranarray(
                [[row_of(wanted[i], f) for f in focal_group] for i in keep], dtype=np.int32)
            level = resp_term[len(fit.params["condition_name"]):]
            dummy = (fit.condition_code
                     == fit.params["condition_levels"].index(level)).astype(np.float64)
            responder_cov = np.asfortranarray(
                _core.group_covariances(
                    np.asfortranarray(fit.kernel * dummy[:, None]), fit.type_code,
                    n_types, n_threads))

    empty = np.asfortranarray(np.zeros((0, 0)))
    blocks = _core.variance_decomposition(
        ct_means=ct_means, group_size=group_size, focal_group=focal_group,
        n_focal=group_size, u=np.asfortranarray(fit.U), se_u=np.asfortranarray(fit.se_U),
        slope_rows=slope_rows, n_kernel=n_types, kernel_cov=kernel_cov,
        responder_rows=responder_rows, n_responder=n_responder,
        responder_keep=responder_keep, responder_cov=responder_cov,
        intercept_rows=intercept_rows,
        toff_var=np.asfortranarray(stats.toff_var) if stats.has_contamination else empty,
        spill_cov=np.asfortranarray(np.zeros(0)), beta_spill=empty, n_spill=0,
        mu_mean=np.asfortranarray(stats.mu_mean), alpha=np.asfortranarray(fit.alpha),
        nb1=nb1, n_threads=n_threads)

    # The core writes focal-major with gene fastest, and R drops the rows the
    # core flagged: a focal with fewer than five cells, or a gene whose blocks
    # do not sum to something finite and positive.
    keep_rows = np.asarray(blocks["keep"]) == 1
    genes = np.tile(np.asarray(fit.genes, dtype=object), n_types)[keep_rows]
    focal = np.repeat(np.asarray(types, dtype=object), n_genes)[keep_rows]
    n_focal_column = np.repeat(group_size, n_genes)[keep_rows]

    def column(name):
        return np.asarray(blocks[name])[keep_rows]

    pct = np.asfortranarray(np.column_stack([
        column("pct_celltype"), column("pct_state"), column("pct_responder"),
        column("pct_spill"), column("pct_residual")]))
    components = np.asfortranarray(np.column_stack([
        column("celltype_offset_sq"), column("v_state_baseline"),
        column("v_state_responder"), column("v_spill"), column("v_disp")]))

    # The aggregates group by focal in the order the surviving rows present
    # them, which is R's group_by() order.
    present = [t for t in types if np.any(focal == t)]
    focal_code = np.ascontiguousarray(
        [present.index(t) for t in focal], dtype=np.int32)
    spec = column("spec")
    aggregates = _core.decomposition_aggregates(
        focal_code=focal_code, n_focals=len(present), pct=pct, components=components,
        spec=np.asfortranarray(spec),
        is_contaminated=np.ascontiguousarray(column("is_contaminated"), dtype=np.int32))

    return Decomposition(
        gene=genes, focal=focal, n_focal=n_focal_column,
        v_state_baseline=column("v_state_baseline"),
        v_state_responder=column("v_state_responder"),
        v_spill=column("v_spill"), v_disp=column("v_disp"),
        celltype_offset_sq=column("celltype_offset_sq"),
        focal_mean=column("focal_mean"), max_other_mean=column("max_other_mean"),
        spec=spec, focal_other_ratio=column("focal_other_ratio"),
        is_contaminated=column("is_contaminated") == 1,
        pct=pct, focal_levels=present,
        mean_table=np.asarray(aggregates["mean"]), specw_table=np.asarray(aggregates["specw"]),
        pooled_table=np.asarray(aggregates["pooled"]),
        total_ss=np.asarray(aggregates["total_ss"]), n_rows=np.asarray(aggregates["n_rows"]),
        n_specific=np.asarray(aggregates["n_specific"]))


# ---------------------------------------------------------------------------
# The observed single frame
# ---------------------------------------------------------------------------


@dataclass
class SingleFrameDecomposition:
    """The observed single-frame decomposition, per (focal, gene).

    The magnitudes here are observed rather than fitted: the cell-type share is
    the between-type sum of squares of log1p(CP10k), and only the split of what
    remains comes from the model. That is why this is the table R reports from
    `varianceDecomposition()`, and why its percentages differ from the
    five-block ones.
    """

    gene: np.ndarray
    focal: np.ndarray
    n_focal: np.ndarray
    pct_celltype: np.ndarray
    pct_spatial: np.ndarray
    pct_responder: np.ndarray | None
    pct_spill: np.ndarray
    pct_residual: np.ndarray
    ss_lineage: np.ndarray
    ss_within: np.ndarray
    denom: np.ndarray
    focal_levels: list[str]
    focal_table: np.ndarray     # focals x blocks, sum-of-squares weighted

    def __len__(self) -> int:
        return self.gene.size

    def __repr__(self) -> str:
        return (f"SingleFrameDecomposition({len(self)} gene-focal rows, "
                f"{len(self.focal_levels)} focal types)")

    @property
    def block_names(self) -> list[str]:
        names = ["Cell type %", "Spatial %"]
        if self.pct_responder is not None:
            names.append("Responder spatial %")
        return names + ["Spillover %", "Residual %"]

    def to_dataframe(self):
        import pandas as pd
        frame = pd.DataFrame({"focal": self.focal, "gene": self.gene,
                              "Cell type %": self.pct_celltype,
                              "Spatial %": self.pct_spatial})
        if self.pct_responder is not None:
            frame["Responder spatial %"] = self.pct_responder
        frame["Spillover %"] = self.pct_spill
        frame["Residual %"] = self.pct_residual
        frame["SS_lineage"] = self.ss_lineage
        frame["SS_within"] = self.ss_within
        frame["denom"] = self.denom
        frame["n_focal"] = self.n_focal
        return frame

    def by_focal(self):
        import pandas as pd
        frame = pd.DataFrame(self.focal_table, columns=self.block_names)
        frame.insert(0, "focal", self.focal_levels)
        return frame


def single_frame_decomposition(fit, counts, decomposition=None, *, n_threads=1):
    """The observed single-frame decomposition.

    `decomposition` is a `pace_decompose()` result, which supplies the
    within-cell-type split; it is computed if not given.
    """
    if decomposition is None:
        decomposition = pace_decompose(fit, counts, n_threads=n_threads)

    types = list(fit.cell_types)
    n_types = len(types)
    genes = np.asarray(fit.genes, dtype=object)
    n_genes = genes.size
    group_size = np.bincount(fit.type_code, minlength=n_types).astype(np.int32)

    csc = counts.tocsc()
    frame = _core.single_frame_statistics(
        {"indptr": np.ascontiguousarray(csc.indptr, dtype=np.int32),
         "indices": np.ascontiguousarray(csc.indices, dtype=np.int32),
         "data": np.asfortranarray(csc.data, dtype=np.float64),
         "shape": (int(csc.shape[0]), int(csc.shape[1]))},
        n_count=np.asfortranarray(np.exp(fit.offset)), group=fit.type_code,
        n_groups=n_types, n_threads=n_threads)

    # The fit's within-type components, as n_types x n_genes with NaN where the
    # decomposition dropped a row. Note the layout: these inputs are indexed
    # (type, gene) with the TYPE varying fastest, where the decomposition's own
    # output runs the other way round.
    has_responder = bool(np.any(decomposition.v_state_responder != 0))
    position = {(f, g): i for i, (f, g) in
                enumerate(zip(decomposition.focal, decomposition.gene))}

    def spread(values):
        out = np.full((n_types, n_genes), np.nan, order="F")
        for t, focal in enumerate(types):
            for j, gene in enumerate(genes):
                row = position.get((focal, gene))
                if row is not None:
                    out[t, j] = values[row]
        return out

    empty = np.asfortranarray(np.zeros((0, 0)))
    shares = _core.single_frame_shares(
        focal_mean=np.asfortranarray(frame["focal_mean"]),
        within_ss=np.asfortranarray(frame["within_ss"]),
        global_mean=np.asfortranarray(frame["global_mean"]),
        group_size=group_size,
        v_state=spread(decomposition.v_state_baseline),
        v_responder=spread(decomposition.v_state_responder) if has_responder else empty,
        v_spill=spread(decomposition.v_spill), v_disp=spread(decomposition.v_disp),
        n_threads=n_threads)

    # The outputs are grouped by type, so each kept type owns a contiguous run
    # of n_genes rows.
    kept = np.flatnonzero(np.asarray(shares["keep"]) == 1)
    if kept.size == 0:
        raise ValueError("no cell type has enough cells for the single-frame decomposition")
    take = np.concatenate([np.arange(k * n_genes, (k + 1) * n_genes) for k in kept])

    def column(name):
        values = np.asarray(shares[name])
        return values[take] if values.size else None

    present = [types[k] for k in kept]
    pct_responder = column("pct_responder") if has_responder else None
    pct_columns = [column("pct_celltype"), column("pct_spatial")]
    if pct_responder is not None:
        pct_columns.append(pct_responder)
    pct_columns += [column("pct_spill"), column("pct_residual")]
    pct = np.asfortranarray(np.column_stack(pct_columns))
    denom = column("denom")

    focal_code = np.ascontiguousarray(np.repeat(np.arange(kept.size), n_genes), dtype=np.int32)
    focal_table = np.asarray(_core.single_frame_focal_blocks(
        focal_code=focal_code, n_focals=kept.size, pct=pct,
        n_blocks=pct.shape[1], denom=np.asfortranarray(denom)))

    return SingleFrameDecomposition(
        gene=np.tile(genes, kept.size), focal=np.repeat(np.asarray(present, dtype=object), n_genes),
        n_focal=np.repeat(group_size[kept], n_genes),
        pct_celltype=column("pct_celltype"), pct_spatial=column("pct_spatial"),
        pct_responder=pct_responder, pct_spill=column("pct_spill"),
        pct_residual=column("pct_residual"), ss_lineage=column("ss_lineage"),
        ss_within=column("ss_within"), denom=denom,
        focal_levels=present, focal_table=focal_table)


# ---------------------------------------------------------------------------
# Pairwise spatial variance
# ---------------------------------------------------------------------------

MIN_CELLS_FOR_PAIRS = 50     # R skips a focal with fewer, as too few to attribute


def pair_variance(fit, counts, single_frame=None, *, block="spatial", resp_term=None,
                  n_threads=1):
    """For each focal type, how much of its variance each neighbour contributes.

    This is the pair heatmap's table, and its definition is particular: the
    focal's spatial share comes from the OBSERVED single frame, and that one
    number is then split across the focal's neighbours by a Pratt attribution
    renormalised over the off-diagonal pairs alone. So a row sums to the
    focal's observed spatial share, and the own-type pair is excluded rather
    than shown as a large diagonal.

    Parameters
    ----------
    block
        `"spatial"` for the baseline neighbour effect, or `"responder"` for its
        interaction with a condition.

    Returns
    -------
    A dict with `focal`, `neighbour` and `val` arrays, off-diagonal pairs only.
    """
    if block not in ("spatial", "responder"):
        raise ValueError('block must be "spatial" or "responder"')
    if single_frame is None:
        single_frame = single_frame_decomposition(fit, counts, n_threads=n_threads)
    if block == "responder" and single_frame.pct_responder is None:
        raise ValueError('this fit has no responder block; use block="spatial"')

    types = list(fit.cell_types)
    n_types = len(types)
    slab = n_types * n_types
    column = single_frame.block_names.index(
        "Spatial %" if block == "spatial" else "Responder spatial %")
    share = dict(zip(single_frame.focal_levels, single_frame.focal_table[:, column]))

    # Sigma is the covariance of the focal's own kernel columns over its own
    # cells, which is that focal's slice; for the responder block it is the
    # same with the condition indicator as a row scale.
    kernel = fit.kernel
    if block == "responder":
        if resp_term is None:
            resp_term = f"{fit.params['condition_name']}{fit.params['condition_levels'][-1]}"
        level = resp_term[len(fit.params["condition_name"]):]
        dummy = (fit.condition_code
                 == fit.params["condition_levels"].index(level)).astype(np.float64)
        kernel = kernel * dummy[:, None]
    covariance = np.asarray(_core.group_covariances(
        np.asfortranarray(kernel), fit.type_code, n_types, n_threads))

    row_of = {name: i for i, name in enumerate(fit.row_names)}
    sizes = np.bincount(fit.type_code, minlength=n_types)

    focal_out, neighbour_out, value_out = [], [], []
    for f, focal in enumerate(types):
        if sizes[f] < MIN_CELLS_FOR_PAIRS or focal not in share:
            continue
        wanted = [(neighbour,
                   f"{focal}::{neighbour}" if block == "spatial"
                   else f"{focal}::{resp_term}:{neighbour}") for neighbour in types]
        keep = [(neighbour, row_of[name]) for neighbour, name in wanted if name in row_of]
        if not keep:
            continue
        neighbours = [n for n, _ in keep]
        rows = [r for _, r in keep]
        columns = [types.index(n) for n in neighbours]
        sigma = covariance[f * slab:(f + 1) * slab].reshape(n_types, n_types, order="F")
        pratt = _core.pair_variance_pratt(
            np.asfortranarray(sigma[np.ix_(columns, columns)]),
            np.asfortranarray(fit.U[rows, :]))

        # off-diagonal only, then renormalised over what is left
        v_pair = np.asarray(pratt["v_pair"])
        off = [i for i, neighbour in enumerate(neighbours) if neighbour != focal]
        total = v_pair[off].sum()
        if not off or total == 0:
            continue
        for i in off:
            focal_out.append(focal)
            neighbour_out.append(neighbours[i])
            value_out.append(share[focal] * v_pair[i] / total)

    return {"focal": np.asarray(focal_out, dtype=object),
            "neighbour": np.asarray(neighbour_out, dtype=object),
            "val": np.asarray(value_out)}
