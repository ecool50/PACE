"""Per-pair driver tables: which genes mediate one focal-neighbour relationship.

For a pair, the genes called at lfsr < 0.05 are ranked by

    driver_score = b^2 * spec^2 * focal_mean

with `b` the shrunken slope, `spec` the gene's specificity to the focal type
and `focal_mean` its mean expression there. Squaring the slope makes the score
direction-blind, so a gene pushed down near the neighbour ranks alongside one
pushed up; multiplying by specificity and level is what stops a gene that is
barely expressed in the focal type from topping the table on a noisy slope.

The ranking is filtered on lfsr alone. There is no specificity cut at this
layer, deliberately: `spec` already enters the score.

Mirrors `pace_top_drivers()` in R.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import _core

LFSR_THRESHOLD = 0.05


@dataclass
class PairDrivers:
    """One focal-neighbour pair's ranked genes."""

    focal: str
    neighbour: str
    term: str
    gene: np.ndarray
    driver_score: np.ndarray
    driver_score_spec4: np.ndarray
    b_clean: np.ndarray
    spec: np.ndarray
    focal_mean: np.ndarray
    mu_bar: np.ndarray
    alpha: np.ndarray
    lfsr: np.ndarray
    sd_shrunk: np.ndarray
    v_s: np.ndarray
    v_resid: np.ndarray
    v_total: np.ndarray
    r2_s: np.ndarray
    u_raw: np.ndarray | None
    v_rxs: np.ndarray | None
    r2_rxs: np.ndarray | None
    status: str
    expected_false_sign: float
    false_sign_rate: float

    def __len__(self) -> int:
        return self.gene.size

    def __repr__(self) -> str:
        return (f"PairDrivers({self.focal} <- {self.neighbour}, {len(self)} genes, "
                f"{self.status})")

    def to_dataframe(self):
        import pandas as pd
        columns = {
            "rank": np.arange(1, len(self) + 1), "gene": self.gene,
            "driver_score": self.driver_score,
            "driver_score_spec4": self.driver_score_spec4, "b_clean": self.b_clean}
        if self.u_raw is not None:
            columns["u_raw"] = self.u_raw
        columns.update({"spec": self.spec, "focal_mean": self.focal_mean, "R2_S": self.r2_s})
        if self.r2_rxs is not None:
            columns["R2_RxS"] = self.r2_rxs
        columns.update({"mu_bar": self.mu_bar, "alpha": self.alpha, "V_S": self.v_s})
        if self.v_rxs is not None:
            columns["V_RxS"] = self.v_rxs
        columns.update({"V_resid": self.v_resid, "V_total": self.v_total,
                        "lfsr": self.lfsr, "sd_shrunk": self.sd_shrunk})
        return pd.DataFrame(columns)


@dataclass
class DriverTables:
    """Every pair's table, keyed "<focal>_<neighbour>" as R keys them."""

    pairs: dict[str, PairDrivers]

    def __len__(self) -> int:
        return len(self.pairs)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            key = f"{key[0]}_{key[1]}"
        return self.pairs[key]

    def __repr__(self) -> str:
        significant = sum(p.status == "significant" for p in self.pairs.values())
        return f"DriverTables({len(self)} pairs, {significant} with three or more calls)"

    def to_dataframe(self):
        """Every pair stacked, with focal and neighbour columns."""
        import pandas as pd
        frames = []
        for pair in self.pairs.values():
            if not len(pair):
                continue
            frame = pair.to_dataframe()
            frame.insert(0, "neighbour", pair.neighbour)
            frame.insert(0, "focal", pair.focal)
            frames.append(frame)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)


def _expected_false_sign(lfsr):
    """How many of these calls should have the wrong sign (Stephens 2017)."""
    called = lfsr[np.isfinite(lfsr)]
    return (float(called.sum()),
            float(called.mean()) if called.size else float("nan"))


def pace_drivers(fit, slopes, decomposition, stats, *, pairs=None, resp_term=None,
                 n_threads=1):
    """Rank the genes mediating each focal-neighbour pair.

    Parameters
    ----------
    fit
        A `PaceFit`.
    slopes
        A `pace_shrink()` result.
    decomposition
        A `pace_decompose()` result, for `spec` and `focal_mean`.
    stats
        A `fit_statistics()` result, for the mean fitted means.
    pairs
        Which (focal, neighbour) pairs to report. None does every ordered pair
        of distinct types, as R does; a subset is a reporting convenience.
    resp_term
        For a condition fit, which condition term's interaction to rank. None
        takes the condition's last level, as R does.

    Returns
    -------
    DriverTables
    """
    types = list(fit.cell_types)
    n_types = len(types)
    if pairs is None:
        pairs = [(f, n) for f in types for n in types if f != n]

    has_responder = fit.params.get("condition_name") is not None
    if has_responder and resp_term is None:
        resp_term = f"{fit.params['condition_name']}{fit.params['condition_levels'][-1]}"

    # var(N_t) among a focal's own cells, which is the diagonal of that focal's
    # kernel covariance; and with the condition indicator as a row scale, the
    # same for the responder product. R reads these out of Z per focal, which
    # comes to the same numbers from the columns the kernel already gives.
    slab = n_types * n_types
    kernel_cov = np.asarray(_core.group_covariances(
        np.asfortranarray(fit.kernel), fit.type_code, n_types, n_threads))
    var_n = np.array([np.diag(kernel_cov[f * slab:(f + 1) * slab].reshape(
        n_types, n_types, order="F")) for f in range(n_types)])
    if has_responder:
        level = resp_term[len(fit.params["condition_name"]):]
        dummy = (fit.condition_code
                 == fit.params["condition_levels"].index(level)).astype(np.float64)
        responder_cov = np.asarray(_core.group_covariances(
            np.asfortranarray(fit.kernel * dummy[:, None]), fit.type_code, n_types, n_threads))
        var_rn = np.array([np.diag(responder_cov[f * slab:(f + 1) * slab].reshape(
            n_types, n_types, order="F")) for f in range(n_types)])
    else:
        var_rn = np.zeros((n_types, n_types))

    # The slopes once by (focal, neighbour, term), and the decomposition's spec
    # and focal_mean once by focal: R splits both for the same reason, that
    # rescanning per pair is hundreds of passes over the same table.
    slope_index = {}
    for i, key in enumerate(zip(slopes.focal, slopes.neighbour, slopes.term)):
        slope_index.setdefault(key, []).append(i)
    block_index = {}
    for i, focal in enumerate(decomposition.focal):
        block_index.setdefault(focal, []).append(i)

    gene_row = {gene: i for i, gene in enumerate(fit.genes)}
    alpha = np.maximum(fit.alpha, 0.0)
    row_of_term = {name: i for i, name in enumerate(fit.row_names)}

    tables = {}
    for focal, neighbour in pairs:
        term = neighbour if not has_responder else f"{resp_term}:{neighbour}"
        rows = slope_index.get((focal, neighbour, term), [])
        blocks = block_index.get(focal, [])
        if not rows or not blocks:
            continue

        # distinct by gene, first occurrence, on both sides, then an inner join
        pair_genes, keep_slope = np.unique(slopes.gene[rows], return_index=True)
        order = np.argsort(keep_slope, kind="stable")
        slope_rows = np.asarray(rows)[keep_slope[order]]
        block_genes = decomposition.gene[blocks]
        _, keep_block = np.unique(block_genes, return_index=True)
        block_rows = np.asarray(blocks)[np.sort(keep_block)]
        block_of_gene = {decomposition.gene[r]: r for r in block_rows}

        matched = [(s, block_of_gene[slopes.gene[s]]) for s in slope_rows
                   if slopes.gene[s] in block_of_gene]
        if not matched:
            continue
        slope_rows = np.array([s for s, _ in matched])
        block_rows = np.array([b for _, b in matched])
        genes = slopes.gene[slope_rows]
        columns = np.array([gene_row[g] for g in genes])

        focal_index = types.index(focal)
        neighbour_index = types.index(neighbour)
        u_raw = (np.asfortranarray(
            fit.U[row_of_term[f"{focal}::{neighbour}"], columns]) if has_responder
            else np.asfortranarray(np.zeros(0)))
        scores = _core.driver_scores(
            estimate_shrunk=np.asfortranarray(slopes.estimate_shrunk[slope_rows]),
            spec=np.asfortranarray(decomposition.spec[block_rows]),
            focal_mean=np.asfortranarray(decomposition.focal_mean[block_rows]),
            mu_bar=np.asfortranarray(stats.mu_mean[focal_index, columns]),
            alpha=np.asfortranarray(alpha[columns]), u_raw=u_raw,
            var_n=float(var_n[focal_index, neighbour_index]),
            var_rn=float(var_rn[focal_index, neighbour_index]),
            has_responder=has_responder)

        # lfsr alone, then by driver score descending; the sort is stable, as
        # R's arrange() is, so ties keep the order they came in.
        lfsr = slopes.lfsr[slope_rows]
        called = np.flatnonzero(lfsr < LFSR_THRESHOLD)
        ranked = called[np.argsort(-np.asarray(scores["driver_score"])[called], kind="stable")]
        expected, rate = _expected_false_sign(lfsr[ranked])

        def pick(values):
            return np.asarray(values)[ranked]

        tables[f"{focal}_{neighbour}"] = PairDrivers(
            focal=focal, neighbour=neighbour, term=term, gene=genes[ranked],
            driver_score=pick(scores["driver_score"]),
            driver_score_spec4=pick(scores["driver_score_spec4"]),
            b_clean=pick(slopes.estimate_shrunk[slope_rows]),
            spec=pick(decomposition.spec[block_rows]),
            focal_mean=pick(decomposition.focal_mean[block_rows]),
            mu_bar=pick(stats.mu_mean[focal_index, columns]), alpha=pick(alpha[columns]),
            lfsr=lfsr[ranked], sd_shrunk=pick(slopes.sd_shrunk[slope_rows]),
            v_s=pick(scores["V_S"]), v_resid=pick(scores["V_resid"]),
            v_total=pick(scores["V_total"]), r2_s=pick(scores["R2_S"]),
            u_raw=pick(u_raw) if has_responder else None,
            v_rxs=pick(scores["V_RxS"]) if has_responder else None,
            r2_rxs=pick(scores["R2_RxS"]) if has_responder else None,
            status="significant" if ranked.size >= 3 else "honestly null",
            expected_false_sign=expected, false_sign_rate=rate)

    return DriverTables(pairs=tables)
