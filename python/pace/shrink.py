"""Multivariate adaptive shrinkage of the neighbour slopes.

The fit gives one slope per (focal type, neighbour type, gene) with a standard
error. Most are noise. mash pools them: for one neighbour, the focal types are
treated as conditions measured on the same genes, so a slope that is consistent
across focals is kept and an isolated one is pulled towards zero, and every
estimate comes back with a local false sign rate.

This mirrors `apply_mashr_shrinkage()` in R, which is where the shape of the
problem comes from: one mash fit per neighbour, the focal types as its
conditions, and a single gene set shared by every slice so that one lfsr cut
means the same thing across them.

The mixture itself is the compiled core's, which reproduces mashr's canonical
covariances exactly. It does not implement the data-driven ones (`cov_pca` and
`cov_ed`), so this is R at `data_driven = FALSE`: on the two canonical cohorts
that costs 11 of 1,637 breast cancer calls and none of 46 melanoma ones.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import _core

# Both from R. An SE above the cap is a BLUP that never converged; below the
# floor, estimate and error are both effectively zero and mash's likelihood
# matrix goes singular on them.
SE_CAP = 50.0
SE_FLOOR = 1e-4


@dataclass
class ShrunkenSlopes:
    """One row per (gene, focal, neighbour) that survived shrinkage."""

    gene: np.ndarray
    focal: np.ndarray
    neighbour: np.ndarray
    term: np.ndarray
    estimate: np.ndarray
    std_error: np.ndarray
    estimate_shrunk: np.ndarray
    sd_shrunk: np.ndarray
    lfsr: np.ndarray

    def __len__(self) -> int:
        return self.gene.size

    def __repr__(self) -> str:
        calls = int(np.sum(self.lfsr < 0.05)) if len(self) else 0
        return (f"ShrunkenSlopes({len(self)} gene-focal-neighbour rows, "
                f"{len(set(self.neighbour.tolist()))} neighbours, "
                f"{calls} calls at lfsr < 0.05)")

    def to_dataframe(self):
        """The same thing as a pandas DataFrame, with R's column names.

        pandas is imported here rather than at the top so it stays optional.
        """
        import pandas as pd
        return pd.DataFrame({
            "gene": self.gene, "focal": self.focal, "neighbour": self.neighbour,
            "term": self.term, "estimate": self.estimate, "std.error": self.std_error,
            "estimate_shrunk": self.estimate_shrunk, "sd_shrunk": self.sd_shrunk,
            "lfsr": self.lfsr,
        })

    def calls(self, threshold=0.05):
        """How many effects are called, and how many of those should have the
        wrong sign.

        lfsr is the posterior probability that a called effect's sign is wrong,
        so summing it over the calls estimates how many of them are sign errors
        (Stephens 2017). Both numbers are only as calibrated as the lfsr.
        """
        called = self.lfsr[np.isfinite(self.lfsr) & (self.lfsr < threshold)]
        return {"n_calls": called.size,
                "expected_false_sign": float(called.sum()),
                "false_sign_rate": float(called.mean()) if called.size else float("nan")}


def _default_resp_term(fit):
    """R's default: the term of the condition's LAST level.

    `.pace_prepare()` builds a term for every level, the reference included,
    and takes the last. The reference is the first level, so the last is always
    one of the fitted dummies.
    """
    name = fit.params.get("condition_name")
    if name is None:
        return None
    return f"{name}{fit.params['condition_levels'][-1]}"


def _null_correlation(bhat, shat):
    """mash's V: the correlation between conditions among the null-looking genes.

    The focal types are estimated in the same per-gene solve, so their errors
    are correlated, and mash's default V = I mis-calibrates lfsr. This is
    mashr's `estimate_null_correlation_simple` with its own defaults: the genes
    whose largest absolute z is under 2, then the correlation of their z
    scores. None means the estimate is not usable and independence stands.
    """
    z = bhat / shat
    with np.errstate(invalid="ignore"):
        spread = np.nanstd(np.where(np.isfinite(z), z, np.nan), axis=0, ddof=1)
    live = np.flatnonzero(np.isfinite(spread) & (spread > 0))
    if live.size < 2:
        return None

    z_live = z[:, live]
    nullish = np.max(np.abs(z_live), axis=1) < 2.0
    if int(nullish.sum()) < live.size:        # mashr raises here; R falls back
        return None
    covariance = np.cov(z_live[nullish, :], rowvar=False, ddof=1)
    scale = np.sqrt(np.diag(covariance))
    if not np.all(np.isfinite(scale)) or np.any(scale <= 0):
        return None
    correlation = covariance / np.outer(scale, scale)
    if not np.all(np.isfinite(correlation)):
        return None
    if not np.allclose(correlation, correlation.T, atol=1e-8):
        return None
    if np.min(np.linalg.eigvalsh(correlation)) <= 1e-8:
        return None

    v = np.eye(bhat.shape[1])
    v[np.ix_(live, live)] = correlation
    return np.asfortranarray(v)


def pace_shrink(fit, *, resp_term=None, null_correlation=True, verbose=False):
    """Shrink a fit's neighbour slopes.

    Parameters
    ----------
    fit
        A `PaceFit` from `pace_fit()` or `fit_anndata()`.
    resp_term
        For a condition fit, which condition term's interaction to shrink. The
        slopes shrunk are then `<resp_term>:<neighbour>`, so what comes back is
        the responder layer and NOT the baseline proximity effect. None takes
        R's default, the term of the condition's last level. Ignored when the
        fit has no condition.
    null_correlation
        Estimate mash's V from the null-looking genes rather than assuming the
        focal types are independent. True is R's default. The published breast
        cancer counts were produced with False, so pass False to match them.

    Returns
    -------
    ShrunkenSlopes
    """
    focals = list(fit.cell_types)
    neighbours = list(fit.cell_types)
    if fit.params.get("condition_name") is None:
        resp_term = None
    elif resp_term is None:
        resp_term = _default_resp_term(fit)

    # (focal, term) -> row of U, from the names the fit already carries.
    row_of = {}
    for row, name in enumerate(fit.row_names):
        focal, _, term = name.partition("::")
        row_of[(focal, term)] = row
    genes = np.asarray(fit.genes, dtype=object)

    # Pass one: build each neighbour's slice and note the genes it can use.
    slices = []
    for neighbour in neighbours:
        term = neighbour if resp_term is None else f"{resp_term}:{neighbour}"
        rows = [row_of.get((focal, term)) for focal in focals]
        present = [i for i, row in enumerate(rows) if row is not None]
        if len(present) < 2:
            if verbose:
                print(f"  [mash] term '{term}' has fewer than two focal columns, skipping")
            continue
        columns = [focals[i] for i in present]
        bhat = np.asfortranarray(np.column_stack([fit.U[rows[i], :] for i in present]))
        shat = np.asfortranarray(np.column_stack([fit.se_U[rows[i], :] for i in present]))

        keep = (np.all(np.isfinite(bhat), axis=1) & np.all(np.isfinite(shat), axis=1)
                & np.all(shat > 0, axis=1) & np.all(shat <= SE_CAP, axis=1)
                # every focal's BLUP collapsed to about zero with about no
                # error: a degenerate boundary fit, and mash chokes on it
                & ~np.all((np.abs(bhat) < SE_FLOOR) & (shat < SE_FLOOR), axis=1))
        if int(keep.sum()) < 5:
            if verbose:
                print(f"  [mash] term '{term}' has fewer than five well-fit genes, skipping")
            continue
        slices.append({"neighbour": neighbour, "term": term, "columns": columns,
                       "bhat": bhat, "shat": shat, "usable": set(genes[keep].tolist())})

    if not slices:
        return _empty()

    # Every slice is shrunk over the same genes: mash calibrates its null
    # against the genes it is given, so slices filtered apart would carry lfsr
    # on different bases and one cut across them would not compare.
    common = set.intersection(*[s["usable"] for s in slices])
    if verbose:
        for s in slices:
            dropped = len(s["usable"] - common)
            if dropped:
                print(f"  [mash] term '{s['term']}': {dropped} genes dropped to keep "
                      "the gene set common across neighbours")
    take = np.array([gene in common for gene in genes])
    if int(take.sum()) < 5:
        return _empty()
    kept_genes = genes[take]

    pieces = []
    for s in slices:
        bhat = np.asfortranarray(s["bhat"][take, :])
        shat = np.asfortranarray(np.maximum(s["shat"][take, :], SE_FLOOR))
        v = _null_correlation(bhat, shat) if null_correlation else None
        out = _core.mash(bhat, shat, **({"v": v} if v is not None else {}))
        n_genes, n_focals = bhat.shape
        pieces.append(ShrunkenSlopes(
            # column-major, as R unrolls it: gene varies fastest within focal
            gene=np.tile(kept_genes, n_focals),
            focal=np.repeat(np.asarray(s["columns"], dtype=object), n_genes),
            neighbour=np.full(n_genes * n_focals, s["neighbour"], dtype=object),
            term=np.full(n_genes * n_focals, s["term"], dtype=object),
            estimate=bhat.ravel(order="F"),
            std_error=shat.ravel(order="F"),
            estimate_shrunk=np.asarray(out["posterior_mean"]).ravel(order="F"),
            sd_shrunk=np.asarray(out["posterior_sd"]).ravel(order="F"),
            lfsr=np.asarray(out["lfsr"]).ravel(order="F")))
        if verbose:
            calls = int(np.sum(pieces[-1].lfsr < 0.05))
            print(f"  [mash] {s['term']}: {n_genes} genes shrunk; "
                  f"gene-focal calls (lfsr<0.05) = {calls}")

    return _concatenate(pieces)


def _empty():
    text = np.array([], dtype=object)
    number = np.array([], dtype=float)
    return ShrunkenSlopes(text, text.copy(), text.copy(), text.copy(),
                          number, number.copy(), number.copy(), number.copy(),
                          number.copy())


def _concatenate(pieces):
    fields = ("gene", "focal", "neighbour", "term", "estimate", "std_error",
              "estimate_shrunk", "sd_shrunk", "lfsr")
    return ShrunkenSlopes(**{
        field: np.concatenate([getattr(piece, field) for piece in pieces])
        for field in fields})
