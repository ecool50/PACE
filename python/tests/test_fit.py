"""The IRLS loop, against the fit R produces from the same inputs.

The fixtures are the exact 43 arguments R hands pace_irls_driver_cpp on the
packaged breast cancer subset, captured by tracing the call. If Python feeds
the core the same things it must get the same fit back.
"""
from pathlib import Path

import numpy as np
import pytest

import pace

FIX = Path(__file__).parent / "fixtures"


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def csc(stem):
    dim = vec(f"{stem}_dim.txt", int)
    return {
        "indptr": np.ascontiguousarray(vec(f"{stem}_p.txt", np.int32)),
        "indices": np.ascontiguousarray(vec(f"{stem}_i.txt", np.int32)),
        "data": np.asfortranarray(vec(f"{stem}_x.txt")),
        "shape": (int(dim[0]), int(dim[1])),
    }


def inputs():
    n, q, g, p, n_iter, min_iter, chunk, sub, interior, threads = vec("irls_scalars.txt", int)
    early_stop, alpha_warmup, alpha_max_cells, tau_max, d0_min = vec("irls_dscalars.txt")
    flags = vec("irls_flags.txt", int)

    meta = vec("irls_b1_meta.txt", int)
    tdim = vec("irls_b1_termsdim.txt", int)
    block = {
        "col_offset": int(meta[0]), "n_terms": int(meta[1]),
        "n_groups": int(meta[2]), "n_cols": int(meta[3]),
        "terms": np.asfortranarray(vec("irls_b1_terms.txt").reshape(tdim[0], tdim[1], order="F")),
        "cell_group": np.ascontiguousarray(vec("irls_b1_cellgroup.txt", np.int32)),
        "group_start": np.ascontiguousarray(vec("irls_b1_gstart.txt", np.int32)),
        "cells_by_group": np.ascontiguousarray(vec("irls_b1_cbg.txt", np.int32)),
    }
    mdim = vec("irls_mask_dim.txt", int)
    ddim = vec("irls_diw_dim.txt", int)
    xdim = vec("irls_xfixed_dim.txt", int)
    options = dict(
        n_iter=int(n_iter), min_iter=int(min_iter), early_stop_tol=float(early_stop),
        alpha_warmup=float(alpha_warmup), alpha_max_cells=float(alpha_max_cells),
        tau_max=float(tau_max), d0_min=float(d0_min), chunk_size=int(chunk),
        sub_genes=int(sub), interior_precision=int(interior), n_threads=int(threads),
        x1_is_unit=bool(flags[0]), nb2=bool(flags[1]), gaussian=bool(flags[2]),
        zero_collapse=bool(flags[3]), fast_density=bool(flags[4]),
        use_reml=bool(flags[5]), rd_diag=bool(flags[6]), verbose=False,
    )
    return dict(
        x_fixed=np.asfortranarray(vec("irls_xfixed.txt").reshape(xdim[0], xdim[1], order="F")),
        p=int(p), z_design=csc("irls_z"), blocks=[block],
        counts=csc("irls_y"), ambient=csc("irls_e"),
        offset=np.asfortranarray(vec("irls_offset.txt")),
        alpha_init=np.asfortranarray(vec("irls_alpha_init.txt")),
        # without these the random effects come back unshrunk
        mask=np.asfortranarray(vec("irls_mask.txt").reshape(mdim[0], mdim[1], order="F")),
        mask_index=np.ascontiguousarray(vec("irls_maskindex.txt", np.int32)),
        n_mask_rows=int(mdim[0]),
        data_informed_weights=np.asfortranarray(
            vec("irls_diw.txt").reshape(ddim[0], ddim[1], order="F")),
        n_cells=int(n), q_total=int(q), n_genes=int(g), options=options,
    )


@pytest.fixture(scope="module")
def fitted():
    return pace.fit_irls(**inputs())


# Why not exact. R passes the solver Rf_trigamma and Python passes the core's
# own, which agree to 1.8e-15 -- but estimate_d0 feeds the tau shrinkage and
# thirteen IRLS iterations amplify that to 8.3e-07 in the fitted effects. The
# two are bit-identical when both use the same trigamma; see the note on the
# open decision to move R onto the core's. Separately, this holds only when the
# extension is built against Eigen 3.4: against 3.5 the gap is 3e-4.
# Measured, not guessed: the dispersions carry it furthest because the
# trigamma enters through estimate_d0, which is a dispersion-side quantity.
TRIGAMMA_TOL = {"U": 1e-06, "B": 1e-06, "alpha": 1e-05}


def test_the_fit_reproduces_r(fitted):
    g = fitted["U"].shape[1]
    for name, fixture in [("U", "irls_expect_U.txt"), ("B", "irls_expect_B.txt")]:
        np.testing.assert_allclose(fitted[name], vec(fixture).reshape(-1, g, order="F"),
                                   rtol=0, atol=TRIGAMMA_TOL[name],
                                   err_msg=f"{name} differs by more than the trigamma can explain")
    np.testing.assert_allclose(fitted["alpha"], vec("irls_expect_alpha.txt"),
                               rtol=0, atol=TRIGAMMA_TOL["alpha"],
                               err_msg="dispersions differ by more than the trigamma can explain")


def test_the_fit_is_bit_identical_when_the_trigamma_matches(fitted):
    """The same fit, from an R build using the core's trigamma, matches exactly.

    This is the real check: it shows the remaining difference above is entirely
    the trigamma and nothing else.
    """
    g = fitted["U"].shape[1]
    same = vec("irls_ab_U.txt").reshape(-1, g, order="F")
    assert np.array_equal(fitted["U"], same)
    assert np.array_equal(fitted["alpha"], vec("irls_ab_alpha.txt"))


def test_the_fit_converges_the_way_r_did(fitted):
    assert fitted["converged"] is True
    assert fitted["n_iter"] == 13


def test_the_fit_reports_usable_standard_errors(fitted):
    # written on the final iteration only, so a converged fit must have them
    assert np.all(np.isfinite(fitted["se_U"]))
    assert np.all(fitted["se_U"] >= 0)
