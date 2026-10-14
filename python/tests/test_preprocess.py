"""The stages between the kernels and the fit, against the R package.

These are small functions, but each one changes what the solver sees, and a
silent difference here would show up only as a slightly different fit.
"""
from pathlib import Path

import numpy as np

import pace

FIX = Path(__file__).parent / "fixtures"


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def kernel():
    x = np.ascontiguousarray(vec("kern_x.txt"))
    y = np.ascontiguousarray(vec("kern_y.txt"))
    cell_type = vec("kern_type.txt").astype(np.int32)
    n_types = int(vec("kern_ntypes.txt", int)[0])
    group = np.zeros(x.size, dtype=np.int32)
    out = pace.neighbour_kernels(x, y, cell_type, n_types, group, eps=90.0)
    return out["k_bio"], cell_type, group, n_types


def test_dropping_sparse_pairs_matches_r():
    k_bio, cell_type, _, n_types = kernel()
    out = pace.drop_sparse_kernel(k_bio, cell_type, 30.0)
    expected = vec("pre_drop.txt").reshape(k_bio.shape, order="F")
    np.testing.assert_allclose(out["kernel"], expected, rtol=1e-14, atol=0)
    assert out["n_dropped"] == int(vec("pre_drop_n.txt", int)[0])
    # dropping means zeroing, so it can only remove signal
    assert np.all((out["kernel"] == 0) | (out["kernel"] == k_bio))


def test_centring_within_groups_matches_r():
    k_bio, cell_type, group, n_types = kernel()
    dropped = pace.drop_sparse_kernel(k_bio, cell_type, 30.0)["kernel"]
    out = pace.centre_within_groups(dropped, group, cell_type, 1, n_types)
    np.testing.assert_allclose(out, vec("pre_centre.txt").reshape(out.shape, order="F"),
                               rtol=1e-12, atol=1e-12)
    # centring leaves each (image, type) group with mean zero
    for t in range(n_types):
        rows = cell_type == t
        if rows.sum() > 1:
            np.testing.assert_allclose(out[rows].mean(axis=0), 0, atol=1e-9)


def test_standardise_matches_r():
    k_bio = kernel()[0]
    out = pace.standardise(np.asfortranarray(k_bio[:, 0].copy()))
    np.testing.assert_allclose(out["values"], vec("pre_std.txt"), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(out["sd"], vec("pre_std_sd.txt")[0], rtol=1e-12)


def test_anchor_mask_matches_r():
    dim = vec("pre_det_dim.txt", int)
    det = np.asfortranarray(vec("pre_det.txt").reshape(dim[0], dim[1], order="F"))
    out = pace.anchor_mask(det, 0.2, 0.1)
    np.testing.assert_allclose(out["mask"], vec("pre_anchor.txt").reshape(det.shape, order="F"),
                               rtol=0, atol=0)
    np.testing.assert_array_equal(out["n_anchor"], vec("pre_anchor_n.txt", np.int32))
    # a mask is an indicator, and a type never anchors a gene it owns
    assert set(np.unique(out["mask"])) <= {0.0, 1.0}
    assert out["mask"].sum() == out["n_anchor"].sum()


def test_data_informed_weights_match_r():
    dim = vec("pre_det_dim.txt", int)
    det = np.asfortranarray(vec("pre_det.txt").reshape(dim[0], dim[1], order="F"))
    focal_of_row = np.ascontiguousarray(vec("pre_focalrow.txt", np.int32))
    scale_of_row = np.asfortranarray(vec("pre_scalerow.txt"))
    out = pace.data_informed_weights(det, focal_of_row, scale_of_row)
    ddim = vec("pre_diw_dim.txt", int)
    np.testing.assert_allclose(out, vec("pre_diw.txt").reshape(ddim[0], ddim[1], order="F"),
                               rtol=1e-14, atol=0)
    # weights are a prior scale: positive, and capped at one
    assert out.min() > 0 and out.max() <= 1.0
