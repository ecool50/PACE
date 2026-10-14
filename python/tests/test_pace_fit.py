"""`pace_fit()` end to end, against the fit R produces from the same raw data.

test_fit.py checks the solver when it is handed R's own captured inputs. This
checks the layer above: starting from nothing but counts, coordinates and cell
type labels, Python builds every design, kernel, offset, anchor mask and prior
weight itself. Agreement here means the orchestration matches R's as well as
the numerics do.

The reference is an R build using the core's trigamma, because that is the one
remaining difference between the two languages. Against R's stock trigamma the
fit still agrees, but only to the tolerance test_fit.py documents.
"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

import pace
from pace.fit import pace_fit

FIX = Path(__file__).parent / "fixtures"


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def lines(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


def _raw_inputs():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    return counts, vec("raw_x.txt"), vec("raw_y.txt"), lines("raw_celltype.txt")


@pytest.fixture(scope="module")
def fit():
    counts, x, y, cell_type = _raw_inputs()
    return pace_fit(counts, x, y, cell_type, n_threads=1)


def test_the_random_effects_are_bit_identical_to_r(fit):
    reference = vec("irls_ab_U.txt").reshape(fit.U.shape, order="F")
    assert np.array_equal(fit.U, reference)


def test_the_dispersions_are_bit_identical_to_r(fit):
    assert np.array_equal(fit.alpha, vec("irls_ab_alpha.txt"))


def test_the_fit_agrees_with_rs_stock_trigamma_to_the_documented_tolerance(fit):
    """R's released build uses Rf_trigamma, the core uses its own.

    Thirteen IRLS iterations amplify a 1.8e-15 difference in the trigamma into
    this, so the tolerance is the one test_fit.py justifies, not a guess.
    """
    np.testing.assert_allclose(fit.U, vec("raw_expect_U.txt").reshape(fit.U.shape, order="F"),
                               rtol=0, atol=1e-06)
    np.testing.assert_allclose(fit.alpha, vec("raw_expect_alpha.txt"), rtol=0, atol=1e-05)


def test_the_random_effect_rows_are_named_the_way_r_names_them(fit):
    assert fit.row_names == lines("raw_expect_rownames.txt")


def test_the_fit_converges_the_way_r_did(fit):
    n_iter, converged = vec("raw_expect_meta.txt", int)
    assert fit.n_iter == int(n_iter)
    assert fit.converged is bool(converged)


def test_the_contamination_loadings_match_rs(fit):
    """One loading per cell, against R's `percell_bleed_rho`.

    R's reference here is its stock build, so this carries the same trigamma
    difference as the fit itself: 1.1e-06 absolute on loadings that run from
    1e-04 to 3.7.
    """
    reference = vec("raw_expect_rho.txt")
    assert fit.rho.size == reference.size
    np.testing.assert_allclose(fit.rho, reference, rtol=1e-05, atol=0)


def test_a_one_level_condition_is_refused(fit):
    counts, x, y, cell_type = _raw_inputs()
    with pytest.raises(ValueError, match="at least two"):
        pace_fit(counts, x, y, cell_type, condition=["PD"] * x.size)


# --- several images -------------------------------------------------------
# The packaged subset is one image, so it cannot tell a global neighbour search
# from a per-image one. These split it on the median x and hold both settings to
# R. The two R fits differ by 0.167 in U, so the pair discriminates: taking the
# wrong branch is not a rounding error. Tolerances are the trigamma floor.

@pytest.fixture(scope="module")
def two_image_inputs():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    return counts, vec("raw_x.txt"), vec("raw_y.txt"), lines("raw_celltype.txt"), \
        lines("two_image_label.txt")


@pytest.mark.parametrize("kernel_per_image,tag", [(False, "global"), (True, "perimage")])
def test_both_neighbour_scopes_match_r_on_two_images(two_image_inputs, kernel_per_image, tag):
    counts, x, y, cell_type, image = two_image_inputs
    fit = pace_fit(counts, x, y, cell_type, image=image,
                   kernel_per_image=kernel_per_image, n_threads=1)
    reference = vec(f"two_image_{tag}_U.txt").reshape(fit.U.shape, order="F")
    np.testing.assert_allclose(fit.U, reference, rtol=0, atol=1e-06)
    np.testing.assert_allclose(fit.alpha, vec(f"two_image_{tag}_alpha.txt"), rtol=0, atol=1e-05)


def test_the_default_neighbour_scope_is_global_like_rs(two_image_inputs):
    """R's `kernel_per_image` defaults to FALSE even with several images."""
    counts, x, y, cell_type, image = two_image_inputs
    fit = pace_fit(counts, x, y, cell_type, image=image, n_threads=1)
    reference = vec("two_image_global_U.txt").reshape(fit.U.shape, order="F")
    np.testing.assert_allclose(fit.U, reference, rtol=0, atol=1e-06)
    assert fit.params["kernel_per_image"] is False


# --- what must be refused rather than guessed ----------------------------

def test_a_missing_cell_type_is_refused_not_turned_into_a_type(two_image_inputs):
    """str(nan) is "nan", which would otherwise be a cell type like any other."""
    counts, x, y, cell_type, _ = two_image_inputs
    labels = list(cell_type)
    labels[7] = np.nan
    with pytest.raises(ValueError, match="every cell needs a cell type"):
        pace_fit(counts, x, y, labels)


def test_an_explicit_type_order_is_honoured(two_image_inputs):
    """R's `types` argument, which is how a locked row order is reproduced."""
    counts, x, y, cell_type, _ = two_image_inputs
    reversed_types = sorted(set(cell_type))[::-1]
    fit = pace_fit(counts, x, y, cell_type, types=reversed_types, n_threads=1)
    assert fit.cell_types == reversed_types
    assert fit.row_names[0] == f"{reversed_types[0]}::(Intercept)"


def test_an_unknown_label_is_refused(two_image_inputs):
    counts, x, y, cell_type, _ = two_image_inputs
    with pytest.raises(ValueError, match="not in `types`"):
        pace_fit(counts, x, y, cell_type, types=sorted(set(cell_type))[1:])


# --- sparse matrices that are not in canonical form ----------------------
# scipy will hand over a CSC holding duplicate or unsorted entries, and
# `tocsc()` on one passes it straight through. The core reads stored entries one
# at a time, so a duplicated (cell, gene) counts twice. Such a matrix can be
# elementwise identical to a good one, hold none but non-negative integers, and
# give the same library sizes, so nothing else here would catch it. R is safe
# only because Matrix refuses the object outright.

def _split_stored_entries(matrix):
    """The same matrix, with every value >= 2 stored as two entries."""
    coo = matrix.tocoo()
    split = coo.data >= 2
    rows = np.concatenate([coo.row, coo.row[split]])
    cols = np.concatenate([coo.col, coo.col[split]])
    data = np.concatenate([np.where(split, coo.data - 1, coo.data), np.ones(split.sum())])
    order = np.lexsort((rows, cols))
    indptr = np.searchsorted(cols[order], np.arange(matrix.shape[1] + 1))
    return csc_matrix((data[order], rows[order], indptr), shape=matrix.shape)


def _shuffle_within_columns(matrix):
    """The same matrix, with each column's row indices out of order."""
    out = matrix.copy()
    rng = np.random.default_rng(0)
    for column in range(out.shape[1]):
        lo, hi = out.indptr[column], out.indptr[column + 1]
        order = rng.permutation(hi - lo)
        out.indices[lo:hi] = out.indices[lo:hi][order]
        out.data[lo:hi] = out.data[lo:hi][order]
    out.has_sorted_indices = False
    out.has_canonical_format = False
    return out


@pytest.mark.parametrize("damage", [_split_stored_entries, _shuffle_within_columns],
                         ids=["duplicate-entries", "unsorted-indices"])
def test_a_non_canonical_matrix_gives_the_same_fit(fit, damage):
    counts, x, y, cell_type = _raw_inputs()
    broken = damage(counts)
    assert not broken.has_canonical_format
    assert np.array_equal(np.asarray(broken.todense()), np.asarray(counts.todense()))
    refit = pace_fit(broken, x, y, cell_type, n_threads=1)
    assert np.array_equal(refit.U, fit.U)
    assert np.array_equal(refit.alpha, fit.alpha)


def test_normalising_does_not_touch_the_callers_matrix(fit):
    counts, x, y, cell_type = _raw_inputs()
    broken = _split_stored_entries(counts)
    before = broken.nnz, broken.data.copy()
    pace_fit(broken, x, y, cell_type, n_threads=1)
    assert broken.nnz == before[0]
    assert np.array_equal(broken.data, before[1])


# --- a condition cohort --------------------------------------------------
# `~ 1 + condition * (types)` over a `~ 1 + condition` fixed design. Two levels
# pins the columns; THREE pins their order, because model.matrix() varies the
# condition level fastest within each type and a two-level fit cannot show it.
# Tolerances are the trigamma floor: these references are R's stock build.

@pytest.mark.parametrize("tag,name,n_terms", [("cond2", "Arm", 18), ("cond3", "Arm3", 27)])
def test_a_condition_fit_matches_r(tag, name, n_terms):
    counts, x, y, cell_type = _raw_inputs()
    fit = pace_fit(counts, x, y, cell_type, condition=lines(f"{tag}_label.txt"),
                   condition_name=name, n_threads=1)
    dim = vec(f"{tag}_U_dim.txt", int)
    assert fit.U.shape == (int(dim[0]), int(dim[1]))
    assert len(fit.params["term_names"]) == n_terms
    np.testing.assert_allclose(
        fit.U, vec(f"{tag}_U.txt").reshape(fit.U.shape, order="F"), rtol=0, atol=1e-06)
    np.testing.assert_allclose(
        fit.B, vec(f"{tag}_B.txt").reshape(fit.B.shape, order="F"), rtol=0, atol=1e-06)
    np.testing.assert_allclose(fit.alpha, vec(f"{tag}_alpha.txt"), rtol=0, atol=1e-05)


@pytest.mark.parametrize("tag,name", [("cond2", "Arm"), ("cond3", "Arm3")])
def test_the_condition_rows_are_named_the_way_r_names_them(tag, name):
    """The interaction order is R's, level innermost, not the obvious one."""
    counts, x, y, cell_type = _raw_inputs()
    fit = pace_fit(counts, x, y, cell_type, condition=lines(f"{tag}_label.txt"),
                   condition_name=name, n_threads=1)
    assert fit.row_names == lines(f"{tag}_rownames.txt")
    assert fit.params["fixed_names"] == lines(f"{tag}_fixednames.txt")


def test_the_condition_terms_take_their_prefix_from_the_name():
    counts, x, y, cell_type = _raw_inputs()
    fit = pace_fit(counts, x, y, cell_type, condition=lines("cond2_label.txt"),
                   condition_name="Responder", n_threads=1)
    assert "Responder::PR" not in fit.row_names
    assert fit.params["fixed_names"] == ["(Intercept)", "ResponderPR"]
    assert any(name.endswith("::ResponderPR:Tumour") for name in fit.row_names)


def test_the_reference_level_is_the_first_and_is_dropped():
    counts, x, y, cell_type = _raw_inputs()
    flipped = pace_fit(counts, x, y, cell_type, condition=lines("cond2_label.txt"),
                       condition_name="Arm", condition_levels=["PR", "PD"], n_threads=1)
    assert flipped.params["fixed_names"] == ["(Intercept)", "ArmPD"]
    # the same model, reparameterised: the intercept absorbs the other arm, so
    # the fit is not the same numbers, only the same shape
    assert flipped.U.shape == (144, 278)


def test_a_missing_condition_is_refused():
    counts, x, y, cell_type = _raw_inputs()
    labels = list(lines("cond2_label.txt"))
    labels[3] = np.nan
    with pytest.raises(ValueError, match="every cell needs a condition"):
        pace_fit(counts, x, y, cell_type, condition=labels, condition_name="Arm")


def test_an_unknown_condition_label_is_refused():
    counts, x, y, cell_type = _raw_inputs()
    with pytest.raises(ValueError, match="not a level"):
        pace_fit(counts, x, y, cell_type, condition=lines("cond2_label.txt"),
                 condition_name="Arm", condition_levels=["PD"])


# --- what the audit of the condition path turned up --------------------
# Each of these produced an all-NaN fit, or a different fit from R's, in
# silence. They are grouped because they share a cause: something that is
# empty, or chosen by the wrong rule, and nothing saying so.

def test_a_cell_type_with_no_cells_does_not_poison_the_fit():
    """An empty type's detection rate is NaN from the core; R fills zero.

    Reachable through the documented route: `types` is how a row order is
    matched to R or shared across cohorts, so a type absent here is normal.
    """
    counts, x, y, cell_type = _raw_inputs()
    fit = pace_fit(counts, x, y, cell_type, types=sorted(set(cell_type)) + ["Ghost"],
                   n_threads=1)
    assert np.all(np.isfinite(fit.U))
    assert np.all(np.isfinite(fit.B))
    assert fit.converged
    ghost = [i for i, name in enumerate(fit.row_names) if name.startswith("Ghost::")]
    assert np.array_equal(fit.U[ghost], np.zeros_like(fit.U[ghost]))


def test_an_empty_cell_type_barely_moves_the_real_rows(fit):
    """It adds all-zero rows, which shift the pooled variance components a little."""
    counts, x, y, cell_type = _raw_inputs()
    padded = pace_fit(counts, x, y, cell_type, types=sorted(set(cell_type)) + ["Ghost"],
                      n_threads=1)
    index = {name: i for i, name in enumerate(padded.row_names)}
    rows = [index[name] for name in fit.row_names]
    np.testing.assert_allclose(padded.U[rows], fit.U, rtol=0, atol=1e-05)


def test_a_condition_level_no_cell_carries_is_refused():
    """It would add an all-zero column and make every coefficient NaN."""
    counts, x, y, cell_type = _raw_inputs()
    with pytest.raises(ValueError, match="no cell carries the condition level"):
        pace_fit(counts, x, y, cell_type, condition=lines("cond2_label.txt"),
                 condition_name="Arm", condition_levels=["PD", "PR", "SD"])


def test_a_numeric_condition_is_refused_rather_than_factored():
    """R would model it continuously; applying contrasts instead is a different fit."""
    counts, x, y, cell_type = _raw_inputs()
    with pytest.raises(ValueError, match="numeric or boolean condition"):
        pace_fit(counts, x, y, cell_type, condition=np.ones(x.size, dtype=int),
                 condition_name="Arm")


def test_inferring_the_reference_arm_warns():
    """The reference arm changes the fit, and byte order is not R's rule."""
    counts, x, y, cell_type = _raw_inputs()
    with pytest.warns(UserWarning, match="reference arm"):
        pace_fit(counts, x, y, cell_type, condition=lines("cond2_label.txt"),
                 condition_name="Arm", n_threads=1)


def test_naming_the_levels_does_not_warn():
    counts, x, y, cell_type = _raw_inputs()
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter("error")
        pace_fit(counts, x, y, cell_type, condition=lines("cond2_label.txt"),
                 condition_name="Arm", condition_levels=["PD", "PR"], n_threads=1)


def test_the_reference_arm_changes_the_fit_not_just_its_labels():
    """Why the warning above exists: this is not a permutation."""
    counts, x, y, cell_type = _raw_inputs()
    labels = ["Responder" if v == "PR" else "non-responder"
              for v in lines("cond2_label.txt")]
    byte_order = pace_fit(counts, x, y, cell_type, condition=labels,
                          condition_name="Resp",
                          condition_levels=["Responder", "non-responder"], n_threads=1)
    collated = pace_fit(counts, x, y, cell_type, condition=labels,
                        condition_name="Resp",
                        condition_levels=["non-responder", "Responder"], n_threads=1)
    assert np.abs(byte_order.U - collated.U).max() > 1.0


def test_the_cores_warnings_reach_the_caller(monkeypatch):
    """R raises each of these; dropping them is what made a NaN fit silent."""
    import pace.fit as module
    real = module._core.fit_irls

    def noisy(**kwargs):
        out = real(**kwargs)
        out["warnings"] = ["non-finite genes at iteration 2"]
        return out

    monkeypatch.setattr(module._core, "fit_irls", noisy)
    counts, x, y, cell_type = _raw_inputs()
    with pytest.warns(UserWarning, match="non-finite genes"):
        pace_fit(counts, x, y, cell_type, n_threads=1)


def test_a_diverged_fit_raises_instead_of_returning_nan(monkeypatch):
    """R's final pass refuses a non-finite eta; Python has no final pass.

    The three inputs known to cause this are refused earlier now, so this
    forces the state directly: the point is that there is a backstop at all,
    for the routes nobody has found yet.
    """
    import pace.fit as module
    real = module._core.fit_irls

    def diverged(**kwargs):
        out = real(**kwargs)
        out["B"] = np.asfortranarray(np.full_like(out["B"], np.nan))
        return out

    monkeypatch.setattr(module._core, "fit_irls", diverged)
    counts, x, y, cell_type = _raw_inputs()
    with pytest.raises(FloatingPointError, match="finite coefficients"):
        pace_fit(counts, x, y, cell_type, n_threads=1)


def test_rows_with_no_data_behind_them_are_counted(fit):
    counts, x, y, cell_type = _raw_inputs()
    assert fit.params["n_structurally_zero"] == 0
    padded = pace_fit(counts, x, y, cell_type, types=sorted(set(cell_type)) + ["Ghost"],
                      n_threads=1)
    assert padded.params["n_structurally_zero"] == 10


# --- streaming the ambient field -----------------------------------------

def test_streaming_the_ambient_field_changes_no_number(fit):
    """`ambient_mode` is a memory knob, so it must not touch the answer.

    Cache mode forms its product with the same routine the streamed path uses
    per chunk, so the two differ only in chunking, which is exact.
    """
    counts, x, y, cell_type = _raw_inputs()
    streamed = pace_fit(counts, x, y, cell_type, ambient_mode="stream", n_threads=1)
    assert np.array_equal(streamed.U, fit.U)
    assert np.array_equal(streamed.alpha, fit.alpha)
    assert np.array_equal(streamed.rho, fit.rho)


def test_a_streamed_condition_fit_also_agrees():
    counts, x, y, cell_type = _raw_inputs()
    shared = dict(condition=lines("cond2_label.txt"), condition_name="Arm",
                  condition_levels=["PD", "PR"], n_threads=1)
    assert np.array_equal(pace_fit(counts, x, y, cell_type, ambient_mode="stream", **shared).U,
                          pace_fit(counts, x, y, cell_type, ambient_mode="cache", **shared).U)


def test_the_readouts_work_off_a_streamed_fit():
    """Streaming must not cost the fit the pieces the readouts rebuild from."""
    counts, x, y, cell_type = _raw_inputs()
    streamed = pace_fit(counts, x, y, cell_type, genes=lines("raw_genes.txt"),
                        ambient_mode="stream", n_threads=1)
    stats = pace.fit_statistics(streamed, counts)
    np.testing.assert_allclose(
        stats.mu_mean,
        vec("stats_mu_mean.txt").reshape(stats.mu_mean.shape, order="F"),
        rtol=1e-06, atol=1e-06)


def test_an_unknown_ambient_mode_is_refused():
    counts, x, y, cell_type = _raw_inputs()
    with pytest.raises(ValueError, match="ambient_mode"):
        pace_fit(counts, x, y, cell_type, ambient_mode="lazy")
