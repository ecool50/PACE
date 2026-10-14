"""`fit_anndata()`, the reader for callers holding an AnnData.

It must agree with the array path exactly, since it only rearranges inputs,
and it must refuse the things R never has to think about: `counts(spe)` names
the assay it wants, where `adata.X` is whatever preprocessing left behind.
"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

import pace

anndata = pytest.importorskip("anndata")
pd = pytest.importorskip("pandas")

FIX = Path(__file__).parent / "fixtures"


def vec(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def lines(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


@pytest.fixture(scope="module")
def arrays():
    dim = vec("raw_y_dim.txt", int)
    counts = csc_matrix(
        (vec("raw_y_x.txt"), vec("raw_y_i.txt", np.int32), vec("raw_y_p.txt", np.int32)),
        shape=(int(dim[0]), int(dim[1])))
    genes = [f"g{i}" for i in range(counts.shape[1])]
    return (counts, vec("raw_x.txt"), vec("raw_y.txt"), lines("raw_celltype.txt"),
            lines("two_image_label.txt"), genes)


@pytest.fixture(scope="module")
def adata(arrays):
    counts, x, y, cell_type, panel, genes = arrays
    out = anndata.AnnData(
        X=counts.tocsr().astype(np.float32),
        obs=pd.DataFrame({"cellType": pd.Categorical(cell_type), "panel": panel},
                         index=[f"c{i}" for i in range(counts.shape[0])]),
        var=pd.DataFrame(index=genes))
    out.obsm["spatial"] = np.column_stack([x, y])
    return out


def test_it_agrees_with_the_array_path_exactly(adata, arrays):
    counts, x, y, cell_type, _, genes = arrays
    from_arrays = pace.pace_fit(counts, x, y, cell_type, genes=genes, n_threads=1)
    from_adata = pace.fit_anndata(adata, "cellType", n_threads=1)
    assert np.array_equal(from_adata.U, from_arrays.U)
    assert np.array_equal(from_adata.alpha, from_arrays.alpha)
    assert from_adata.genes == from_arrays.genes
    assert from_adata.row_names == from_arrays.row_names


def test_it_is_still_bit_identical_to_r(adata):
    fit = pace.fit_anndata(adata, "cellType", n_threads=1)
    assert np.array_equal(fit.U, vec("irls_ab_U.txt").reshape(fit.U.shape, order="F"))


def test_the_gene_names_come_from_var_names(adata):
    fit = pace.fit_anndata(adata, "cellType", n_threads=1)
    assert fit.genes == list(adata.var_names)


def test_an_image_key_is_honoured(adata):
    fit = pace.fit_anndata(adata, "cellType", image_key="panel", n_threads=1)
    reference = vec("two_image_global_U.txt").reshape(fit.U.shape, order="F")
    np.testing.assert_allclose(fit.U, reference, rtol=0, atol=1e-06)


def test_a_dense_x_works_too(adata, arrays):
    counts = arrays[0]
    dense = adata.copy()
    dense.X = np.asarray(counts.todense(), dtype=np.float32)
    assert np.array_equal(pace.fit_anndata(dense, "cellType", n_threads=1).U,
                          pace.fit_anndata(adata, "cellType", n_threads=1).U)


def test_a_layer_can_hold_the_counts(adata):
    moved = adata.copy()
    moved.layers["counts"] = moved.X.copy()
    moved.X = moved.X.copy()
    moved.X.data = np.log1p(moved.X.data)       # the usual state of X
    fit = pace.fit_anndata(moved, "cellType", layer="counts", n_threads=1)
    assert np.array_equal(fit.U, pace.fit_anndata(adata, "cellType", n_threads=1).U)


def test_the_adata_is_not_modified(adata):
    before_obs = list(adata.obs.columns)
    before_obsm = list(adata.obsm)
    pace.fit_anndata(adata, "cellType", n_threads=1)
    assert list(adata.obs.columns) == before_obs
    assert list(adata.obsm) == before_obsm


# --- what it refuses -----------------------------------------------------

def test_normalised_counts_are_refused_rather_than_fitted(adata):
    """The one check R never needs, and the costliest to get wrong."""
    logged = adata.copy()
    logged.X = logged.X.copy()
    logged.X.data = np.log1p(logged.X.data)
    with pytest.raises(ValueError, match="not integer-valued"):
        pace.fit_anndata(logged, "cellType", n_threads=1)


def test_negative_values_are_refused(adata):
    scaled = adata.copy()
    scaled.X = scaled.X.copy()
    scaled.X.data = scaled.X.data - 5.0
    with pytest.raises(ValueError, match="negative values"):
        pace.fit_anndata(scaled, "cellType", n_threads=1)


def test_missing_coordinates_are_named_not_guessed(adata):
    bare = adata.copy()
    del bare.obsm["spatial"]
    with pytest.raises(KeyError, match="spatial"):
        pace.fit_anndata(bare, "cellType", n_threads=1)


def test_three_dimensional_coordinates_are_refused_not_projected(adata):
    stacked = adata.copy()
    stacked.obsm["spatial"] = np.column_stack(
        [stacked.obsm["spatial"], np.zeros(stacked.n_obs)])
    with pytest.raises(ValueError, match="two-dimensional"):
        pace.fit_anndata(stacked, "cellType", n_threads=1)


def test_a_missing_label_stays_missing_rather_than_becoming_a_type(adata):
    """`.astype(str)` would make this a cell type called "nan"."""
    holed = adata.copy()
    labels = holed.obs["cellType"].astype(object).to_numpy()
    labels[11] = np.nan
    holed.obs["cellType"] = labels
    with pytest.raises(ValueError, match="every cell needs a cell type"):
        pace.fit_anndata(holed, "cellType", n_threads=1)


def test_an_unknown_obs_key_is_named(adata):
    with pytest.raises(KeyError, match="celltype"):
        pace.fit_anndata(adata, "celltype", n_threads=1)


def test_an_unknown_layer_is_named(adata):
    with pytest.raises(KeyError, match="raw"):
        pace.fit_anndata(adata, "cellType", layer="raw", n_threads=1)


def test_pace_fit_points_an_anndata_at_the_right_function(adata):
    with pytest.raises(TypeError, match="fit_anndata"):
        pace.pace_fit(adata, [0.0], [0.0], ["Tumour"])


def test_duplicated_gene_names_are_refused(adata):
    """R refuses these: two genes would share one row of the output."""
    repeated = adata.copy()
    names = list(repeated.var_names)
    names[4] = names[0]
    repeated.var.index = pd.Index(names)
    with pytest.raises(ValueError, match="var_names must be unique"):
        pace.fit_anndata(repeated, "cellType", n_threads=1)


def test_non_finite_coordinates_name_the_obsm_key(adata):
    holed = adata.copy()
    coordinates = holed.obsm["spatial"].copy()
    coordinates[13, 1] = np.nan
    holed.obsm["spatial"] = coordinates
    with pytest.raises(ValueError, match="non-finite coordinates, first at cell 13"):
        pace.fit_anndata(holed, "cellType", n_threads=1)


def test_a_backed_object_is_named_rather_than_mis_diagnosed(adata, tmp_path):
    """Backed sparse X is not a scipy matrix, and would look like a shape error."""
    path = tmp_path / "backed.h5ad"
    adata.write_h5ad(path)
    backed = anndata.read_h5ad(path, backed="r")
    with pytest.raises(ValueError, match="backed"):
        pace.fit_anndata(backed, "cellType", n_threads=1)


def test_gene_names_are_strings_even_when_var_names_are_not(adata):
    numbered = adata.copy()
    numbered.var.index = pd.Index(range(numbered.n_vars))
    fit = pace.fit_anndata(numbered, "cellType", n_threads=1)
    assert all(isinstance(name, str) for name in fit.genes)


def test_a_missing_obsm_key_is_reported_before_the_counts_are_converted(adata):
    """The cheap lookup should not wait on the expensive conversion."""
    bare = adata.copy()
    del bare.obsm["spatial"]
    bare.X = bare.X.copy()
    bare.X.data = np.log1p(bare.X.data)       # would also fail the counts check
    with pytest.raises(KeyError, match="spatial"):
        pace.fit_anndata(bare, "cellType", n_threads=1)


def test_the_counts_arrays_are_not_written_through(adata):
    """The core is handed adata.X's own index arrays when X is already CSC."""
    held = adata.copy()
    held.X = held.X.tocsc()
    before = (held.X.data.copy(), held.X.indices.copy(), held.X.indptr.copy())
    pace.fit_anndata(held, "cellType", n_threads=1)
    assert np.array_equal(held.X.data, before[0])
    assert np.array_equal(held.X.indices, before[1])
    assert np.array_equal(held.X.indptr, before[2])


# --- a condition cohort through AnnData ----------------------------------

def test_a_condition_key_matches_r(adata):
    conditioned = adata.copy()
    conditioned.obs["Arm"] = lines("cond2_label.txt")
    fit = pace.fit_anndata(conditioned, "cellType", condition_key="Arm", n_threads=1)
    np.testing.assert_allclose(
        fit.U, vec("cond2_U.txt").reshape(fit.U.shape, order="F"), rtol=0, atol=1e-06)
    assert fit.row_names == lines("cond2_rownames.txt")


def test_the_obs_key_becomes_the_term_prefix(adata):
    """As the colData column name does in R."""
    conditioned = adata.copy()
    conditioned.obs["Responder"] = lines("cond2_label.txt")
    fit = pace.fit_anndata(conditioned, "cellType", condition_key="Responder", n_threads=1)
    assert fit.params["fixed_names"] == ["(Intercept)", "ResponderPR"]


def test_a_categorical_condition_keeps_its_category_order(adata):
    """R honours a factor's level order; pandas categories are the same idea."""
    conditioned = adata.copy()
    conditioned.obs["Arm"] = pd.Categorical(
        lines("cond2_label.txt"), categories=["PR", "PD"])
    fit = pace.fit_anndata(conditioned, "cellType", condition_key="Arm", n_threads=1)
    assert fit.params["condition_levels"] == ["PR", "PD"]
    assert fit.params["fixed_names"] == ["(Intercept)", "ArmPD"]


def test_an_unused_category_is_dropped_as_r_drops_it(adata):
    conditioned = adata.copy()
    conditioned.obs["Arm"] = pd.Categorical(
        lines("cond2_label.txt"), categories=["PD", "PR", "CR"])
    fit = pace.fit_anndata(conditioned, "cellType", condition_key="Arm", n_threads=1)
    assert fit.params["condition_levels"] == ["PD", "PR"]
