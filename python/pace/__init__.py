"""PACE: Proximity-Associated Changes in Expression.

This package binds the same compiled core the R package uses, so results do
not depend on which language you call it from.
"""

from .fit import PaceFit, fit_anndata, pace_fit
from .shrink import ShrunkenSlopes, pace_shrink
from .stats import FitStatistics, fit_statistics, fitted_means
from .decompose import (Decomposition, SingleFrameDecomposition, pace_decompose,
                        pair_variance, single_frame_decomposition)
from .drivers import DriverTables, PairDrivers, pace_drivers
from .plots import (plot_decomposition, plot_drivers, plot_pair_heatmap,
                    plot_proximity, plot_shrinkage)
from ._core import (
    ambient_field,
    anchor_mask,
    centre_within_groups,
    data_informed_weights,
    drop_sparse_kernel,
    group_column_means,
    same_type_fraction,
    driver_scores,
    fit_irls,
    mash,
    neighbour_kernels,
    pair_variance_pratt,
    sparse_product,
    standardise,
    random_design_block,
)

__all__ = [
    "mash",
    "driver_scores",
    "pair_variance_pratt",
    "neighbour_kernels",
    "ambient_field",
    "random_design_block",
    "fit_irls",
    "sparse_product",
    "drop_sparse_kernel",
    "centre_within_groups",
    "anchor_mask",
    "data_informed_weights",
    "standardise",
    "group_column_means",
    "same_type_fraction",
    "pace_fit",
    "fit_anndata",
    "PaceFit",
    "pace_shrink",
    "ShrunkenSlopes",
    "fit_statistics",
    "fitted_means",
    "FitStatistics",
    "pace_decompose",
    "Decomposition",
    "single_frame_decomposition",
    "SingleFrameDecomposition",
    "pair_variance",
    "pace_drivers",
    "DriverTables",
    "PairDrivers",
    "plot_decomposition",
    "plot_pair_heatmap",
    "plot_drivers",
    "plot_proximity",
    "plot_shrinkage",
]
__version__ = "0.0.1"
