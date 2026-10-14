"""A full PACE run in Python, on the breast cancer Xenium crop the R package ships.

Builds an AnnData from the packaged crop, fits, and draws the figures. Every
number printed here is produced by the Python package alone.

    python examples/demo.py
"""
import gzip
from pathlib import Path

import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

import pace

FIX = Path(__file__).parent.parent / "tests" / "fixtures"
OUT = Path(__file__).parent


def column(name, dtype=float):
    return np.loadtxt(FIX / (name + ".gz"), dtype=dtype, ndmin=1)


def labels(name):
    with gzip.open(FIX / (name + ".gz"), "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


# --- the data, as a user would hold it -------------------------------------
shape = column("raw_y_dim.txt", int)
counts = csc_matrix(
    (column("raw_y_x.txt"), column("raw_y_i.txt", np.int32), column("raw_y_p.txt", np.int32)),
    shape=(int(shape[0]), int(shape[1])))
genes = labels("raw_genes.txt")
cell_type = labels("raw_celltype.txt")

adata = anndata.AnnData(
    X=counts.tocsr(),
    obs=pd.DataFrame({"cell_type": pd.Categorical(cell_type)},
                     index=[f"cell{i}" for i in range(counts.shape[0])]),
    var=pd.DataFrame(index=genes))
adata.obsm["spatial"] = np.column_stack([column("raw_x.txt"), column("raw_y.txt")])
print(adata)
print("cells per type:")
print(adata.obs["cell_type"].value_counts().sort_index().to_string(), "\n")

# --- the pipeline ----------------------------------------------------------
fit = pace.fit_anndata(adata, "cell_type", n_threads=1)
print(fit, "\n")

slopes = pace.pace_shrink(fit, null_correlation=False)
print(slopes)
print("  ", slopes.calls(), "\n")

stats = pace.fit_statistics(fit, counts)
blocks = pace.pace_decompose(fit, counts, stats=stats)
frame = pace.single_frame_decomposition(fit, counts, blocks)
drivers = pace.pace_drivers(fit, slopes, blocks, stats)
pairs = pace.pair_variance(fit, counts, frame)
print(frame, "\n", drivers, "\n")

print("variance by focal type (observed, %):")
print(frame.by_focal().round(2).to_string(index=False), "\n")

strongest = np.argsort(pairs["val"])[::-1][:5]
print("strongest spatial pairs:")
for i in strongest:
    print(f"  {pairs['focal'][i]:<16s} <- {pairs['neighbour'][i]:<16s} {pairs['val'][i]:5.2f}%")

print("\nmacrophage drivers near tumour:")
print(drivers["Macrophage", "Tumour"].to_dataframe()
      .loc[:, ["rank", "gene", "driver_score", "b_clean", "spec", "lfsr"]]
      .round(4).to_string(index=False))

# --- the figures -----------------------------------------------------------
figure = plt.figure(figsize=(16.5, 9))
pace.plot_decomposition(frame, ax=figure.add_subplot(2, 3, 1), zoom=False,
                        title="Variance by focal type")
small = figure.add_subplot(2, 3, 2)
names = [n for n in frame.block_names if n not in ("Cell type %", "Residual %")]
width = 0.8 / len(names)
positions = np.arange(len(frame.focal_levels))
for slot, name in enumerate(names):
    small.bar(positions + slot * width - 0.4 + width / 2,
              frame.focal_table[:, frame.block_names.index(name)], width,
              label=name.replace(" %", ""), color=pace.plots.BLOCK_COLOURS[name])
small.set_xticks(positions, frame.focal_levels, rotation=45, ha="right")
small.set_ylabel("% of observed variance")
small.set_title("the spatial and spillover blocks, rescaled", loc="left")
small.legend(frameon=False, fontsize=8)
small.spines[["top", "right"]].set_visible(False)

pace.plot_pair_heatmap(pairs, cell_types=fit.cell_types, ax=figure.add_subplot(2, 3, 3))
pace.plot_drivers(drivers, "Macrophage", "Tumour", ax=figure.add_subplot(2, 3, 4), n=8)
pace.plot_proximity(fit, counts, ["MRC1", "APOC1"], "Macrophage", "Tumour",
                    ax=figure.add_subplot(2, 3, 5))
pace.plot_shrinkage(slopes, ax=figure.add_subplot(2, 3, 6))
figure.tight_layout()
figure.savefig(OUT / "pace_demo.png", dpi=110)
print("\nwrote", OUT / "pace_demo.png")
