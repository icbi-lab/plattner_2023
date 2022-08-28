# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:.conda-2020-organoids-scanpy]
#     language: python
#     name: conda-env-.conda-2020-organoids-scanpy-py
# ---

# %%
import scvelo as scv
import scanpy as sc
from multiprocessing import Pool
import anndata
import numpy as np
import pandas as pd
from matplotlib import rcParams

rcParams["axes.grid"] = False

# %%
# set this to the data directory where you extracted the data from zenodo
data_dir = "/data/projects/2017/Organoids-ICBI/zenodo/scrnaseq/"

# %%
adata = sc.read_h5ad(f"{data_dir}/04_cell_types/adata_cell_type.h5ad")


# %%
def _read_scvelo(organoid):
    adata = scv.read_loom(
        f"{data_dir}/05_scvelo/organoids-batch1-{organoid}/velocyto/organoids-batch1-{organoid}.loom"
    )
    adata.obs["organoid"] = organoid
    adata.var_names_make_unique()
    return adata


with Pool(16) as p:
    adatas = p.map(_read_scvelo, adata.obs["organoid"].unique())

for tmp_adata in adatas:
    tmp_adata.obs_names = [
        x.split(":")[1].replace("x", "") + "-" + organoid
        for x, organoid in zip(tmp_adata.obs_names, tmp_adata.obs["organoid"])
    ]

# %%
adata_scvelo = anndata.concat(adatas)

# %%
# relabel cells by organoid as batch id in `obs_names`.
adata.obs_names = [
    x.split("-")[0] + "-" + organoid
    for x, organoid in zip(adata.obs_names, adata.obs["organoid"])
]

# %%
# This only affects ~10 cells in CRC03 and this only arises by matching the barcode names
# between anndata_scvelo and anndata.
adata = adata[~adata.obs_names.duplicated(), :]

# %%
# Number of common cells in transcriptomics and scvelo anndata objects.
len(set(adata.obs_names)), len(set(adata_scvelo.obs_names)), len(
    set(adata.obs_names) & set(adata_scvelo.obs_names)
)

# %% [markdown]
# ### Preprocess anndata scvelo object

# %%
adata_scvelo = scv.utils.merge(adata_scvelo, adata)

# %%
scv.pp.filter_and_normalize(adata_scvelo)

# %%
scv.pp.moments(adata_scvelo, n_pcs=30, n_neighbors=30)

# %%
scv.tl.velocity(adata_scvelo)

# %%
scv.tl.velocity_graph(adata_scvelo)

# %% [markdown]
# ## Plots

# %%
scv.pl.proportions(adata_scvelo, groupby="organoid")

# %%
sc.set_figure_params(figsize=(10, 10))
rcParams["axes.grid"] = False

# %%
ax = sc.pl.embedding(
    adata,
    basis="umap_uncorrected",
    color="cell_type",
    show=False,
    legend_loc="None",
    size=45,
    alpha=0.3,
)
scv.pl.velocity_embedding_stream(
    adata_scvelo,
    basis="umap_uncorrected",
    color="cell_type",
    #     arrow_color="white",
    legend_loc="right margin",
    ax=ax,
    alpha=0,
)

# %%
ax = sc.pl.embedding(
    adata,
    basis="umap",
    color="cell_type",
    show=False,
    legend_loc="None",
    size=45,
    alpha=0.3,
)
scv.pl.velocity_embedding_stream(
    adata_scvelo,
    basis="umap",
    color="cell_type",
    #     arrow_color="white",
    legend_loc="right margin",
    ax=ax,
    alpha=0,
)
ax.get_figure().savefig("figures/umap_scvelo_cell_types.svg", dpi=600, bbox_inches="tight")
