# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: 'SSH apollo-15 apollo-15: 2020-organoids-scanpy'
#     language: ''
#     name: rik_ssh_apollo_15_apollo152020organoidsscanpy
# ---

# %%
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import anndata
import infercnvpy as cnv

sc.set_figure_params(figsize=(5, 5))
from scanpy_helpers.annotation import AnnotationHelper

# %%
scvi.settings.seed = 0

# %%
organoids = ["CRC02", "CRC03", "CRC03_2", "CRC04", "CRC13", "CRC26", "CRC26LM"]
adatas = {
    organoid: sc.read_h5ad(
        f"../data/scrnaseq/01_qc_and_filtering/{organoid}/{organoid}.qc.h5ad"
    )
    for organoid in organoids
}
doublets = {
    organoid: pd.read_csv(
        f"../data/scrnaseq/02_solo/{organoid}/{organoid}.is_doublet.csv",
        header=None,
        names=["barcode", "is_doublet"],
    )
    for organoid in organoids
}

# %%
for org, adata in adatas.items(): 
    print(f"{org} {adata.shape[0]}")

# %%
adatas["CRC02"].obs_names

# %%
for adata, doublet in zip(adatas.values(), doublets.values()):
    adata.obs.loc[doublet["barcode"], "is_doublet"] = doublet[
        "is_doublet"
    ].values.astype(str)

# %%
for organoid, adata in adatas.items():
    adata.obs["organoid"] = organoid.split("_")[0]

# %%
adata = anndata.concat(adatas.values(), index_unique="-")

# %%
adata

# %%
sc.pp.highly_variable_genes(
    adata, n_top_genes=5000, flavor="seurat_v3", batch_key="organoid"
)

# %%
adata_scvi = adata.copy()

# %%
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

# %%
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# %%
adata.obs

# %%
sc.pl.umap(
    adata, color=["n_genes", "pct_counts_mito", "organoid", "leiden", "is_doublet"]
)

# %%
sc.pl.umap(adata, color=["organoid"])

# %%
ah = AnnotationHelper(sheet="organoids")

# %% [markdown]
# ## scVI

# %%
adata_scvi = adata_scvi[adata_scvi.obs["is_doublet"] == "False", :].copy()

# %%
adata_scvi_hvg = adata_scvi[:, adata_scvi.var["highly_variable"]].copy()

# %%
scvi.data.setup_anndata(adata_scvi, batch_key="organoid")
scvi.data.view_anndata_setup(adata_scvi)
model = scvi.model.SCVI(adata_scvi)
model.train(max_epochs=700, early_stopping=True)
model.to_device("cpu")

# %%
scvi.data.setup_anndata(adata_scvi_hvg, batch_key="organoid")
scvi.data.view_anndata_setup(adata_scvi_hvg)
model_hvg = scvi.model.SCVI(adata_scvi_hvg)
model_hvg.train(max_epochs=700, early_stopping=True)
model_hvg.to_device("cpu")

# %% [markdown]
# ### Prepare to save models and latent representations

# %%
adata_scvi.raw = adata[adata_scvi.obs_names, :]

# %%
adata_scvi.obsm["X_scVI"] = model_hvg.get_latent_representation()
adata_scvi.obsm["X_scVI_all_genes"] = model.get_latent_representation()

# %%
sc.pp.neighbors(adata_scvi, use_rep="X_scVI")

# %%
sc.tl.umap(adata_scvi)

# %%
sc.tl.leiden(adata_scvi, resolution=0.5)

# %%
adata.uns["neighbors"]

# %%
adata_scvi.uns["neighbors_uncorrected"] = adata.uns["neighbors"]
adata_scvi.uns["neighbors_uncorrected"]["connectivityes_key"] = "connectivities_uncorrected"
adata_scvi.uns["neighbors_uncorrected"]["distances_key"] = "distances_uncorrected"
adata_scvi.obsp["connectivities_uncorrected"] = adata[adata_scvi.obs_names, :].obsp["connectivities"]
adata_scvi.obsp["distances_uncorrected"] = adata[adata_scvi.obs_names, :].obsp["distances"]
adata_scvi.obsm["X_pca"] = adata[adata_scvi.obs_names, :].obsm["X_pca"]
adata_scvi.obsm["X_umap_uncorrected"] = adata[adata_scvi.obs_names, :].obsm["X_umap"]

# %%
# !rm -rvf ../data/scrnaseq/03_scvi

# %%
# !mkdir -p ../data/scrnaseq/03_scvi
model.save("../data/scrnaseq/03_scvi/scvi_model")
model_hvg.save("../data/scrnaseq/03_scvi/scvi_model_hvg")
adata_scvi.write_h5ad("../data/scrnaseq/03_scvi/adata_integrated.h5ad")

# %%
sc.pl.umap(adata_scvi, color=["n_genes", "pct_counts_mito", "organoid", "leiden"])

# %%
sc.pl.embedding(adata_scvi, basis="scVI", color="leiden")

# %%
sc.pl.umap(
    adata_scvi, color=["EPCAM", "SLC26A3", "CPE", "TOP2A"], ncols=2, cmap="inferno"
)
