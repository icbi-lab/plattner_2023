# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:conda-test-decoupler]
#     language: python
#     name: conda-env-conda-test-decoupler-py
# ---

# %%
import decoupler as dc

# visualization:
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

# %%
# set this to the data directory where you extracted the data from zenodo
data_dir = "/data/projects/2017/Organoids-ICBI/zenodo/scrnaseq/"

# %%
adata = sc.read_h5ad(
    f"{data_dir}/04_cell_types/adata_cell_type.h5ad"
)
adata

# %%
sc.pl.umap(adata, color="cell_type", title="UMAP", frameon=False, save="_celltypes.pdf")

# %%
# PROGENy model:
model = dc.get_progeny(organism="human", top=1000)
model

# %%
# Activity interference with Multivariate Linear Model:
dc.run_mlm(
    mat=adata,
    net=model,
    source="source",
    target="target",
    weight="weight",
    verbose=True,
    use_raw=True,
)

# %%
# Rename cell types (delete "immature"):
new_cell_types = [
    "Enterocyte/Goblet-like",
    "Enteroendocrine-like",
    "Goblet-like",
    "M-cell like",
    "Stem-like",
    "Stem/TA-like",
    "TA-like",
    "dividing",
]
adata.rename_categories("cell_type", new_cell_types)

# %%
# Remove M- and Goblet cells for comparison (are not present in all organoids):
adata_subset_plot = adata[(adata.obs.cell_type != "M-cell like") & (adata.obs.cell_type != "Goblet-like")]

# %%
sc.pl.umap(
    adata_subset_plot, color="cell_type", title="UMAP", frameon=False, save="_celltypes_new.pdf"
)

# %%
# Split into organoids:
acts_old = dc.get_acts(adata_subset_plot, obsm_key="progeny")

for org in ["CRC02", "CRC03", "CRC04", "CRC13", "CRC26", "CRC26LM"]:
    ax4_dict = sc.pl.matrixplot(
        acts_old[acts_old.obs["organoid"] == org],
        var_names=acts_old.var.index,
        groupby="cell_type",
        vmin=-2,
        vmax=2,
        cmap="bwr",
        swap_axes=True,
        title = org,
        save=(f"{org}_progeny_heatmap.pdf")
    )

