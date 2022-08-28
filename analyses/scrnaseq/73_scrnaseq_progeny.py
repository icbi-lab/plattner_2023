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
adata.obsm["mlm_estimate"]

# %%
# extract activities:
acts = dc.get_acts(adata, obsm_key="mlm_estimate")
sc.pp.scale(acts)

# %%
# Visualize WNT:
sc.pl.umap(acts, color="WNT", vcenter=0, cmap="coolwarm")

# %%
# Split into organoids:
for org in ["CRC02", "CRC03", "CRC04", "CRC13", "CRC26", "CRC26LM"]:
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
        1, 4, figsize=(24, 6)
    )  # , gridspec_kw={'wspace':0.9})
    ax1_dict = sc.pl.umap(
        adata[adata.obs["organoid"] == org],
        color="cell_type",
        title=org,
        ax=ax1,
        show=False,
        legend_loc=False,
    )
    ax2_dict = sc.pl.umap(
        acts[acts.obs["organoid"] == org],
        color="WNT",
        vcenter=0,
        cmap="coolwarm",
        ax=ax2,
        show=False,
    )
    ax3_dict = sc.pl.violin(
        acts[acts.obs["organoid"] == org],
        keys="WNT",
        groupby="cell_type",
        rotation=90,
        ax=ax3,
        show=False,
    )
    ax4_dict = sc.pl.matrixplot(
        acts[acts.obs["organoid"] == org],
        var_names=acts.var.index,
        groupby="cell_type",
        vmin=-2,
        vmax=2,
        cmap="bwr",
        swap_axes=True,
        ax=ax4,
        show=False,
    )

# %%
mean_acts = dc.summarize_acts(acts, groupby="cell_type", min_std=0)
mean_acts

# %%
sns.clustermap(
    mean_acts, xticklabels=mean_acts.columns, vmin=-2, vmax=2, cmap="coolwarm"
)
plt.show()

# %%
sc.pl.matrixplot(
    acts,
    var_names=acts.var.index,
    groupby="organoid",
    vmin=-2,
    vmax=2,
    cmap="bwr",
    swap_axes=True,
    dendrogram=True
)

# %%
sc.pl.matrixplot(
    acts,
    var_names=acts.var.index,
    groupby="cell_type",
    vmin=-2,
    vmax=2,
    cmap="bwr",
    swap_axes=True,
    dendrogram=True,
)

# %%
acts_old = dc.get_acts(adata, obsm_key="progeny")

# %%
sc.pl.umap(acts_old, color="WNT", vcenter=0, cmap="coolwarm")

# %%
mean_acts_old = dc.summarize_acts(acts_old, groupby="cell_type", min_std=0)
mean_acts_old

# %%
sns.clustermap(
    mean_acts_old, xticklabels=mean_acts_old.columns, vmin=-2, vmax=2, cmap="coolwarm"
)
plt.show()

# %%
# Split into organoids:
for org in ["CRC02", "CRC03", "CRC04", "CRC13", "CRC26", "CRC26LM"]:
    # fig, (ax1, ax2, ax3, ax4) = plt.subplots(1,4, figsize = (24,6))
    fig, (ax1, ax2, ax4) = plt.subplots(1, 3, figsize=(18, 6))
    plt.gcf().subplots_adjust(bottom=0.42)
    ax1_dict = sc.pl.umap(
        adata[adata.obs["organoid"] == org],
        color="cell_type",
        title=org,
        ax=ax1,
        show=False,
        legend_loc=False,
    )
    ax2_dict = sc.pl.umap(
        acts_old[acts_old.obs["organoid"] == org],
        color="WNT",
        vcenter=0,
        cmap="coolwarm",
        ax=ax2,
        show=False,
    )
    # ax3_dict = sc.pl.violin(acts_old[acts_old.obs['organoid']==org], keys='WNT', groupby='cell_type', rotation=90, ax=ax3, show=False)
    ax4_dict = sc.pl.matrixplot(
        acts_old[acts_old.obs["organoid"] == org],
        var_names=acts.var.index,
        groupby="cell_type",
        vmin=-2,
        vmax=2,
        cmap="bwr",
        swap_axes=True,
        ax=ax4,
        show=False,
    )
    plt.savefig(f"{org}_progeny_WNT.pdf")

# %%
adata.obs.cell_type

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


# %%
