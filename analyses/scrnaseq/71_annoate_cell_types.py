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
import scanpy as sc
import scvi
import progeny
import dorothea
from scanpy_helpers.annotation import AnnotationHelper
from scanpy_helpers import de
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
import altair as alt

sc.set_figure_params(figsize=(5, 5))

# %%
# set this to the data directory where you extracted the data from zenodo
data_dir = "/data/projects/2017/Organoids-ICBI/zenodo/scrnaseq/"
res_dir = "../../results"

# %%
# !mkdir -p {res_dir}

# %%
adata = sc.read_h5ad(f"{data_dir}/03_scvi/adata_integrated.h5ad")
model_all_genes = scvi.model.SCVI.load(f"{data_dir}/03_scvi/scvi_model/", adata)

# %%
adata.shape

# %%
adata.obs.groupby("organoid").size()

# %%
sc.pl.dotplot(
    adata,
    var_names={
        "Enterocyte": ["NDRG1", "EMP1", "TFF2", "SDCBP2", "CEACAM7"],
        "Enteroendocrine": ["MYH7B", "NOTUM", "DEFA6", "APCDD1", "NKD1"],
        "Goblet-like": ["SLC26A3", "CEACAM7", "CEACAM6", "TFF1", "MUC17"],
        "M-cell": ["MT2A", "CAV1", "HS3ST1", "KITLG", "CALD1"],
        "Stem-like": ["PRSS2", "SMOC2", "PLA2G2A", "MYC", "MXD4"],
        "Stem/TA-like": ["CDC20", "CCNB1", "PTTG1", "BIRC5", "CDKN3"],
        "TA-like": ["HIST1H1B", "CDK1", "UBE2C", "SPC25", "HJURP"],
        "dividing": ["UBE2C", "CDK1", "PLK1", "CENPA", "CCNB1"],
    },
    groupby="organoid",
    swap_axes=True,
)

# %%
sc.pl.dotplot(
    adata,
    var_names=[
        x.strip()
        for x in "LGR5, ASCL2, AXIN1, AXIN2, CD44, EPHB2, SMAD2, SMAD3, TFF1, TFF3, CD274".split(
            ","
        )
    ],
    groupby="organoid",
    swap_axes=True,
)

# %%
sc.pl.umap(adata, color=["organoid", "leiden"], wspace=0.5)

# %%
sc.pl.embedding(adata, "umap_uncorrected", color=["organoid", "leiden"], wspace=0.5)

# %% [markdown]
# ## scVI DE

# %%
try:
    de_res = pd.read_csv("tmp/de_res.csv", index_col=0)
except IOError:
    de_res = de.scvi(adata, model_all_genes, groupby="leiden")
    de_res.to_csv("tmp/de_res.csv")

# %%
de_res["score"] = de_res["lfc_mean"] * (de_res["raw_normalized_mean1"] > 0.5)
de_res["leiden"] = [x.split()[0] for x in de_res["comparison"]]

# %%
de_res["gene_symbol"] = de_res.index

# %%
de.de_res_to_anndata(
    adata,
    de_res,
    groupby="leiden",
    pval_col="proba_not_de",
    pval_adj_col="proba_not_de",
    lfc_col="lfc_mean",
)

# %%
sc.pl.rank_genes_groups_dotplot(adata, swap_axes=True, n_genes=5)

# %% [markdown]
# ## Annotate cell-types

# %%
ah = AnnotationHelper(sheet="organoids")

# %%
ah.plot_umap(adata)

# %%
adata.shape

# %%
sc.pl.umap(adata, color="BIRC5", cmap="inferno")

# %%
ah.annotate_cell_types(
    adata,
    {
        "dividing": [10],
        "M-cell like": [12],
        "Enteroendocrine-like (immature)": [3],
        "Stem/TA-like": [2, 8],
        "Stem-like": [4, 7, 9],
        "TA-like": [0, 5],
        "Enterocyte/Goblet-like (immature)": [11, 1, 6],
    },
)

# %%
adata_sec_prec = adata[adata.obs["cell_type"] == "Enterocyte/Goblet-like (immature)", :]

# %%
ah.reprocess_adata_subset_scvi(adata_sec_prec, leiden_res=0.5)

# %%
sc.pl.umap(adata_sec_prec, color="leiden")

# %%
ah.plot_umap(adata_sec_prec, filter_cell_type=["Goblet", "Enterocyte"])

# %%
ah.annotate_cell_types(
    adata_sec_prec,
    {
        "Enterocyte/Goblet-like (immature)": list(range(12)) + [13],
        "Goblet-like": [12],
    },
)

# %%
ah.integrate_back(adata, adata_sec_prec)

# %% [markdown]
# ## Cell-type overview plots

# %%
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": (300)}):
    sc.pl.umap(
        adata,
        color="cell_type",
        show=True,
        frameon=False,
        add_outline=True,
        title="",
        size=30,
        save="_integrated_cell_type.pdf",
    )

# %%
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": (300)}):
    sc.pl.umap(
        adata,
        color="organoid",
        show=True,
        frameon=False,
        add_outline=True,
        title="",
        size=30,
        save="_integrated_organoid.pdf",
    )

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    fig = sc.pl.embedding(
        adata,
        "umap_uncorrected",
        color="organoid",
        show=True,
        frameon=False,
        add_outline=True,
        title="",
        size=30,
        return_fig=True,
        legend_loc="on data",
        legend_fontoutline=3,
    )
    fig.savefig(f"figures/umap_uncorrected_cell_type.svg", dpi=600, bbox_inches="tight")

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    fig = sc.pl.embedding(
        adata,
        "umap_uncorrected",
        color="cell_type",
        show=True,
        frameon=False,
        add_outline=True,
        title="",
        size=30,
        return_fig=True,
    )
    fig.savefig(f"figures/umap_uncorrected_cell_type.svg", dpi=600, bbox_inches="tight")


# %%
with plt.rc_context({"figure.figsize": (3, 3)}):
    fig = sc.pl.embedding(
        adata,
        "umap",
        color=["LGR5", "AXIN2", "TFF3", "FABP1"],
        show=True,
        frameon=False,
        add_outline=True,
        size=5,
        cmap="coolwarm",
        ncols=2,
        return_fig=True,
    )
    fig.savefig(f"figures/umap_cell_type_markers.svg", dpi=600, bbox_inches="tight")

# %%
cell_counts = (
    adata.obs.groupby(["organoid", "cell_type"]).size().reset_index(name="count")
)

# %%
chart = (
    alt.Chart(cell_counts)
    .mark_bar()
    .encode(
        y=alt.Y("sum(count)", stack="normalize", title="cell-type fraction"),
        x="organoid",
        color=alt.Color(
            "cell_type",
            scale=alt.Scale(
                domain=list(adata.obs["cell_type"].cat.categories),
                range=adata.uns["cell_type_colors"],
            ),
        ),
    )
)
chart

# %%
chart.save("figures/cell_type_proportions.svg", dpi=300)

# %% [markdown]
# ## Plot genes of interest

# %%
genes_of_interest = [
    "LGR5",
    "MYC",
    "KRT20",
    "FABP1",
    "SPINK4",
    "AURKA",
    "PRKCA",
    "STAT3",
]

# %%
sc.pl.stacked_violin(
    adata,
    var_names=genes_of_interest,
    groupby="organoid",
    cmap="bwr",
    swap_axes=True,
    normalize="var",
    save="organoid.svg",
)

# %% [markdown]
# ## scVI DE of cell-types

# %%
try:
    de_res_cell_type = pd.read_csv("tmp/de_res_cell_type.csv", index_col=0)
except IOError:
    de_res_cell_type = de.scvi(adata, model_all_genes, groupby="cell_type")
    de_res_cell_type.to_csv("tmp/de_res_cell_type.csv")

# %%
de_res_cell_type["score"] = de_res_cell_type["lfc_mean"] * (
    de_res_cell_type["raw_normalized_mean1"] > 0.5
)
de_res_cell_type["leiden"] = [x.split()[0] for x in de_res_cell_type["comparison"]]

# %%
de_res_cell_type["gene_symbol"] = de_res_cell_type.index

# %%
de.de_res_to_anndata(
    adata,
    de_res_cell_type.rename(columns={"leiden": "cell_type"}),
    groupby="cell_type",
    pval_col="proba_not_de",
    pval_adj_col="proba_not_de",
    lfc_col="lfc_mean",
)

# %%
fig = sc.pl.rank_genes_groups_dotplot(adata, swap_axes=True, n_genes=5, return_fig=True)
fig.savefig("figures/dotplot_cell_type_markers.svg", bbox_inches="tight")

# %%
sc.pl.dotplot(
    adata,
    swap_axes=True,
    groupby="cell_type",
    var_names=["NDRG1", "NKD1", "TFF1", "MT2A", "MXD4", "MYC", "CDK1"],
)

# %% [markdown]
# ## Write results

# %%
del adata.uns["rank_genes_groups"]

# %%
adata.write_h5ad(f"{res_dir}/adata_cell_type.h5ad")

# %%
de_res.to_csv(f"{res_dir}/de_result.csv")
de_res_cell_type.to_csv(f"{res_dir}/de_result_cell_type.csv")
