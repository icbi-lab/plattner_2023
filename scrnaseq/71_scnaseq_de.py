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
adata = sc.read_h5ad("../data/scrnaseq/03_scvi/adata_integrated.h5ad")
model_all_genes = scvi.model.SCVI.load("../data/scrnaseq/03_scvi/scvi_model/", adata)

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
        return_fig = True,
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
        return_fig = True
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
        return_fig=True
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

# %%
sc.pl.dotplot(
    adata,
    var_names=genes_of_interest,
    groupby="organoid",
    swap_axes=True,
    save="organoid.pdf",
)

# %%
sc.pl.stacked_violin(
    adata,
    var_names=genes_of_interest,
    groupby="cell_type",
    cmap="bwr",
    swap_axes=True,
    normalize="var",
    save="cell_type.pdf",
)

# %%
sc.pl.dotplot(
    adata,
    var_names=genes_of_interest,
    groupby="cell_type",
    swap_axes=True,
    save="cell_type.pdf",
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
fig =sc.pl.rank_genes_groups_dotplot(adata, swap_axes=True, n_genes=5, return_fig=True)
fig.savefig("figures/dotplot_cell_type_markers.svg", bbox_inches="tight")

# %%
sc.pl.dotplot(
    adata,
    swap_axes=True,
    groupby="cell_type",
    var_names=["NDRG1", "NKD1", "TFF1", "MT2A", "MXD4", "MYC", "CDK1"],
)

# %% [markdown]
# ## Progeny

# %%
progeny_model = progeny.load_model(organism="Human", top=1000)

# %%
progeny.run(
    adata,  # Data to use
    progeny_model,  # PROGENy network
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # Pathways with less than 5 targets will be ignored
)

# %%
ad_progeny = progeny.extract(adata)

# %%
sc.pl.matrixplot(
    ad_progeny,
    var_names=ad_progeny.var.index,
    groupby="organoid",
    vmin=-2,
    vmax=2,
    cmap="bwr",
    swap_axes=True,
    dendrogram=True,
    save="_progeny_organoids.pdf",
)

# %%
sc.pl.matrixplot(
    ad_progeny,
    var_names=ad_progeny.var.index,
    groupby="cell_type",
    vmin=-2,
    vmax=2,
    cmap="bwr",
    swap_axes=True,
    dendrogram=True,
    save="_progeny_cell_types.pdf",
)

# %%
for pw in ["MAPK", "PI3K", "WNT"]:
    tmp_ad = sc.AnnData(
        pd.DataFrame(ad_progeny[:, pw].X, columns=[pw], index=ad_progeny.obs_names)
        .assign(
            cell_type=ad_progeny.obs["cell_type"], organoid=ad_progeny.obs["organoid"]
        )
        .pivot_table(values=pw, columns="cell_type", index="organoid", aggfunc=np.mean)
    )
    tmp_ad.obs = tmp_ad.obs.reset_index()
    tmp_ad.obs["organoid"] = pd.Categorical(tmp_ad.obs["organoid"])

    sc.pl.matrixplot(
        tmp_ad,
        groupby="organoid",
        var_names=tmp_ad.var_names,
        vmin=-2,
        vmax=2,
        cmap="bwr",
        swap_axes=False,
        dendrogram=False,
        title=pw,
        save=f"_progeny_organoid_cell_type_{pw}.pdf",
    )

# %%
sc.set_figure_params(figsize=(3, 3))
sc.pl.umap(
    ad_progeny, color=ad_progeny.var.index, vmin=-2, vmax=2, cmap="coolwarm", ncols=4
)

# %%
sc.set_figure_params(figsize=(3, 3))
sc.pl.embedding(
    ad_progeny,
    basis="X_umap_uncorrected",
    color=ad_progeny.var.index,
    vmin=-2,
    vmax=2,
    cmap="coolwarm",
    ncols=4,
)

# %% [markdown]
# ### More pathways from Uhlitz et al
#  * YAP, LGR5, WNT
#  
# `Lgr5_ISC-Merlos`, `Wnt targets` and `EYR-Gregorieff` (YAP/EGFR) from `Signatures_Single_cells.xlsx`; 
#
# `CORDENONSI_YAP_CONSERVED_SIGNATURE` from YAP_targets.txt. 

# %%
sig_sheet = pd.read_excel("../tables/signatures/Signatures_Single_cells.xlsx").drop(
    0, axis="index"
)

# %%
signatures = {
    "YAP": pd.read_csv("../tables/signatures/YAP_targets.txt", skiprows=2, header=None)[
        0
    ].tolist(),
    "LGR5": sig_sheet["Lgr5_ISC-Merlos"].dropna().tolist(),
    "WNT": [x.upper() for x in sig_sheet["Wnt targets"].dropna().tolist()],
    "YAP_EGFR": sig_sheet["EYR-Gregorieff"].dropna().tolist(),
}

# %%
for key, genes in signatures.items():
    print(key)
    sc.tl.score_genes(adata, signatures[key], score_name=key)

# %%
for pw in signatures:
    tmp_ad = sc.AnnData(
        pd.DataFrame(adata.obs[pw], columns=[pw], index=adata.obs_names)
        .assign(
            cell_type=adata.obs["cell_type"], organoid=adata.obs["organoid"]
        )
        .pivot_table(values=pw, columns="cell_type", index="organoid", aggfunc=np.mean)
    )
    tmp_ad.obs = tmp_ad.obs.reset_index()
    tmp_ad.obs["organoid"] = pd.Categorical(tmp_ad.obs["organoid"])

    sc.pl.matrixplot(
        tmp_ad,
        groupby="organoid",
        var_names=tmp_ad.var_names,
        vmin=-1,
        vmax=1,
        cmap="bwr",
        swap_axes=False,
        dendrogram=False,
        title=pw,
        save=f"_signatures_organoid_cell_type_{pw}.pdf",
    )

# %%
sc.pl.matrixplot(adata, 

# %% [markdown]
# ## Dorothea

# %%
regulons = dorothea.load_regulons(
    [
        "A",
        "B",
    ],  # Which levels of confidence to use (A most confident, E least confident)
    organism="Human",  # If working with mouse, set to Mouse
)

# %%
dorothea.run(
    adata,  # Data to use
    regulons,  # Dorothea network
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # TF with less than 5 targets will be ignored
)

# %%
ad_dorothea = dorothea.extract(adata)

# %%
sc.pl.matrixplot(
    ad_dorothea,
    var_names=ad_dorothea.var.index,
    groupby="organoid",
    vmin=-2,
    vmax=2,
    cmap="bwr",
    swap_axes=True,
    dendrogram=True,
)

# %%
sc.pl.matrixplot(
    ad_dorothea,
    var_names=ad_dorothea.var.index,
    groupby="cell_type",
    vmin=-2,
    vmax=2,
    cmap="bwr",
    swap_axes=True,
    dendrogram=True,
)

# %% [markdown]
# ## Export for cellxgene

# %%
adata_cellxgene = sc.AnnData(
    X=sp.hstack(
        [
            adata.raw.X,
            sp.csr_matrix(np.clip(ad_progeny.X, -2, 2)),
            sp.csr_matrix(np.clip(ad_dorothea.X, -2, 2)),
        ]
    ).tocsc(),
    var=pd.DataFrame(
        index=np.hstack(
            [
                adata.raw.var.index,
                [f"PW:{x}" for x in ad_progeny.var.index],
                [f"TF:{x}" for x in ad_dorothea.var.index],
            ]
        )
    ),
    obs=adata.obs.loc[:, ~adata.obs.columns.str.startswith("_")],
    obsm={
        "X_umap": adata.obsm["X_umap"],
        "X_umap_uncorrected": adata.obsm["X_umap_uncorrected"],
    },
)

# %% [markdown]
# ## Write results

# %%
del adata.uns["rank_genes_groups"]

# %%
adata.write_h5ad("../results/71_scrnaseq_de/adata_cell_type.h5ad")

# %%
adata_cellxgene.write_h5ad("../results/71_scrnaseq_de/adata_cellxgene.h5ad")

# %%
de_res.to_csv("../results/71_scrnaseq_de/de_result.csv")
de_res_cell_type.to_csv("../results/71_scrnaseq_de/de_result_cell_type.csv")
