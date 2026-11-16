# 10_subcluster.py — subcluster one broad cell type with its own scVI model
# Run in Jupyter Lab (rnaseq environment) from the folder containing ERC_annotated.h5ad
# Run once per cell type: set cell_type below

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

############################################
# Settings
############################################

cell_type = "ExcitatoryNeuron"      # or "InhibitoryNeuron"

seed = 42
n_hvg = 2000
min_nuclei_hvg = 10                 # HVG candidates: detected in >= 10 nuclei of this cell type
hvg_biotypes = ["protein_coding", "lncRNA"]
n_latent = 20
n_neighbors = 30
resolutions = [0.2, 0.3, 0.5, 0.8]
marker_res = 0.3                    # resolution used for the marker test; change and rerun section 5 if needed

# Donor screen (same rules as 09a)
min_nuclei_donor, breadth_min_pct, top2_max_pct, min_donors_per_dx = 3, 25, 50, 3

panels = {
    "ExcitatoryNeuron": {
        "Pan-excitatory": ["SLC17A7", "SATB2"],
        "Layer II":       ["RELN", "CALB1"],
        "Upper layers":   ["CUX2", "CBLN2", "LAMP5"],
        "L3/L5 boundary": ["RORB"],
        "Layer V":        ["BCL11B", "FEZF2", "TSHZ2", "ETV1"],
        "Layer VI":       ["FOXP2", "TLE4", "SEMA3E"],
        "Layer VIb":      ["CCN2", "NR4A2"],
    },
    "InhibitoryNeuron": {
        "Pan-inhibitory": ["GAD1", "GAD2"],
        "MGE":            ["LHX6", "PVALB", "SST", "CHODL", "NOS1"],
        "CGE":            ["ADARB2", "VIP", "LAMP5", "PAX6", "SNCG", "NDNF", "RELN"],
        "Other":          ["CXCL14", "NPY", "CCK", "CALB2", "CA8"],
    },
}[cell_type]

scvi.settings.seed = seed
out = f"10_subcluster_{cell_type}"
os.makedirs(out, exist_ok=True)

############################################
# 1. Subset
############################################

full = sc.read_h5ad("ERC_annotated.h5ad")
adata = full[full.obs["cell_type"] == cell_type].copy()
del full
adata.obs["Sample_ID"] = adata.obs["Sample_ID"].cat.remove_unused_categories()
print(adata)
print("Nuclei per Dx:", adata.obs["Dx"].value_counts().to_dict())
print("Donors:", adata.obs["Sample_ID"].nunique(),
      "| nuclei per donor min/median/max:",
      adata.obs["Sample_ID"].value_counts().agg(["min", "median", "max"]).tolist())

# Figure labels: gene symbols (repeated symbols get -1, -2, ...)
adata.var["label"] = ad.utils.make_index_unique(pd.Index(adata.var["gene_name"].astype(str))).values

############################################
# 2. Variable genes within this cell type
############################################

n_cells = np.asarray((adata.layers["counts"] > 0).sum(axis=0)).ravel()
candidate = (n_cells >= min_nuclei_hvg) & adata.var["gene_biotype"].isin(hvg_biotypes).values
print(f"HVG candidates (protein_coding + lncRNA, >= {min_nuclei_hvg} nuclei):", int(candidate.sum()))

hvg_input = adata[:, candidate].copy()
span_used = 0.5
try:
    sc.pp.highly_variable_genes(hvg_input, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts",
                                batch_key="Sample_ID", span=span_used, subset=False)
except ValueError as e:
    print("HVG loess failed at span 0.5 (", e, "); retrying at span 0.8")
    span_used = 0.8
    sc.pp.highly_variable_genes(hvg_input, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts",
                                batch_key="Sample_ID", span=span_used, subset=False)

adata.var["highly_variable_sub"] = False
adata.var.loc[hvg_input.var_names, "highly_variable_sub"] = hvg_input.var["highly_variable"].values
adata.var["use_for_scvi_sub"] = (adata.var["highly_variable_sub"]
                                 & ~adata.var["is_mt"].astype(bool) & ~adata.var["is_rb"].astype(bool))
print("HVGs:", int(adata.var["highly_variable_sub"].sum()),
      "| removed MT:", int((adata.var["highly_variable_sub"] & adata.var["is_mt"].astype(bool)).sum()),
      "| removed ribosomal:", int((adata.var["highly_variable_sub"] & adata.var["is_rb"].astype(bool)).sum()),
      "| used for scVI:", int(adata.var["use_for_scvi_sub"].sum()))
adata.var.loc[adata.var["use_for_scvi_sub"], ["gene_name", "gene_biotype"]].to_csv(f"{out}/scvi_genes.csv")

############################################
# 3. scVI on this cell type
############################################

ad_s = adata[:, adata.var["use_for_scvi_sub"]].copy()
scvi.model.SCVI.setup_anndata(ad_s, layer="counts", batch_key="Sample_ID")
model = scvi.model.SCVI(ad_s, n_latent=n_latent, gene_likelihood="nb")
model.train(max_epochs=400, early_stopping=True, accelerator="auto")
hist = model.history
print("Epochs trained:", len(hist["elbo_train"]))
adata.obsm["X_scVI_sub"] = model.get_latent_representation()
model.save(f"{out}/scvi_model", overwrite=True)

############################################
# 4. Neighbors, UMAP, Leiden
############################################

sc.pp.neighbors(adata, use_rep="X_scVI_sub", n_neighbors=n_neighbors, random_state=seed)
sc.tl.umap(adata, random_state=seed)
keys = []
for r in resolutions:
    k = f"sub_leiden_{r}"
    sc.tl.leiden(adata, resolution=r, key_added=k, flavor="igraph", n_iterations=2,
                 directed=False, random_state=seed)
    keys.append(k)
    print(f"Resolution {r}: {adata.obs[k].nunique()} subclusters")

############################################
# 5. Summary + donor screen per resolution; markers at marker_res
############################################

obs = adata.obs
donor_dx = obs.groupby("Sample_ID", observed=True)["Dx"].first()

def summarize(k):
    tab = pd.crosstab(obs[k], obs["Sample_ID"])
    dx = pd.crosstab(obs[k], obs["Dx"])
    frac = tab / tab.sum(axis=0)
    share = frac.div(frac.sum(axis=1), axis=0)
    has = tab >= min_nuclei_donor
    p = share.where(share > 0)
    s = pd.DataFrame({
        "nuclei": tab.sum(axis=1), "Control": dx.get("Control", 0), "ASD": dx.get("ASD", 0),
        "median_nFeature": obs.groupby(k, observed=True)["nFeature_RNA"].median(),
        "median_percent_mt": obs.groupby(k, observed=True)["percent.mt"].median().round(2),
        "median_doublet_score": obs.groupby(k, observed=True)["scDblFinder.score"].median().round(3),
        "donors_ge3": has.sum(axis=1),
        "donors_ge3_pct": (100 * has.sum(axis=1) / tab.shape[1]).round(1),
        "donors_ge3_Control": has.loc[:, donor_dx.index[donor_dx == "Control"]].sum(axis=1),
        "donors_ge3_ASD": has.loc[:, donor_dx.index[donor_dx == "ASD"]].sum(axis=1),
        "top2_corrected_pct": (100 * share.apply(lambda r: r.nlargest(2).sum(), axis=1)).round(1),
        "entropy_norm": (-(p * np.log(p)).sum(axis=1) / np.log(tab.shape[1])).round(3),
    })
    flag = ((s["donors_ge3_pct"] < breadth_min_pct) | (s["top2_corrected_pct"] > top2_max_pct) |
            (s["donors_ge3_Control"] < min_donors_per_dx) | (s["donors_ge3_ASD"] < min_donors_per_dx))
    s["screen"] = np.where(flag, "CHECK IDENTITY", "pass")
    s.index.name = "subcluster"
    return s, share

summaries = {}
for k in keys:
    s, share = summarize(k)
    s.to_csv(f"{out}/summary_{k}.csv")
    summaries[k] = (s, share)
    print(f"{k}: {len(s)} subclusters, {int((s['screen'] != 'pass').sum())} flagged")

mk = f"sub_leiden_{marker_res}"
sc.tl.rank_genes_groups(adata, groupby=mk, method="wilcoxon", layer="lognorm",
                        use_raw=False, pts=True, key_added="sub_markers")
res = sc.get.rank_genes_groups_df(adata, group=None, key="sub_markers").rename(
    columns={"group": "subcluster", "names": "gene_id"})
res.insert(2, "gene_name", adata.var.loc[res["gene_id"], "gene_name"].values)
top = (res[res["logfoldchanges"] > 0].sort_values(["subcluster", "scores"], ascending=[True, False])
       .groupby("subcluster", observed=True).head(50))
top.to_csv(f"{out}/top50_markers_{mk}.csv", index=False)

s = summaries[mk][0].copy()
s["top5_markers"] = top.groupby("subcluster", observed=True)["gene_name"].apply(lambda g: ", ".join(g.head(5)))
s.to_csv(f"{out}/summary_{mk}.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
print(f"\nSummary at {mk}:")
print(s.drop(columns="top5_markers").to_string())
print("\nTop 5 markers:")
print(s["top5_markers"].to_string())

############################################
# 6. Figures
############################################

sym = adata.var["gene_name"].astype(str)
panel_labels = {}
for grp, genes in panels.items():
    labs = []
    for g in genes:
        hits = adata.var_names[sym == g]
        if len(hits) == 0:
            print("Marker not found:", g); continue
        labs.append(adata.var.loc[hits[0], "label"])
    if labs:
        panel_labels[grp] = labs

with PdfPages(f"{out}/subcluster_figures.pdf") as pdf:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(hist["elbo_train"], label="train"); ax.plot(hist["elbo_validation"], label="validation")
    ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend(); ax.set_title(f"scVI training ({cell_type})")
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    for k in keys:
        fig = sc.pl.umap(adata, color=k, legend_loc="on data", legend_fontsize=8, title=k,
                         frameon=False, show=False, return_fig=True)
        fig.set_size_inches(8, 7); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig = sc.pl.umap(adata, color=["Sample_ID", "Dx"], legend_loc=None, frameon=False,
                     show=False, return_fig=True)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
    fig = sc.pl.umap(adata, color=["nFeature_RNA", "percent.mt", "scDblFinder.score", "scanvi_conf"],
                     ncols=2, frameon=False, show=False, return_fig=True)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    for k in keys:
        dp = sc.pl.dotplot(adata, var_names=panel_labels, groupby=k, gene_symbols="label",
                           layer="lognorm", standard_scale="var", show=False, return_fig=True)
        dp.make_figure(); dp.fig.suptitle(f"{k}: marker panel")
        pdf.savefig(dp.fig, bbox_inches="tight"); plt.close(dp.fig)

    dp = sc.pl.rank_genes_groups_dotplot(adata, key="sub_markers", n_genes=3, gene_symbols="label",
                                         layer="lognorm", standard_scale="var", dendrogram=False,
                                         show=False, return_fig=True)
    dp.make_figure(); dp.fig.suptitle(f"{mk}: top 3 markers per subcluster")
    pdf.savefig(dp.fig, bbox_inches="tight"); plt.close(dp.fig)

    s_mk, share = summaries[mk]
    order = list(donor_dx.sort_values().index)
    fig, ax = plt.subplots(figsize=(14, 0.35 * share.shape[0] + 2))
    im = ax.imshow(share[order].values, aspect="auto", cmap="Reds", vmin=0, vmax=0.5)
    ax.set_yticks(range(share.shape[0]))
    ax.set_yticklabels([f"{c} ({'CHECK' if v != 'pass' else 'ok'})" for c, v in zip(share.index, s_mk["screen"])], fontsize=8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{d} ({donor_dx[d]})" for d in order], rotation=90, fontsize=6)
    ax.set_title(f"{mk}: size-corrected donor share"); fig.colorbar(im, ax=ax, fraction=0.02)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    for m in ["nFeature_RNA", "percent.mt", "scDblFinder.score"]:
        fig, ax = plt.subplots(figsize=(12, 3.5))
        sc.pl.violin(adata, keys=m, groupby=mk, stripplot=False, ax=ax, show=False,
                     log=(m == "nFeature_RNA"))
        ax.set_title(f"{m} by subcluster ({mk})")
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

############################################
# 7. Save
############################################

adata.obs[keys].assign(UMAP_sub_1=adata.obsm["X_umap"][:, 0], UMAP_sub_2=adata.obsm["X_umap"][:, 1]) \
    .to_csv(f"{out}/subclusters_and_umap.csv")
adata.write_h5ad(f"ERC_{cell_type}_sub.h5ad", compression="gzip")

with open(f"{out}/subcluster_settings.txt", "w") as f:
    f.write(f"Cell type: {cell_type} ({adata.n_obs} nuclei, {adata.obs['Sample_ID'].nunique()} donors)\n")
    f.write(f"scvi-tools {scvi.__version__}; scanpy {sc.__version__}\n")
    f.write(f"HVGs: seurat_v3, n_top_genes={n_hvg}, batch_key=Sample_ID, span={span_used}, raw counts; "
            f"candidates protein_coding + lncRNA detected in >= {min_nuclei_hvg} nuclei; MT and RPS/RPL removed after\n")
    f.write(f"Genes used for scVI: {int(adata.var['use_for_scvi_sub'].sum())}\n")
    f.write(f"scVI: n_latent={n_latent}, gene_likelihood=nb, other model settings default; "
            f"max_epochs=400, early stopping; epochs trained {len(hist['elbo_train'])}\n")
    f.write(f"Neighbors: X_scVI_sub, n_neighbors={n_neighbors}; UMAP default; Leiden igraph, n_iterations=2, "
            f"resolutions {resolutions}; seed {seed}\n")
    f.write(f"Marker test at {mk}: Wilcoxon on log-normalized values vs rest\n")
print(f"\nSaved: ERC_{cell_type}_sub.h5ad and {out}/")
