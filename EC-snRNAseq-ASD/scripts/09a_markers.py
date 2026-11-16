# 09a_markers.py — evidence for annotating leiden_0.3 clusters
# Run in Jupyter Lab (rnaseq environment) from the folder containing ERC_clustered.h5ad

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

key = "leiden_0.3"
out = "09_markers"
os.makedirs(out, exist_ok=True)

############################################
# 1. Load
############################################

adata = sc.read_h5ad("ERC_clustered.h5ad")
print(adata)
print(f"{key}: {adata.obs[key].nunique()} clusters")

sym = adata.var["gene_name"].astype(str)

# Figure labels: gene symbols only (repeated symbols get -1, -2, ... so each label is unique)
adata.var["label"] = ad.utils.make_index_unique(pd.Index(sym)).values

############################################
# 2. Marker test: each cluster vs all others (Wilcoxon, log-normalized values)
############################################

sc.tl.rank_genes_groups(adata, groupby=key, method="wilcoxon", layer="lognorm",
                        use_raw=False, pts=True, key_added="markers")

res = sc.get.rank_genes_groups_df(adata, group=None, key="markers")
res = res.rename(columns={"group": "cluster", "names": "gene_id"})
res.insert(2, "gene_name", sym.loc[res["gene_id"]].values)

top = (res[res["logfoldchanges"] > 0]
       .sort_values(["cluster", "scores"], ascending=[True, False])
       .groupby("cluster", observed=True).head(50))
top.to_csv(f"{out}/top50_markers_per_cluster.csv", index=False)

############################################
# 3. Summary table: one row per cluster
############################################

obs = adata.obs
tab = pd.crosstab(obs[key], obs["Sample_ID"])
dx = pd.crosstab(obs[key], obs["Dx"])
summary = pd.DataFrame({
    "nuclei": tab.sum(axis=1),
    "Control": dx["Control"],
    "ASD": dx["ASD"],
    "n_samples": (tab > 0).sum(axis=1),
    "top_sample": tab.idxmax(axis=1),
    "top_sample_pct": (100 * tab.max(axis=1) / tab.sum(axis=1)).round(1),
    "median_nFeature": obs.groupby(key, observed=True)["nFeature_RNA"].median(),
    "median_nCount": obs.groupby(key, observed=True)["nCount_RNA"].median(),
    "median_percent_mt": obs.groupby(key, observed=True)["percent.mt"].median().round(2),
    "median_doublet_score": obs.groupby(key, observed=True)["scDblFinder.score"].median().round(3),
    "top5_markers": top.groupby("cluster", observed=True)["gene_name"]
                       .apply(lambda g: ", ".join(g.head(5))),
})
summary.index.name = "cluster"

############################################
# 3b. Donor-restriction screen (donor-size corrected)
############################################

min_nuclei_donor = 3      # a donor "has" the cluster if it contributes >= 3 nuclei
breadth_min_pct = 25      # flag if fewer than 25% of donors have the cluster
top2_max_pct = 50         # flag if the top 2 donors hold > 50% of the size-corrected share
min_donors_per_dx = 3     # flag if fewer than 3 donors per Dx group have the cluster

n_donors = tab.shape[1]
donor_dx = obs.groupby("Sample_ID", observed=True)["Dx"].first()

# Fraction of each donor's own nuclei that fall in each cluster, rescaled to sum to 1 per cluster
frac = tab / tab.sum(axis=0)
share = frac.div(frac.sum(axis=1), axis=0)

has = tab >= min_nuclei_donor
top2 = share.apply(lambda r: r.sort_values(ascending=False).head(2).sum(), axis=1)
p = share.where(share > 0)
entropy = -(p * np.log(p)).sum(axis=1) / np.log(n_donors)

screen = pd.DataFrame({
    "donors_ge3": has.sum(axis=1),
    "donors_ge3_pct": (100 * has.sum(axis=1) / n_donors).round(1),
    "donors_ge3_Control": has.loc[:, donor_dx.index[donor_dx == "Control"]].sum(axis=1),
    "donors_ge3_ASD": has.loc[:, donor_dx.index[donor_dx == "ASD"]].sum(axis=1),
    "top_donor_corrected": share.idxmax(axis=1),
    "top_donor_corrected_pct": (100 * share.max(axis=1)).round(1),
    "top2_corrected_pct": (100 * top2).round(1),
    "entropy_norm": entropy.round(3),
})
screen["flag_breadth"] = screen["donors_ge3_pct"] < breadth_min_pct
screen["flag_dominance"] = screen["top2_corrected_pct"] > top2_max_pct
screen["flag_group_coverage"] = ((screen["donors_ge3_Control"] < min_donors_per_dx) |
                                 (screen["donors_ge3_ASD"] < min_donors_per_dx))
screen["screen"] = np.where(screen[["flag_breadth", "flag_dominance", "flag_group_coverage"]].any(axis=1),
                            "CHECK IDENTITY", "pass")

summary = summary.join(screen)
summary.to_csv(f"{out}/cluster_summary.csv")
share.round(4).to_csv(f"{out}/donor_share_corrected_by_cluster.csv")

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
print(summary.drop(columns=["top5_markers"]).to_string())
print("\nTop 5 markers per cluster:")
print(summary["top5_markers"].to_string())
print(f"\nScreen settings: donor has cluster if >= {min_nuclei_donor} nuclei; flag if < {breadth_min_pct}% of donors, "
      f"top-2 corrected share > {top2_max_pct}%, or < {min_donors_per_dx} donors in either Dx group")
print("Clusters to check:", list(summary.index[summary["screen"] == "CHECK IDENTITY"]))

############################################
# 4. Figures
############################################

panel = {
    "Neuron":        ["SNAP25", "RBFOX3", "SYT1"],
    "Excitatory":    ["SLC17A7", "SATB2", "RORB"],
    "Inhibitory":    ["GAD1", "GAD2", "LHX6", "ADARB2"],
    "Astrocyte":     ["AQP4", "GFAP", "SLC1A2"],
    "Oligodendro":   ["MOBP", "PLP1", "MBP"],
    "OPC":           ["PDGFRA", "VCAN"],
    "Microglia":     ["CX3CR1", "P2RY12", "CSF1R"],
    "Endothelial":   ["CLDN5", "FLT1"],
    "Vascular/mural": ["PDGFRB", "DCN"],
    "Immune":        ["PTPRC"],
    "Choroid plexus": ["TTR"],
    "Dentate gyrus": ["PROX1"],
}
panel_labels = {}
for grp, genes in panel.items():
    labs = []
    for g in genes:
        hits = adata.var_names[sym == g]
        if len(hits) == 0:
            print("Marker not found:", g)
            continue
        if len(hits) > 1:
            print(f"{g} matches {len(hits)} IDs; using {hits[0]}")
        labs.append(adata.var.loc[hits[0], "label"])
    if labs:
        panel_labels[grp] = labs

with PdfPages(f"{out}/marker_figures.pdf") as pdf:

    # Reference UMAP
    fig = sc.pl.umap(adata, color=key, legend_loc="on data", legend_fontsize=8,
                     title=key, frameon=False, show=False, return_fig=True)
    fig.set_size_inches(8, 7)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    # Donor makeup of each cluster (size-corrected share), donors grouped by Dx
    order = list(donor_dx.sort_values().index)
    fig, ax = plt.subplots(figsize=(14, 0.35 * share.shape[0] + 2))
    im = ax.imshow(share[order].values, aspect="auto", cmap="Reds", vmin=0, vmax=0.5)
    ax.set_yticks(range(share.shape[0]))
    ax.set_yticklabels([f"{c}  ({'CHECK' if s == 'CHECK IDENTITY' else 'ok'})"
                        for c, s in zip(share.index, summary["screen"])], fontsize=8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{d} ({donor_dx[d]})" for d in order], rotation=90, fontsize=6)
    ax.set_title("Size-corrected donor share per cluster (rows sum to 1; colour capped at 0.5)")
    fig.colorbar(im, ax=ax, fraction=0.02)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    # Canonical and off-target marker panel
    dp = sc.pl.dotplot(adata, var_names=panel_labels, groupby=key, gene_symbols="label",
                       layer="lognorm", standard_scale="var", show=False, return_fig=True)
    dp.make_figure()
    pdf.savefig(dp.fig, bbox_inches="tight"); plt.close(dp.fig)

    # Top 3 data-driven markers per cluster
    dp = sc.pl.rank_genes_groups_dotplot(adata, key="markers", n_genes=3, gene_symbols="label",
                                         layer="lognorm", standard_scale="var", dendrogram=False,
                                         show=False, return_fig=True)
    dp.make_figure()
    pdf.savefig(dp.fig, bbox_inches="tight"); plt.close(dp.fig)

    # Quality metrics per cluster
    for m in ["nFeature_RNA", "nCount_RNA", "percent.mt", "scDblFinder.score"]:
        fig, ax = plt.subplots(figsize=(12, 3.5))
        sc.pl.violin(adata, keys=m, groupby=key, stripplot=False, ax=ax, show=False,
                     log=(m in ["nFeature_RNA", "nCount_RNA"]))
        ax.set_title(f"{m} by cluster")
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

print("Saved: 09_markers/ (cluster_summary.csv, top50_markers_per_cluster.csv, marker_figures.pdf)")
