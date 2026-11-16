# 08_cluster.py — neighbors, UMAP, Leiden clustering on the scVI latent space
# Run in Jupyter Lab (rnaseq environment) from the folder containing ERC_scvi.h5ad

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

seed = 42
resolutions = [0.3, 0.5, 0.8, 1.0]
out = "08_cluster"
os.makedirs(out, exist_ok=True)

############################################
# 1. Load
############################################

adata = sc.read_h5ad("ERC_scvi.h5ad")
print(adata)

# Gene labels for all plots: "SYMBOL (Ensembl ID)"
adata.var["label"] = adata.var["gene_name"].astype(str) + " (" + adata.var_names + ")"

# Log-normalized values for plotting only (raw counts stay in X and layers["counts"])
adata.layers["lognorm"] = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4, layer="lognorm")
sc.pp.log1p(adata, layer="lognorm")

adata.uns["Dx_colors"] = ["#4DBBD5", "#E64B35"]   # Control, ASD

############################################
# 2. Neighbors, UMAP, Leiden
############################################

sc.pp.neighbors(adata, use_rep="X_scVI", random_state=seed)
sc.tl.umap(adata, random_state=seed)

keys = []
for r in resolutions:
    key = f"leiden_{r}"
    sc.tl.leiden(adata, resolution=r, key_added=key, flavor="igraph",
                 n_iterations=2, directed=False, random_state=seed)
    keys.append(key)
    print(f"Resolution {r}: {adata.obs[key].nunique()} clusters")

############################################
# 3. Tables
############################################

summary = []
for key in keys:
    tab = pd.crosstab(adata.obs[key], adata.obs["Sample_ID"])
    dx = pd.crosstab(adata.obs[key], adata.obs["Dx"])
    t = pd.DataFrame({
        "nuclei": tab.sum(axis=1),
        "Control": dx["Control"],
        "ASD": dx["ASD"],
        "n_samples": (tab > 0).sum(axis=1),
        "top_sample": tab.idxmax(axis=1),
        "top_sample_pct": (100 * tab.max(axis=1) / tab.sum(axis=1)).round(1),
    })
    t.to_csv(f"{out}/clusters_{key}.csv")
    summary.append(pd.DataFrame({"resolution": key, "n_clusters": [len(t)],
                                 "clusters_top_sample_over_50pct": [(t["top_sample_pct"] > 50).sum()]}))
summary = pd.concat(summary)
summary.to_csv(f"{out}/resolution_summary.csv", index=False)
print(summary.to_string(index=False))

############################################
# 4. Figures
############################################

markers = ["SNAP25", "RBFOX3", "SLC17A7", "SATB2", "GAD1", "GAD2",
           "AQP4", "GFAP", "SLC1A2", "MOBP", "PLP1", "MBP",
           "PDGFRA", "VCAN", "CX3CR1", "P2RY12", "CSF1R", "CLDN5", "FLT1"]
sym = adata.var["gene_name"].astype(str)
marker_labels = []
for g in markers:
    hits = adata.var_names[sym == g]
    if len(hits) == 0:
        print("Marker not found:", g)
        continue
    if len(hits) > 1:
        print(f"{g} matches {len(hits)} IDs; using {hits[0]}")
    marker_labels.append(adata.var.loc[hits[0], "label"])

qc_cols = ["nCount_RNA", "nFeature_RNA", "percent.mt", "scDblFinder.score"]
adata.obs["log10_nCount"] = np.log10(adata.obs["nCount_RNA"])

with PdfPages(f"{out}/cluster_figures.pdf") as pdf:

    # UMAP by cluster at each resolution
    for key in keys:
        fig = sc.pl.umap(adata, color=key, legend_loc="on data", legend_fontsize=7,
                         title=key, frameon=False, show=False, return_fig=True)
        fig.set_size_inches(8, 7)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    # Integration check: sample and Dx
    fig = sc.pl.umap(adata, color="Sample_ID", legend_loc=None, title="Sample_ID",
                     frameon=False, show=False, return_fig=True)
    fig.set_size_inches(8, 7)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig = sc.pl.umap(adata, color="Dx", title="Dx", frameon=False, show=False, return_fig=True)
    fig.set_size_inches(8, 7)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    # Quality check
    fig = sc.pl.umap(adata, color=["log10_nCount", "nFeature_RNA", "percent.mt", "scDblFinder.score"],
                     ncols=2, frameon=False, show=False, return_fig=True)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    for key in keys:
        # Sample makeup of each cluster
        tab = pd.crosstab(adata.obs[key], adata.obs["Sample_ID"], normalize="index")
        fig, ax = plt.subplots(figsize=(12, 5))
        tab.plot(kind="bar", stacked=True, legend=False, ax=ax, width=0.9, colormap="tab20")
        ax.set_title(f"{key}: share of each cluster from each sample")
        ax.set_xlabel("Cluster"); ax.set_ylabel("Fraction of nuclei")
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Dx makeup of each cluster
        dx = pd.crosstab(adata.obs[key], adata.obs["Dx"], normalize="index")[["Control", "ASD"]]
        fig, ax = plt.subplots(figsize=(12, 4))
        dx.plot(kind="bar", stacked=True, ax=ax, width=0.9, color=["#4DBBD5", "#E64B35"])
        ax.axhline(0.5, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(f"{key}: Control vs ASD share of each cluster")
        ax.set_xlabel("Cluster"); ax.set_ylabel("Fraction of nuclei")
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Canonical markers
        dp = sc.pl.dotplot(adata, var_names=marker_labels, groupby=key, gene_symbols="label",
                           layer="lognorm", standard_scale="var", show=False, return_fig=True)
        dp.make_figure()
        dp.fig.suptitle(f"{key}: canonical markers")
        pdf.savefig(dp.fig, bbox_inches="tight"); plt.close(dp.fig)

############################################
# 5. Save
############################################

obs_out = adata.obs[keys].copy()
obs_out[["UMAP_1", "UMAP_2"]] = adata.obsm["X_umap"]
obs_out.to_csv(f"{out}/clusters_and_umap.csv")

with open(f"{out}/cluster_settings.txt", "w") as f:
    f.write(f"Neighbors: use_rep=X_scVI (30 dims), n_neighbors=15 (default), random_state={seed}\n")
    f.write(f"UMAP: scanpy defaults, random_state={seed}\n")
    f.write(f"Leiden: flavor=igraph, n_iterations=2, undirected, random_state={seed}\n")
    for key in keys:
        f.write(f"{key}: {adata.obs[key].nunique()} clusters\n")
    f.write("Plotting values: normalize_total (1e4) + log1p, layer 'lognorm'\n")

adata.write_h5ad("ERC_clustered.h5ad", compression="gzip")
print("Saved: ERC_clustered.h5ad and 08_cluster/ (figures, tables, clusters_and_umap.csv)")
