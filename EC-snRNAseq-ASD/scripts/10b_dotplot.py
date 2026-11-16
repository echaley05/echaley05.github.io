# 10b_dotplot.py — marker dot plot on saved subclusters (no reclustering)
# Run in Jupyter Lab from the same folder as 10_subcluster.py
# Set cell_type, edit the gene list, and rerun as often as you like

import os
import scanpy as sc
import matplotlib.pyplot as plt

############################################
# Settings
############################################

cell_type = "ExcitatoryNeuron"      # or "InhibitoryNeuron"
key = "sub_leiden_0.3"
cluster_rows = True                 # True: order subclusters by similarity and draw a tree

genes = {
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

out = f"10_subcluster_{cell_type}"          # same folder 10_subcluster.py wrote to
os.makedirs(out, exist_ok=True)

############################################
# Load and plot
############################################

adata = sc.read_h5ad(f"ERC_{cell_type}_sub.h5ad")
sym = adata.var["gene_name"].astype(str)

panel = {}
for grp, gl in genes.items():
    found = [g for g in gl if (sym == g).any()]
    missing = [g for g in gl if g not in found]
    if missing:
        print("Not found:", missing)
    if found:
        panel[grp] = [adata.var.loc[adata.var_names[sym == g][0], "label"] for g in found]

if cluster_rows:
    # Tree built from each subcluster's mean position in the subclustering scVI latent space
    sc.tl.dendrogram(adata, groupby=key, use_rep="X_scVI_sub",
                     cor_method="pearson", linkage_method="complete")

dp = sc.pl.dotplot(adata, var_names=panel, groupby=key, gene_symbols="label",
                   layer="lognorm", standard_scale="var", dendrogram=cluster_rows,
                   show=False, return_fig=True)
dp.make_figure()
dp.fig.suptitle(f"{cell_type} — {key}")
dp.fig.savefig(f"{out}/dotplot_{key}.pdf", bbox_inches="tight")
dp.fig.savefig(f"{out}/dotplot_{key}.png", dpi=300, bbox_inches="tight")
plt.show()
print(f"Saved {out}/dotplot_{key}.pdf and .png")
