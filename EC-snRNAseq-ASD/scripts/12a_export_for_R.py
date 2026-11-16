# 12a_export_for_R.py — export labels and coordinates for the R/Seurat object
# Run in Jupyter Lab (rnaseq environment) after 11_subtype_annotate.py
# Counts are NOT exported: R already has them in 05_singlets.rds (same barcodes)

import os
import numpy as np
import pandas as pd
import scanpy as sc

out = "12_export"
os.makedirs(out, exist_ok=True)

adata = sc.read_h5ad("ERC_annotated_subtypes.h5ad", backed="r")
obs = adata.obs
print(f"Kept nuclei: {adata.n_obs}")

############################################
# 1. Labels (one row per kept nucleus)
############################################

labels = obs[["Sample_ID", "Dx", "leiden_0.3", "cell_type", "subtype",
              "sub_cluster_0.3", "scanvi_conf"]].copy()
for c in ["leiden_0.3", "cell_type", "subtype", "sub_cluster_0.3"]:
    labels[c] = labels[c].astype(str)

# Coordinates from each neuron class's own subclustering UMAP (blank for other cell types)
for ct, tag in [("ExcitatoryNeuron", "Exc"), ("InhibitoryNeuron", "Inh")]:
    sub = sc.read_h5ad(f"ERC_{ct}_sub.h5ad", backed="r")
    xy = pd.DataFrame(np.asarray(sub.obsm["X_umap"]), index=sub.obs_names,
                      columns=[f"UMAPsub{tag}_1", f"UMAPsub{tag}_2"])
    sub.file.close()
    assert xy.index.isin(labels.index).all(), f"{ct}: subcluster nuclei missing from the annotated object"
    labels = labels.join(xy)
    print(f"{ct}: {len(xy)} nuclei with subcluster UMAP")

labels.index.name = "cell"
labels.to_csv(f"{out}/labels.csv")

############################################
# 2. Main UMAP and scVI latent space
############################################

umap = pd.DataFrame(np.asarray(adata.obsm["X_umap"]), index=adata.obs_names,
                    columns=["UMAP_1", "UMAP_2"])
umap.index.name = "cell"
umap.to_csv(f"{out}/umap.csv")

lat = np.asarray(adata.obsm["X_scVI"])
latent = pd.DataFrame(lat, index=adata.obs_names, columns=[f"scVI_{i+1}" for i in range(lat.shape[1])])
latent.index.name = "cell"
latent.to_csv(f"{out}/scvi_latent.csv")
adata.file.close()

############################################
# 3. Counts for checking in R
############################################

check = pd.crosstab(labels["subtype"], labels["Dx"])
check.to_csv(f"{out}/check_counts_subtype_by_Dx.csv")
print(check.to_string())

with open(f"{out}/export_settings.txt", "w") as f:
    f.write("Source: ERC_annotated_subtypes.h5ad (kept nuclei only; excluded clusters and scANVI confidence < 0.9 removed)\n")
    f.write(f"Nuclei: {len(labels)}\n")
    f.write("labels.csv: Sample_ID, Dx, leiden_0.3, cell_type, subtype, sub_cluster_0.3, scanvi_conf, "
            "subcluster UMAPs (UMAPsubExc_*, UMAPsubInh_*; blank for other cell types)\n")
    f.write("umap.csv: main UMAP (step 08); scvi_latent.csv: main scVI latent (step 07, "
            f"{lat.shape[1]} dims)\n")

print(f"\nSaved {out}/ (labels.csv, umap.csv, scvi_latent.csv, check_counts_subtype_by_Dx.csv)")
