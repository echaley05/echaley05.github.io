# 06a_build_h5ad.py — run in Jupyter Lab from the same folder as 06_scvi_input/
# Builds ERC_singlets.h5ad from the R export

import numpy as np
import pandas as pd
import scanpy as sc

d = "06_scvi_input"

adata = sc.read_10x_h5(f"{d}/counts.h5")
adata.var_names = adata.var["gene_ids"].astype(str)   # Ensembl IDs (unique)

obs = pd.read_csv(f"{d}/obs.csv", index_col="cell")
var = pd.read_csv(f"{d}/var.csv", index_col="gene_id")

# Same cells and genes, same order
assert obs.index.equals(adata.obs_names), "cell order mismatch"
assert var.index.equals(adata.var_names), "gene order mismatch"

adata.obs = obs
adata.var = var
adata.obs["Sample_ID"] = adata.obs["Sample_ID"].astype("category")
adata.obs["Dx"] = pd.Categorical(adata.obs["Dx"], categories=["Control", "ASD"])

# Checks
print(adata)
print("Raw integer counts:", np.all(adata.X.data == np.round(adata.X.data)))
print("Samples:", adata.obs["Sample_ID"].nunique())
print(adata.obs["Dx"].value_counts())

adata.write_h5ad("ERC_singlets.h5ad", compression="gzip")
