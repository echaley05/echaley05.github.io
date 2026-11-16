# 07_scvi.py — scVI integration on the filtered singlets
# Run in Jupyter Lab (rnaseq environment) from the folder containing ERC_singlets.h5ad

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt

seed = 42
scvi.settings.seed = seed

############################################
# 1. Load; keep raw counts in a layer
############################################

adata = sc.read_h5ad("ERC_singlets.h5ad")
print(adata)

# Genes were already filtered to >= 10 nuclei in 04_apply_filters.R
adata.layers["counts"] = adata.X.copy()

############################################
# 2. Integration genes: 3,000 variable genes, then drop MT and ribosomal genes
############################################

# Candidates for variable-gene selection: protein-coding or lncRNA, detected in at least 100 nuclei
# (all other biotypes stay in the object for downstream analyses)
hvg_biotypes = ["protein_coding", "lncRNA"]
n_cells = np.asarray((adata.X > 0).sum(axis=0)).ravel()
is_biotype = adata.var["gene_biotype"].isin(hvg_biotypes).values
candidate = (n_cells >= 100) & is_biotype
print("Genes detected in >= 100 nuclei:", int((n_cells >= 100).sum()))
print("Candidate genes (protein_coding + lncRNA, >= 100 nuclei):", int(candidate.sum()),
      "| excluded as candidates:", int((~candidate).sum()))

hvg_input = adata[:, candidate].copy()
sc.pp.highly_variable_genes(hvg_input, n_top_genes=3000, flavor="seurat_v3",
                            layer="counts", batch_key="Sample_ID", span=0.5, subset=False)

adata.var["highly_variable"] = False
adata.var.loc[hvg_input.var_names, "highly_variable"] = hvg_input.var["highly_variable"].values

adata.var["use_for_scvi"] = (adata.var["highly_variable"]
                             & ~adata.var["is_mt"].astype(bool)
                             & ~adata.var["is_rb"].astype(bool))

print("Variable genes selected:", int(adata.var["highly_variable"].sum()))
print("Removed (MT):", int((adata.var["highly_variable"] & adata.var["is_mt"].astype(bool)).sum()))
print("Removed (ribosomal):", int((adata.var["highly_variable"] & adata.var["is_rb"].astype(bool)).sum()))
print("Genes used for scVI:", int(adata.var["use_for_scvi"].sum()))
print(adata.var.loc[adata.var["use_for_scvi"], "gene_biotype"].value_counts())

adata.var.loc[adata.var["use_for_scvi"], ["gene_name", "gene_biotype"]].to_csv("scvi_genes.csv")

adata_scvi = adata[:, adata.var["use_for_scvi"]].copy()

############################################
# 3. Set up and train the model
############################################

scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key="Sample_ID")

model = scvi.model.SCVI(adata_scvi, n_latent=30, gene_likelihood="nb")
print(model)

model.train(max_epochs=400, early_stopping=True, accelerator="auto")

############################################
# 4. Training curve (check that it levelled off)
############################################

hist = model.history
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(hist["elbo_train"], label="train")
ax.plot(hist["elbo_validation"], label="validation")
ax.set_xlabel("Epoch")
ax.set_ylabel("ELBO")
ax.legend()
ax.set_title("scVI training")
fig.savefig("scvi_training_curve.pdf", bbox_inches="tight")
plt.show()
print("Epochs trained:", len(hist["elbo_train"]))

############################################
# 5. Save latent space, model, and object
############################################

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

latent_df = pd.DataFrame(latent, index=adata.obs_names,
                         columns=[f"scVI_{i+1}" for i in range(latent.shape[1])])
latent_df.to_csv("scvi_latent.csv")

model.save("scvi_model", overwrite=True)
adata.write_h5ad("ERC_scvi.h5ad", compression="gzip")

with open("scvi_settings.txt", "w") as f:
    f.write(f"scvi-tools version: {scvi.__version__}\n")
    f.write(f"scanpy version: {sc.__version__}\n")
    f.write(f"Input: ERC_singlets.h5ad ({adata.n_obs} nuclei, {adata.n_vars} genes; all biotypes, genes in >= 10 nuclei)\n")
    f.write("HVG candidates: protein_coding + lncRNA genes detected in >= 100 nuclei (other biotypes kept in the object)\n")
    f.write("HVG: seurat_v3, n_top_genes=3000, batch_key=Sample_ID, raw counts, loess span=0.5\n")
    f.write("MT (37) and RPS/RPL (all ^RPS|^RPL matches) genes removed from the HVG list after selection\n")
    f.write(f"Genes used for scVI: {int(adata.var['use_for_scvi'].sum())}\n")
    f.write("setup_anndata: layer=counts, batch_key=Sample_ID, no covariates\n")
    f.write("Model: n_latent=30, gene_likelihood=nb\n")
    f.write("Training: max_epochs=400, early_stopping=True\n")
    f.write(f"Epochs trained: {len(hist['elbo_train'])}\n")
    f.write(f"Seed: {seed}\n")

print("Saved: scvi_latent.csv, scvi_genes.csv, scvi_model/, ERC_scvi.h5ad, scvi_training_curve.pdf, scvi_settings.txt")
