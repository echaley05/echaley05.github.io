# 09b_annotate.py — apply manual annotation, scANVI for cluster 16, confidence filter
# Run in Jupyter Lab (rnaseq environment) from the folder containing ERC_clustered.h5ad and scvi_model/

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

seed = 42
scvi.settings.seed = seed
key = "leiden_0.3"
conf_min = 0.9
out = "09_annotation"
os.makedirs(out, exist_ok=True)

############################################
# Settings: cluster -> cell type
############################################

mapping = {
    "0": "OPC",
    "1": "Microglia",
    "2": "Astrocyte",
    "3": "ExcitatoryNeuron", "4": "ExcitatoryNeuron", "6": "ExcitatoryNeuron", "7": "ExcitatoryNeuron",
    "9": "ExcitatoryNeuron", "10": "ExcitatoryNeuron", "15": "ExcitatoryNeuron", "20": "ExcitatoryNeuron",
    "8": "InhibitoryNeuron", "13": "InhibitoryNeuron", "14": "InhibitoryNeuron",
    "18": "InhibitoryNeuron", "19": "InhibitoryNeuron",
    "11": "Oligodendrocyte",
    "12": "Endothelial",
    "16": "Unknown",          # assigned per nucleus by scANVI
}
excluded_clusters = {"5": "one-donor population", "17": "low-quality nuclei"}

colors = {"ExcitatoryNeuron": "#2CA02C", "InhibitoryNeuron": "#D62728", "Astrocyte": "#FFBE00",
          "Oligodendrocyte": "#1F77B4", "OPC": "#17BECF", "Microglia": "#9467BD", "Endothelial": "#AAAAAA"}

############################################
# 1. Load and apply manual labels
############################################

adata = sc.read_h5ad("ERC_clustered.h5ad")
cl = adata.obs[key].astype(str)

missing = set(cl.unique()) - set(mapping) - set(excluded_clusters)
assert not missing, f"Clusters without a label: {missing}"

adata.obs["manual_label"] = cl.map(mapping).fillna("Excluded")
adata.obs["exclusion_reason"] = cl.map(excluded_clusters).fillna("")
print(adata.obs["manual_label"].value_counts())

############################################
# 2. scANVI from the trained scVI model
############################################

# Labels for training: manual labels; cluster 16 and excluded clusters are "Unknown" (unlabeled)
adata.obs["scanvi_input"] = adata.obs["manual_label"].replace({"Excluded": "Unknown"})

adata_scvi = adata[:, adata.var["use_for_scvi"].astype(bool)].copy()
print("Genes used (same as scVI):", adata_scvi.n_vars)

vae = scvi.model.SCVI.load("scvi_model", adata=adata_scvi)
lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=adata_scvi, labels_key="scanvi_input",
                                         unlabeled_category="Unknown")
lvae.train(accelerator="auto")
print("scANVI epochs trained:", len(lvae.history["train_loss_epoch"]))

soft = lvae.predict(soft=True)
soft.index = adata.obs_names
if "Unknown" in soft.columns:   # some versions include the unlabeled category; drop and renormalise
    soft = soft.drop(columns="Unknown")
    soft = soft.div(soft.sum(axis=1), axis=0)
adata.obs["scanvi_pred"] = soft.idxmax(axis=1).values
adata.obs["scanvi_conf"] = soft.max(axis=1).values
adata.obsm["X_scANVI"] = lvae.get_latent_representation()

############################################
# 3. Final labels
############################################

obs = adata.obs
final = obs["manual_label"].astype(str).copy()
is16 = cl == "16"
final[is16] = obs.loc[is16, "scanvi_pred"]

status = pd.Series("kept", index=obs.index)
status[obs["manual_label"] == "Excluded"] = "excluded_cluster"
low = (obs["scanvi_conf"] < conf_min) & (status == "kept")
status[low] = "low_confidence"

obs["cell_type"] = np.where(status == "kept", final, "Excluded")
obs["annotation_status"] = status.values

############################################
# 4. Validation tables
############################################

labeled = obs["manual_label"].isin(colors.keys())
agree = (obs.loc[labeled, "manual_label"] == obs.loc[labeled, "scanvi_pred"])
print(f"\nAgreement manual vs scANVI (labeled clusters): {100 * agree.mean():.2f}%")
by_type = agree.groupby(obs.loc[labeled, "manual_label"], observed=True).mean().mul(100).round(2)
print(by_type.to_string())

confusion = pd.crosstab(obs.loc[labeled, "manual_label"], obs.loc[labeled, "scanvi_pred"])
confusion.to_csv(f"{out}/confusion_manual_vs_scanvi.csv")

print("\nCluster 16 assigned by scANVI:")
print(obs.loc[is16, "scanvi_pred"].value_counts().to_string())
print("Cluster 16 nuclei below confidence threshold:", int((obs.loc[is16, "scanvi_conf"] < conf_min).sum()))

print(f"\nNuclei by status (confidence threshold {conf_min}):")
print(pd.crosstab(obs["annotation_status"], obs["Dx"]).to_string())
print("\nLow-confidence nuclei by manual label:")
print(obs.loc[low, "manual_label"].value_counts().to_string())

obs[[key, "manual_label", "scanvi_pred", "scanvi_conf", "cell_type", "annotation_status",
     "exclusion_reason", "Sample_ID", "Dx"]].to_csv(f"{out}/label_assignments.csv")

kept = obs["annotation_status"] == "kept"
by_dx = pd.crosstab(obs.loc[kept, "cell_type"], obs.loc[kept, "Dx"])
by_dx["total"] = by_dx.sum(axis=1)
by_dx.to_csv(f"{out}/celltype_counts_by_Dx.csv")
pd.crosstab(obs.loc[kept, "Sample_ID"], obs.loc[kept, "cell_type"]).to_csv(f"{out}/celltype_counts_by_donor.csv")
print("\nFinal cell types (kept nuclei):")
print(by_dx.to_string())

############################################
# 5. Figures
############################################

ad_kept = adata[kept.values].copy()
ad_kept.obs["cell_type"] = pd.Categorical(ad_kept.obs["cell_type"], categories=list(colors))
ad_kept.uns["cell_type_colors"] = list(colors.values())

with PdfPages(f"{out}/annotation_figures.pdf") as pdf:
    fig = sc.pl.umap(ad_kept, color="cell_type", title="Cell type", frameon=False,
                     show=False, return_fig=True)
    fig.set_size_inches(8, 7)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig = sc.pl.umap(adata, color="scanvi_conf", title="scANVI confidence (all nuclei)",
                     frameon=False, show=False, return_fig=True, vmin=0.5, vmax=1)
    fig.set_size_inches(8, 7)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(obs["scanvi_conf"], bins=50, color="grey")
    ax.axvline(conf_min, color="red", linestyle="--")
    ax.set_xlabel("scANVI confidence (max probability)"); ax.set_ylabel("Nuclei")
    ax.set_title(f"Confidence distribution (threshold {conf_min})")
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    prop = pd.crosstab(obs.loc[kept, "Sample_ID"], obs.loc[kept, "cell_type"], normalize="index")[list(colors)]
    donor_dx = obs.groupby("Sample_ID", observed=True)["Dx"].first()
    prop = prop.loc[donor_dx.loc[prop.index].sort_values().index]
    fig, ax = plt.subplots(figsize=(14, 5))
    prop.plot(kind="bar", stacked=True, ax=ax, width=0.9, color=list(colors.values()))
    ax.set_xticklabels([f"{d} ({donor_dx[d]})" for d in prop.index], rotation=90, fontsize=7)
    ax.set_ylabel("Fraction of nuclei"); ax.set_xlabel(None)
    ax.set_title("Cell-type proportions per donor")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

############################################
# 6. Save
############################################

lvae.save("scanvi_model", overwrite=True)
ad_kept.write_h5ad("ERC_annotated.h5ad", compression="gzip")
adata.obs[["cell_type", "annotation_status", "scanvi_pred", "scanvi_conf"]].to_csv(f"{out}/annotation_all_nuclei.csv")

with open(f"{out}/annotation_settings.txt", "w") as f:
    f.write(f"scvi-tools version: {scvi.__version__}\n")
    f.write(f"Clusters: {key}\n")
    f.write("Manual labels: " + "; ".join(f"{k}={v}" for k, v in mapping.items()) + "\n")
    f.write("Excluded clusters: " + "; ".join(f"{k} ({v})" for k, v in excluded_clusters.items()) + "\n")
    f.write("scANVI: SCANVI.from_scvi_model(scvi_model), labels_key=manual labels, "
            "unlabeled_category='Unknown' (cluster 16 + excluded clusters), default training settings\n")
    f.write(f"scANVI epochs trained: {len(lvae.history['train_loss_epoch'])}\n")
    f.write("Cluster 16: each nucleus assigned its scANVI predicted label\n")
    f.write(f"Confidence filter: nuclei with max predicted probability < {conf_min} excluded\n")
    f.write(f"Agreement manual vs scANVI (labeled clusters): {100 * agree.mean():.2f}%\n")
    f.write(f"Nuclei kept: {int(kept.sum())} of {adata.n_obs}\n")
    f.write(f"Seed: {seed}\n")

print(f"\nSaved: ERC_annotated.h5ad ({int(kept.sum())} nuclei), scanvi_model/, {out}/")
