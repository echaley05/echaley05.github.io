# 11_subtype_annotate.py — apply subtype labels to the neuronal subclusters and merge into the main object
# Run in Jupyter Lab (rnaseq environment) from the folder containing ERC_annotated.h5ad,
# ERC_ExcitatoryNeuron_sub.h5ad and ERC_InhibitoryNeuron_sub.h5ad

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

key = "sub_leiden_0.3"
out = "11_subtypes"
os.makedirs(out, exist_ok=True)

############################################
# Settings: subcluster -> subtype
############################################

mapping = {
    "InhibitoryNeuron": {
        "3": "PAX6+",
        "4": "VIP+", "7": "VIP+",
        "1": "PVALB+", "6": "PVALB+",
        "0": "SST+", "2": "SST+",
        "5": "LAMP5+", "8": "LAMP5+",
    },
    "ExcitatoryNeuron": {
        **{c: "Ex-UL" for c in ["0", "8", "4", "2", "9", "1", "3", "13"]},
        **{c: "Ex-DL" for c in ["12", "5", "11", "10", "6", "7"]},
    },
}

# Plot order and colours (glia keep the colours used in 09b)
colors = {
    "Ex-UL": "#2CA02C", "Ex-DL": "#1B6E1B",
    "PVALB+": "#D62728", "SST+": "#FF7F0E", "VIP+": "#E377C2", "LAMP5+": "#8C564B", "PAX6+": "#7F7F00",
    "Astrocyte": "#FFBE00", "Oligodendrocyte": "#1F77B4", "OPC": "#17BECF",
    "Microglia": "#9467BD", "Endothelial": "#AAAAAA",
}
glia = ["Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Endothelial"]

# Donor-restriction screen (same rules as 09a)
min_nuclei_donor, breadth_min_pct, top2_max_pct, min_donors_per_dx = 3, 25, 50, 3

############################################
# 1. Load and assign subtypes
############################################

adata = sc.read_h5ad("ERC_annotated.h5ad")
adata.obs["subtype"] = adata.obs["cell_type"].astype(str)

for ct, m in mapping.items():
    sub = sc.read_h5ad(f"ERC_{ct}_sub.h5ad")
    cl = sub.obs[key].astype(str)

    unmapped = set(cl.unique()) - set(m)
    unused = set(m) - set(cl.unique())
    assert not unmapped, f"{ct}: subclusters without a label: {unmapped}"
    assert not unused, f"{ct}: mapped clusters not in the object: {unused}"
    assert sub.obs_names.isin(adata.obs_names).all(), f"{ct}: nuclei not found in ERC_annotated.h5ad"

    adata.obs.loc[sub.obs_names, "subtype"] = cl.map(m).values
    adata.obs.loc[sub.obs_names, "sub_cluster_0.3"] = cl.values
    print(f"\n{ct}: {sub.n_obs} nuclei")
    print(pd.crosstab(cl, cl.map(m)).to_string())

assert adata.obs["subtype"].isin(colors).all(), \
    f"Unexpected labels: {set(adata.obs['subtype']) - set(colors)}"
adata.obs["subtype"] = pd.Categorical(adata.obs["subtype"], categories=list(colors))
adata.obs["sub_cluster_0.3"] = adata.obs["sub_cluster_0.3"].fillna("").astype(str)
adata.uns["subtype_colors"] = list(colors.values())

obs = adata.obs
print("\nNuclei per subtype and Dx:")
by_dx = pd.crosstab(obs["subtype"], obs["Dx"])
by_dx["total"] = by_dx.sum(axis=1)
print(by_dx.to_string())

############################################
# 2. Tables and donor-restriction screen
############################################

by_dx.to_csv(f"{out}/subtype_counts_by_Dx.csv")
tab = pd.crosstab(obs["subtype"], obs["Sample_ID"])
tab.to_csv(f"{out}/subtype_counts_by_donor.csv")
obs[["cell_type", "subtype", "sub_cluster_0.3", "Sample_ID", "Dx"]].to_csv(f"{out}/subtype_assignments.csv")

n_donors = tab.shape[1]
donor_dx = obs.groupby("Sample_ID", observed=True)["Dx"].first()
frac = tab / tab.sum(axis=0)
share = frac.div(frac.sum(axis=1), axis=0)
has = tab >= min_nuclei_donor
top2 = share.apply(lambda r: r.sort_values(ascending=False).head(2).sum(), axis=1)
p = share.where(share > 0)

screen = pd.DataFrame({
    "nuclei": tab.sum(axis=1),
    "donors_ge3": has.sum(axis=1),
    "donors_ge3_pct": (100 * has.sum(axis=1) / n_donors).round(1),
    "donors_ge3_Control": has.loc[:, donor_dx.index[donor_dx == "Control"]].sum(axis=1),
    "donors_ge3_ASD": has.loc[:, donor_dx.index[donor_dx == "ASD"]].sum(axis=1),
    "top_donor_corrected": share.idxmax(axis=1),
    "top_donor_corrected_pct": (100 * share.max(axis=1)).round(1),
    "top2_corrected_pct": (100 * top2).round(1),
    "entropy_norm": (-(p * np.log(p)).sum(axis=1) / np.log(n_donors)).round(3),
})
screen["screen"] = np.where(
    (screen["donors_ge3_pct"] < breadth_min_pct) | (screen["top2_corrected_pct"] > top2_max_pct) |
    (screen["donors_ge3_Control"] < min_donors_per_dx) | (screen["donors_ge3_ASD"] < min_donors_per_dx),
    "CHECK IDENTITY", "pass")
screen.to_csv(f"{out}/subtype_donor_screen.csv")
share.round(4).to_csv(f"{out}/subtype_donor_share_corrected.csv")

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
print("\nDonor-restriction screen (report only; nothing excluded):")
print(screen.to_string())

############################################
# 3. Figures
############################################

with PdfPages(f"{out}/subtype_figures.pdf") as pdf:

    note = ("Nuclei from excluded clusters (one-donor population, low-quality nuclei) and nuclei with "
            "scANVI confidence < 0.9 are not shown.")

    # UMAP: neuronal subtypes + glia
    fig = sc.pl.umap(adata, color="subtype", title="Subtypes and glia", frameon=False,
                     show=False, return_fig=True)
    fig.set_size_inches(9, 7)
    fig.text(0.5, 0.02, note, ha="center", va="top", fontsize=7, wrap=True)
    pdf.savefig(fig, bbox_inches="tight")
    fig.savefig(f"{out}/umap_subtypes.png", dpi=300, bbox_inches="tight")
    plt.show(); plt.close(fig)

    # UMAP: subtypes on a grey background, one panel per subtype
    fig, axes = plt.subplots(3, 4, figsize=(16, 11))
    for ax, st in zip(axes.ravel(), colors):
        sc.pl.umap(adata, ax=ax, show=False, frameon=False, size=2, color=None)
        m = (adata.obs["subtype"] == st).values
        xy = adata.obsm["X_umap"]
        ax.scatter(xy[m, 0], xy[m, 1], s=2, c=colors[st])
        ax.set_title(f"{st} (n={int(m.sum())})", fontsize=10)
    fig.text(0.5, 0.06, note, ha="center", va="top", fontsize=8)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    # Per-donor proportions, ordered by Dx
    prop = pd.crosstab(obs["Sample_ID"], obs["subtype"], normalize="index")[list(colors)]
    prop = prop.loc[donor_dx.loc[prop.index].sort_values().index]
    fig, ax = plt.subplots(figsize=(14, 5))
    prop.plot(kind="bar", stacked=True, ax=ax, width=0.9, color=list(colors.values()))
    ax.set_xticklabels([f"{d} ({donor_dx[d]})" for d in prop.index], rotation=90, fontsize=7)
    ax.set_ylabel("Fraction of nuclei"); ax.set_xlabel(None)
    ax.set_title("Subtype proportions per donor")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    # Donor makeup of each subtype
    order = list(donor_dx.sort_values().index)
    fig, ax = plt.subplots(figsize=(14, 0.35 * share.shape[0] + 2))
    im = ax.imshow(share[order].values, aspect="auto", cmap="Reds", vmin=0, vmax=0.5)
    ax.set_yticks(range(share.shape[0]))
    ax.set_yticklabels([f"{c}  ({'CHECK' if s == 'CHECK IDENTITY' else 'ok'})"
                        for c, s in zip(share.index, screen["screen"])], fontsize=8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{d} ({donor_dx[d]})" for d in order], rotation=90, fontsize=6)
    ax.set_title("Size-corrected donor share per subtype (rows sum to 1; colour capped at 0.5)")
    fig.colorbar(im, ax=ax, fraction=0.02)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

############################################
# 4. Save
############################################

adata.write_h5ad("ERC_annotated_subtypes.h5ad", compression="gzip")

with open(f"{out}/subtype_settings.txt", "w") as f:
    f.write(f"Subcluster key: {key}\n")
    for ct, m in mapping.items():
        groups = {}
        for c, s in m.items():
            groups.setdefault(s, []).append(c)
        f.write(f"{ct}: " + "; ".join(f"{s} = clusters {', '.join(cs)}" for s, cs in groups.items()) + "\n")
    f.write("Glia keep their cell-type label as subtype\n")
    f.write("UMAPs show kept nuclei only: excluded clusters (5, 17) and nuclei with scANVI confidence < 0.9 are not plotted\n")
    f.write(f"Donor screen: donor has subtype if >= {min_nuclei_donor} nuclei; flag if < {breadth_min_pct}% of donors, "
            f"top-2 corrected share > {top2_max_pct}%, or < {min_donors_per_dx} donors in either Dx group (report only)\n")
    f.write(f"Nuclei: {adata.n_obs}\n")

print(f"\nSaved: ERC_annotated_subtypes.h5ad ({adata.n_obs} nuclei) and {out}/")
