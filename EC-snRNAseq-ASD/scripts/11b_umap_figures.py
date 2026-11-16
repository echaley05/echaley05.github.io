# 11b_umap_figures.py — publication-style UMAP versions from ERC_annotated_subtypes.h5ad
# Run in Jupyter Lab (rnaseq environment) after 11_subtype_annotate.py
# Kept nuclei only: excluded clusters and nuclei with scANVI confidence < 0.9 are not shown

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

############################################
# Settings
############################################

seed = 42
out = "11_subtypes/umap_versions"
os.makedirs(out, exist_ok=True)

widths = {"single": 3.5, "double": 7.2}     # inches (typical journal column widths)
point_size = {"single": 0.15, "double": 0.6}
dpi = 600

groups = {
    "Excitatory": ["Ex-UL", "Ex-DL"],
    "Inhibitory": ["PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+"],
    "Glia":       ["Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Endothelial"],
}
colors = {
    "Ex-UL": "#2CA02C", "Ex-DL": "#1B6E1B",
    "PVALB+": "#D62728", "SST+": "#FF7F0E", "VIP+": "#E377C2", "LAMP5+": "#8C564B", "PAX6+": "#7F7F00",
    "Astrocyte": "#FFBE00", "Oligodendrocyte": "#1F77B4", "OPC": "#17BECF",
    "Microglia": "#9467BD", "Endothelial": "#AAAAAA",
}
glia_grey = "#D9D9D9"
dx_order = ["Control", "ASD"]
note = "Excluded clusters and nuclei with scANVI confidence < 0.9 not shown."

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "legend.fontsize": 7,
    "pdf.fonttype": 42, "ps.fonttype": 42,      # editable text in Illustrator
    "savefig.dpi": dpi,
})

############################################
# Load (coordinates and labels only)
############################################

adata = sc.read_h5ad("ERC_annotated_subtypes.h5ad", backed="r")
df = pd.DataFrame(np.asarray(adata.obsm["X_umap"]), index=adata.obs_names, columns=["x", "y"])
df["subtype"] = adata.obs["subtype"].astype(str).values
df["Dx"] = adata.obs["Dx"].astype(str).values
adata.file.close()

assert df["subtype"].isin(colors).all(), f"Unexpected labels: {set(df['subtype']) - set(colors)}"
rng = np.random.default_rng(seed)
df = df.iloc[rng.permutation(len(df))]          # random draw order
print(f"{len(df)} nuclei")
print(df["subtype"].value_counts().reindex(list(colors)).to_string())

pad = 0.03 * (df[["x", "y"]].max() - df[["x", "y"]].min())
xlim = (df["x"].min() - pad["x"], df["x"].max() + pad["x"])
ylim = (df["y"].min() - pad["y"], df["y"].max() + pad["y"])

############################################
# Helpers
############################################

def clean(ax, arrows=True):
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_visible(False)
    if arrows:   # small UMAP1/UMAP2 axis arrows, bottom-left
        kw = dict(xycoords="axes fraction", textcoords="axes fraction",
                  arrowprops=dict(arrowstyle="-|>", lw=0.6, color="black", mutation_scale=5))
        ax.annotate("", xy=(0.14, 0.0), xytext=(0.0, 0.0), **kw)
        ax.annotate("", xy=(0.0, 0.14), xytext=(0.0, 0.0), **kw)
        ax.text(0.07, -0.025, "UMAP1", transform=ax.transAxes, ha="center", va="top", fontsize=6)
        ax.text(-0.025, 0.07, "UMAP2", transform=ax.transAxes, ha="right", va="center",
                rotation=90, fontsize=6)

def scatter(ax, d, size, cmap=colors):
    ax.scatter(d["x"], d["y"], s=size, c=d["subtype"].map(cmap).values,
               linewidths=0, rasterized=True)

def label_centres(ax, d, which, fontsize=7):
    for st in which:
        s = d[d["subtype"] == st]
        if len(s) == 0:
            continue
        ax.text(s["x"].median(), s["y"].median(), st, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color="black",
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])

def grouped_legend(ax, which_groups, **kw):
    handles, labels = [], []
    for g in which_groups:
        handles.append(Line2D([], [], linestyle="none")); labels.append(g)
        for st in groups[g]:
            handles.append(Line2D([], [], marker="o", linestyle="none", markersize=4,
                                  markerfacecolor=colors[st], markeredgewidth=0))
            labels.append(st)
    leg = ax.legend(handles, labels, frameon=False, handletextpad=0.2, borderaxespad=0,
                    labelspacing=0.3, **kw)
    for t in leg.get_texts():
        if t.get_text() in which_groups:
            t.set_fontweight("bold")
    return leg

def save(fig, name):
    fig.savefig(f"{out}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{out}/{name}.png", dpi=dpi, bbox_inches="tight")
    plt.show(); plt.close(fig)
    print("Saved", name)

all_subtypes = list(colors)
neurons = groups["Excitatory"] + groups["Inhibitory"]

############################################
# V1. Labelled clusters, no legend
############################################

for w, W in widths.items():
    fig, ax = plt.subplots(figsize=(W, W))
    scatter(ax, df, point_size[w])
    label_centres(ax, df, all_subtypes, fontsize=6 if w == "single" else 8)
    clean(ax)
    save(fig, f"V1_labelled_{w}")

############################################
# V2. Grouped legend (Excitatory / Inhibitory / Glia)
############################################

for w, W in widths.items():
    fig, ax = plt.subplots(figsize=(W, W * 0.78))
    scatter(ax, df, point_size[w])
    clean(ax)
    grouped_legend(ax, list(groups), loc="center left", bbox_to_anchor=(1.01, 0.5))
    save(fig, f"V2_legend_{w}")

############################################
# V3. Neuron focus: glia in grey
############################################

cmap_focus = {**colors, **{g: glia_grey for g in groups["Glia"]}}
d_focus = pd.concat([df[df["subtype"].isin(groups["Glia"])],      # glia underneath
                     df[df["subtype"].isin(neurons)]])
for w, W in widths.items():
    fig, ax = plt.subplots(figsize=(W, W * 0.78))
    scatter(ax, d_focus, point_size[w], cmap=cmap_focus)
    clean(ax)
    leg = grouped_legend(ax, ["Excitatory", "Inhibitory"], loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax.text(1.03, 0.08, "Glia (grey)", transform=ax.transAxes, fontsize=6, color="#808080")
    save(fig, f"V3_neuron_focus_{w}")

############################################
# V4. Split by diagnosis, downsampled to equal nuclei (double width only)
############################################

n_min = df["Dx"].value_counts().reindex(dx_order).min()
fig, axes = plt.subplots(1, 2, figsize=(widths["double"], widths["double"] * 0.42))
for ax, dx in zip(axes, dx_order):
    d = df[df["Dx"] == dx]
    d = d.iloc[rng.choice(len(d), size=n_min, replace=False)]
    scatter(ax, d, point_size["double"] * 0.8)
    clean(ax, arrows=(dx == dx_order[0]))
    ax.set_title(f"{dx} (n = {n_min:,} nuclei shown)")
grouped_legend(axes[-1], list(groups), loc="center left", bbox_to_anchor=(1.01, 0.5))
fig.text(0.5, -0.02, f"Each group randomly downsampled to {n_min:,} nuclei (seed {seed}).",
         ha="center", fontsize=6)
save(fig, "V4_split_Dx_double")

############################################
# V5. Small multiples: one panel per subtype (double width only)
############################################

ncol = 4
nrow = int(np.ceil(len(all_subtypes) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(widths["double"], widths["double"] * nrow / ncol * 0.95))
for ax, st in zip(axes.ravel(), all_subtypes):
    ax.scatter(df["x"], df["y"], s=0.1, c="#E5E5E5", linewidths=0, rasterized=True)
    s = df[df["subtype"] == st]
    ax.scatter(s["x"], s["y"], s=0.3, c=colors[st], linewidths=0, rasterized=True)
    clean(ax, arrows=False)
    ax.set_title(f"{st} (n = {len(s):,})", fontsize=7)
for ax in axes.ravel()[len(all_subtypes):]:
    ax.axis("off")
save(fig, "V5_small_multiples_double")

with open(f"{out}/README.txt", "w") as f:
    f.write("UMAP figure versions (from ERC_annotated_subtypes.h5ad; coordinates from step 08)\n")
    f.write(note + "\n")
    f.write(f"Points drawn in random order (seed {seed}); points rasterized at {dpi} dpi, text editable (Type 42)\n")
    f.write(f"Widths: single {widths['single']} in, double {widths['double']} in\n")
    f.write(f"V4: each Dx group downsampled to {n_min} nuclei\n")

print(f"\nAll versions saved in {out}/")
