# Single-nucleus RNA-seq of the entorhinal cortex in autism

Analysis code for single-nucleus RNA sequencing of post-mortem human entorhinal cortex from donors with autism spectrum disorder (ASD) and neurotypical controls.

- **Platform:** Parse Biosciences Evercode WT, single Parse experiment, 46 donors (23 ASD, 23 control, all male)
- **After QC:** 43 donors, 74,057 nuclei, 12 cell types and subtypes
- **Languages:** R (QC, differential expression, enrichment, clinical analyses) and Python (scVI integration, clustering, annotation)

> Data are not included in this repository. See [Data](#data).

## Pipeline

Scripts are numbered in the order they run. Each step reads the previous step's output from the working directory and writes its own output folder.

| Step | Script | Language | What it does |
|---|---|---|---|
| 01 | `01_load_raw.R` | R | Load the raw Parse matrix (`ERC.h5ad`) and donor metadata into Seurat |
| 02 | `02_qc.R` | R | Annotate gene biotypes (GENCODE GTF), compute QC metrics, sample complexity criterion, diagnostic tables and figures |
| 03 | `03_filter_preview.R` | R | Preview sample exclusions and nucleus thresholds by sample and diagnosis |
| 04 | `04_apply_filters.R` | R | Apply sample/nucleus filters; keep genes detected in ≥ 10 nuclei |
| 05 | `05_doublets.R` | R | scDblFinder per sample |
| 06 | `06_export_for_scvi.R`, `06a_build_h5ad.py` | R → Python | Export counts and metadata; build `ERC_singlets.h5ad` |
| 07 | `07_scvi.py` | Python | Variable genes and scVI integration |
| 08 | `08_cluster.py` | Python | Neighbors, UMAP, Leiden clustering |
| 09 | `09a_markers.py`, `09b_annotate.py` | Python | Marker evidence, donor-restriction screen, broad annotation, scANVI check |
| 10 | `10_subcluster.py`, `10b_dotplot.py` | Python | Re-integration and subclustering of excitatory and inhibitory neurons; marker dot plots |
| 11 | `11_subtype_annotate.py`, `11b_umap_figures.py` | Python | Subtype labels; publication UMAPs |
| 12 | `12a_export_for_R.py`, `12b_import_annotations.R` | Python → R | Move labels and embeddings to the Seurat object |
| 14 | `14a_pseudobulk_deseq2.R`, `14b_mast.R`, `14c_volcano.R` | R | ASD vs control differential expression (pseudobulk DESeq2; MAST); volcano plots |
| 15 | `15_GO.R` | R | GO, Reactome and SynGO over-representation and GSEA |
| 16 | `16_risk_genes.R`, `16b_risk_overlap.R` | R | Enrichment of ASD risk genes (SFARI, Satterstrom 2020, Fu 2022) with an expression-matched null; overlap figures |
| 18 | `18_age_by_dx.R`, `18b_age_by_dx_pseudobulk.R` | R | Age × diagnosis interaction (MAST; donor-level DESeq2) |
| 19 | `19_subcluster_MAST.R` | R | MAST within each neuronal subcluster, with per-subcluster breakdowns |
| 20 | `20a_clinical.R`, `20b_clinical_groups.R` | R | Per-donor signatures, patient subgroups, ADI-R and co-occurring conditions |

Step 13 (cell-type composition) is not yet written. `archive/17_hdWGCNA.R` (co-expression modules) was written but not run.

Steps 10 and 11 are run once per neuronal class: set `cell_type` at the top of the script.

## Requirements

- **R ≥ 4.5:** Seurat (v5), SeuratObject, zellkonverter, SingleCellExperiment, scDblFinder, DropletUtils, Matrix, MASS, dplyr, tidyr, ggplot2, ggrepel, MAST, data.table, DESeq2, edgeR, apeglm, clusterProfiler, ReactomePA, org.Hs.eg.db, enrichplot, readxl, pheatmap, cluster
- **Python (conda environment):** scanpy, anndata, scvi-tools, scikit-misc, python-igraph, pandas, numpy, matplotlib

See `environment/` for package lists. Exact versions used are recorded in the `sessionInfo*.txt` and `*_settings.txt` files each step writes.

## Data

The following inputs are needed but are not distributed here:

| File | Used by | Source |
|---|---|---|
| `ERC.h5ad` | 01 | Parse split-pipe output (raw counts) |
| `sample_metadata.tsv` | 01 | Donor metadata (Sample_ID, Dx, Age, PMI, RIN) |
| `gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz` | 02 | GENCODE |
| `clinical_rnaseq_merged.csv` | 20 | Donor clinical data (ADI-R, co-occurring conditions) |
| `risk_genes/` (SFARI CSV, Satterstrom and Fu tables) | 16 | SFARI Gene; supplementary tables of the cited papers |
| `SynGo/ontologies.xlsx` | 15 | SynGO bulk download |

Raw and processed data will be deposited at [repository and accession, once available].

## Reproducibility

- Random seed 42 throughout.
- Settings are at the top of each script; each step writes a settings file recording what was used.

## Citation

[Manuscript citation, once available]
