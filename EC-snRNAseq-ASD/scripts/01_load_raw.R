############################################
# 01_load_raw.R — load ERC.h5ad + donor metadata
# Parse Evercode WT, 46 samples, no filtering
############################################

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(dplyr)
library(tidyr)

############################################
# 1. Metadata
############################################

metadata <- read.delim("sample_metadata.tsv", header = FALSE)
colnames(metadata) <- c("Sample_ID", "Variable", "Value")

metadata_wide <- metadata %>%
  pivot_wider(names_from = Variable, values_from = Value)

numeric_cols <- setdiff(colnames(metadata_wide), c("Sample_ID", "Dx"))
metadata_wide <- metadata_wide %>%
  mutate(across(all_of(numeric_cols), as.numeric),
         Dx = factor(Dx, levels = c("Control", "ASD")))

stopifnot(nrow(metadata_wide) == 46, !anyNA(metadata_wide))

############################################
# 2. Load h5ad
############################################

sce <- readH5AD("ERC.h5ad", reader = "R", X_name = "counts")
counts <- assay(sce, "counts")
rownames(counts) <- rowData(sce)$gene_id
obs <- as.data.frame(colData(sce))

write.csv(as.data.frame(rowData(sce)), "gene_table.csv", row.names = FALSE)

# Find which obs column matches Sample_ID
sample_col <- c("sample", "subject_id")[
  sapply(c("sample", "subject_id"),
         function(x) all(obs[[x]] %in% metadata_wide$Sample_ID))][1]
print(sample_col)

############################################
# 3. Seurat object + metadata
############################################

seu <- CreateSeuratObject(counts = counts, meta.data = obs,
                          min.cells = 0, min.features = 0)
rm(sce, counts); gc()

seu$Sample_ID <- as.character(seu@meta.data[[sample_col]])
seu@meta.data <- seu@meta.data %>%
  left_join(metadata_wide, by = "Sample_ID") %>%
  `rownames<-`(colnames(seu))

table(seu$Sample_ID, seu$Dx)
seu

saveRDS(seu, "01_raw_loaded.rds")
