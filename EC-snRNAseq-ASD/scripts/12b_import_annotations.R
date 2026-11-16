############################################
# 12b_import_annotations.R — add Python labels and coordinates to the Seurat object
# Run in R from the folder containing 05_singlets.rds, gene_table_biotype.csv and 12_export/
# Restart R before running
############################################

library(Seurat)
library(ggplot2)

d <- "12_export"

seu <- readRDS("05_singlets.rds")
if (length(Layers(seu, search = "counts")) > 1) seu <- JoinLayers(seu)
cat("Seurat object:", ncol(seu), "nuclei,", nrow(seu), "genes\n")

labels <- read.csv(file.path(d, "labels.csv"), row.names = "cell", check.names = FALSE)
umap   <- as.matrix(read.csv(file.path(d, "umap.csv"), row.names = "cell"))
latent <- as.matrix(read.csv(file.path(d, "scvi_latent.csv"), row.names = "cell"))

############################################
# 1. Barcode check
############################################

missing <- setdiff(rownames(labels), colnames(seu))
cat("Exported nuclei:", nrow(labels), "| found in R object:", nrow(labels) - length(missing), "\n")
if (length(missing) > 0) stop(length(missing), " exported nuclei not in 05_singlets.rds, e.g. ",
                              paste(head(missing, 3), collapse = ", "))
stopifnot(identical(rownames(umap), rownames(labels)), identical(rownames(latent), rownames(labels)))

############################################
# 2. Keep annotated nuclei only; add metadata
############################################

seu <- subset(seu, cells = rownames(labels))
seu <- seu[, rownames(labels)]                      # same order as the export
cat("Kept:", ncol(seu), "nuclei (excluded clusters and low-confidence nuclei dropped)\n")

stopifnot(all(seu$Sample_ID == labels$Sample_ID), all(seu$Dx == labels$Dx))
labels$Sample_ID <- NULL
labels$Dx <- NULL
seu <- AddMetaData(seu, labels)

subtype_levels <- c("Ex-UL", "Ex-DL", "PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+",
                    "Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Endothelial")
celltype_levels <- c("ExcitatoryNeuron", "InhibitoryNeuron", "Astrocyte", "Oligodendrocyte",
                     "OPC", "Microglia", "Endothelial")
stopifnot(all(seu$subtype %in% subtype_levels), all(seu$cell_type %in% celltype_levels))
seu$subtype   <- factor(seu$subtype, levels = subtype_levels)
seu$cell_type <- factor(seu$cell_type, levels = celltype_levels)
seu$Dx        <- factor(seu$Dx, levels = c("Control", "ASD"))
seu$Sample_ID <- droplevels(factor(seu$Sample_ID))

############################################
# 3. Reductions (Python coordinates)
############################################

colnames(umap)   <- paste0("UMAP_", 1:2)
colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
seu[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = DefaultAssay(seu))
seu[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scVI_", assay = DefaultAssay(seu))

Idents(seu) <- "subtype"

############################################
# 4. Normalize (for plotting; DE/composition methods use counts their own way)
############################################

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

############################################
# 5. Checks against Python
############################################

check_py <- read.csv(file.path(d, "check_counts_subtype_by_Dx.csv"), row.names = 1, check.names = FALSE)
check_r  <- as.data.frame.matrix(table(seu$subtype, seu$Dx))
check_r  <- check_r[rownames(check_py), colnames(check_py)]
cat("\nNuclei per subtype x Dx (R):\n"); print(check_r)
cat("Matches Python:", isTRUE(all(check_r == check_py)), "\n")
stopifnot(all(check_r == check_py))

pal <- c("Ex-UL" = "#2CA02C", "Ex-DL" = "#1B6E1B",
         "PVALB+" = "#D62728", "SST+" = "#FF7F0E", "VIP+" = "#E377C2", "LAMP5+" = "#8C564B", "PAX6+" = "#7F7F00",
         "Astrocyte" = "#FFBE00", "Oligodendrocyte" = "#1F77B4", "OPC" = "#17BECF",
         "Microglia" = "#9467BD", "Endothelial" = "#AAAAAA")
seu@misc$subtype_colors <- pal

# Marker check: rownames are Ensembl IDs, so map symbols -> IDs and label axes with symbols
gene_table <- read.csv("gene_table_biotype.csv")
markers <- c("SLC17A7", "SATB2", "CUX2", "FEZF2", "GAD1", "PVALB", "SST", "VIP", "LAMP5", "PAX6",
             "AQP4", "MOBP", "PDGFRA", "CX3CR1", "CLDN5")
ids <- gene_table$gene_id[match(markers, gene_table$gene_name)]
keep <- !is.na(ids) & ids %in% rownames(seu)
if (any(!keep)) cat("Markers not found:", markers[!keep], "\n")
ids <- ids[keep]; names(ids) <- markers[keep]

pdf(file.path(d, "import_check.pdf"), width = 10, height = 7)
print(DimPlot(seu, reduction = "umap", group.by = "subtype", cols = pal, raster = TRUE,
              label = TRUE, repel = TRUE) + ggtitle("Subtype (should match Python UMAP)"))
print(DotPlot(seu, features = unname(ids), group.by = "subtype") +
        scale_x_discrete(labels = names(ids)) + RotatedAxis() + ggtitle("Marker check"))
dev.off()

############################################
# 6. Save
############################################

saveRDS(seu, "12_annotated.rds")
writeLines(capture.output(sessionInfo()), file.path(d, "sessionInfo_12b.txt"))
cat("\nSaved 12_annotated.rds (", ncol(seu), " nuclei, ", nrow(seu), " genes) and ", d, "/import_check.pdf\n", sep = "")
