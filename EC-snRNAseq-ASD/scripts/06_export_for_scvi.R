############################################
# 06_export_for_scvi.R — export singlets for scVI in Python
# Counts -> 10x-style HDF5; cell and gene info -> CSV
############################################

library(Seurat)
library(DropletUtils)

out <- "06_scvi_input"
dir.create(out, showWarnings = FALSE)

seu <- readRDS("05_singlets.rds")
gene_table <- read.csv("gene_table_biotype.csv")

counts <- LayerData(seu, layer = "counts")
g <- gene_table[match(rownames(counts), gene_table$gene_id), ]
stopifnot(identical(g$gene_id, rownames(counts)))

############################################
# 1. Counts
############################################

write10xCounts(file.path(out, "counts.h5"), counts,
               barcodes = colnames(counts),
               gene.id = g$gene_id, gene.symbol = g$gene_name,
               type = "HDF5", version = "3", genome = "GRCh38", overwrite = TRUE)

############################################
# 2. Cell info (obs) and gene info (var)
############################################

obs <- data.frame(cell = colnames(seu), seu@meta.data, check.names = FALSE)
write.csv(obs, file.path(out, "obs.csv"), row.names = FALSE)

var <- data.frame(gene_id = g$gene_id, gene_name = g$gene_name,
                  gene_biotype = g$gene_biotype,
                  is_mt = grepl("^MT-", g$gene_name),
                  is_rb = grepl("^RPS|^RPL", g$gene_name))
write.csv(var, file.path(out, "var.csv"), row.names = FALSE)

cat("Nuclei:", ncol(counts), "| Genes:", nrow(counts), "\n")
cat("MT genes:", sum(var$is_mt), "| RPS/RPL genes:", sum(var$is_rb), "\n")
