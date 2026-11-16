############################################
# 04_apply_filters.R — apply QC filters, gene detection filter, save filtered object
############################################

library(Seurat)
library(dplyr)

out <- "04_filtered"
dir.create(out, showWarnings = FALSE)

seu <- readRDS("02_qc_metrics.rds")
gene_table <- read.csv("gene_table_biotype.csv")

############################################
# Filter settings (same as 03_filter_preview.R)
############################################

exclude_samples <- c("77_4KGA_ERC", "57_RMGG_ERC", "95_K47B_ERC")
nFeature_min <- 500
nFeature_max <- 15000
nCount_min   <- 1000
nCount_max   <- 200000
mt_max       <- 10      # percent.mt on all 37 MT genes
rb_max       <- 1       # percent.rb on all ^RPS|^RPL genes (1,775)
min_nuclei_per_gene <- 10   # gene kept if detected in >= 10 retained nuclei

############################################
# 1. Filter nuclei
############################################

md <- seu@meta.data
md$keep <- !(md$Sample_ID %in% exclude_samples) &
  md$nFeature_RNA >= nFeature_min & md$nFeature_RNA <= nFeature_max &
  md$nCount_RNA   >= nCount_min   & md$nCount_RNA   <= nCount_max &
  md$percent.mt   <= mt_max &
  md$percent.rb   <  rb_max

before <- md %>% dplyr::count(Sample_ID, Dx, name = "nuclei_before")

############################################
# 2. Gene filter: detected in >= 10 retained nuclei (all biotypes eligible)
############################################

counts <- LayerData(seu, layer = "counts")[, md$keep]
n_nuclei_detected <- Matrix::rowSums(counts > 0)
keep_genes <- rownames(counts)[n_nuclei_detected >= min_nuclei_per_gene]
cat("Genes before:", nrow(counts), "| removed (<", min_nuclei_per_gene, "nuclei):",
    nrow(counts) - length(keep_genes), "| kept:", length(keep_genes), "\n")

kept_biotypes <- table(gene_table$gene_biotype[match(keep_genes, gene_table$gene_id)])
print(sort(kept_biotypes, decreasing = TRUE)[1:10])

############################################
# 3. Rebuild object
############################################

# Keep QC-time values (calculated with all 37 MT genes) under new names
md_f <- md[md$keep, setdiff(names(md), "keep")]
names(md_f)[names(md_f) == "nCount_RNA"]   <- "nCount_RNA_qc"
names(md_f)[names(md_f) == "nFeature_RNA"] <- "nFeature_RNA_qc"

counts <- counts[keep_genes, ]
rm(seu); invisible(gc())

seu <- CreateSeuratObject(counts = counts, meta.data = md_f,
                          min.cells = 0, min.features = 0)
rm(counts); invisible(gc())
seu

cat("Nuclei retained:", ncol(seu), "\n")
print(table(seu$Dx))
cat("Samples retained:", n_distinct(seu$Sample_ID), "\n")

############################################
# 4. Save outputs
############################################

per_sample <- before %>%
  left_join(seu@meta.data %>% dplyr::count(Sample_ID, name = "nuclei_after"), by = "Sample_ID") %>%
  mutate(nuclei_after = coalesce(nuclei_after, 0L),
         pct_retained = round(100 * nuclei_after / nuclei_before, 1),
         excluded_sample = Sample_ID %in% exclude_samples) %>%
  arrange(Dx, Sample_ID)
write.csv(per_sample, file.path(out, "nuclei_per_sample_before_after.csv"), row.names = FALSE)

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
writeLines(c(
  paste("Date:", Sys.Date()),
  "Input: 02_qc_metrics.rds (all 78,899 genes; QC metrics on all genes)",
  paste("Excluded samples:", paste(exclude_samples, collapse = ", ")),
  paste0("nFeature_RNA: ", fmt(nFeature_min), " to ", fmt(nFeature_max)),
  paste0("nCount_RNA: ", fmt(nCount_min), " to ", fmt(nCount_max)),
  paste0("percent.mt <= ", mt_max, " (all 37 MT genes: 13 protein-coding, 2 rRNA, 22 tRNA)"),
  paste0("percent.rb < ", rb_max, " (all genes matching ^RPS|^RPL, incl. pseudogenes)"),
  paste0("Gene filter after nucleus filtering: detected in >= ", min_nuclei_per_gene, " retained nuclei"),
  "nCount_RNA_qc / nFeature_RNA_qc / percent.mt / percent.rb = values used for filtering",
  "nCount_RNA / nFeature_RNA = recalculated on the final gene set",
  paste("Samples retained:", n_distinct(seu$Sample_ID)),
  paste("Nuclei retained:", fmt(ncol(seu))),
  paste("Genes retained:", fmt(nrow(seu)))
), file.path(out, "filters_used.txt"))

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
saveRDS(seu, "04_filtered.rds")
