############################################
# 05_doublets.R — scDblFinder per sample, remove doublets
############################################

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)
library(ggplot2)

out <- "05_doublets"
dir.create(out, showWarnings = FALSE)

seu <- readRDS("04_filtered.rds")

dbr  <- 0.0706   # expected doublet rate, same for every sample
seed <- 42

############################################
# 1. Check: each round-1 well holds only one sample
############################################

well_tab <- table(as.character(seu$bc1_well), seu$Sample_ID)
cat("Round-1 wells per sample:", range(colSums(well_tab > 0)), "\n")
cat("Max samples in any round-1 well (should be 1):", max(rowSums(well_tab > 0)), "\n")
if (max(rowSums(well_tab > 0)) > 1)
  stop("A round-1 well contains more than one sample; per-sample doublet detection would miss doublets.")

############################################
# 2. scDblFinder, each sample separately
############################################

sce <- SingleCellExperiment(list(counts = LayerData(seu, layer = "counts")),
                            colData = DataFrame(Sample_ID = seu$Sample_ID))
set.seed(seed)
sce <- scDblFinder(sce, samples = "Sample_ID", dbr = dbr)

seu$scDblFinder.score <- sce$scDblFinder.score
seu$scDblFinder.class <- as.character(sce$scDblFinder.class)
rm(sce); invisible(gc())

############################################
# 3. Tables
############################################

calls <- data.frame(cell = colnames(seu), Sample_ID = seu$Sample_ID, Dx = seu$Dx,
                    score = seu$scDblFinder.score, class = seu$scDblFinder.class)
write.csv(calls, file.path(out, "doublet_calls.csv"), row.names = FALSE)

per_sample <- seu@meta.data %>% group_by(Sample_ID, Dx) %>%
  summarise(nuclei = n(), doublets = sum(scDblFinder.class == "doublet"),
            pct_doublet = round(100 * doublets / nuclei, 1), .groups = "drop") %>%
  arrange(Dx, Sample_ID)
write.csv(per_sample, file.path(out, "doublets_per_sample.csv"), row.names = FALSE)

cat("\nOverall doublet rate:", round(100 * mean(seu$scDblFinder.class == "doublet"), 2), "%\n")
print(table(seu$Dx, seu$scDblFinder.class))
cat("Per-sample doublet %: range", range(per_sample$pct_doublet), "\n")

############################################
# 4. Figures
############################################

dx_cols <- c(Control = "#4DBBD5", ASD = "#E64B35")
th <- theme_bw(base_size = 10)

pdf(file.path(out, "doublet_figures.pdf"), width = 14, height = 7)

print(ggplot(per_sample, aes(reorder(Sample_ID, pct_doublet), pct_doublet, fill = Dx)) +
        geom_col() +
        geom_hline(yintercept = 100 * dbr, linetype = "dashed") +
        scale_fill_manual(values = dx_cols) +
        labs(title = "Doublets per sample (dashed = expected rate)", x = NULL, y = "% doublets") +
        th + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

print(seu@meta.data %>%
        select(scDblFinder.class, nFeature_RNA, nCount_RNA) %>%
        tidyr::pivot_longer(-scDblFinder.class, names_to = "metric", values_to = "value") %>%
        ggplot(aes(scDblFinder.class, value, fill = scDblFinder.class)) +
        geom_violin(scale = "width") + geom_boxplot(width = 0.1, outlier.shape = NA) +
        scale_y_log10() + facet_wrap(~ metric, scales = "free_y") +
        labs(title = "Doublets vs singlets", x = NULL, y = NULL) + th +
        theme(legend.position = "none"))

dev.off()

############################################
# 5. Remove doublets and save
############################################

seu <- subset(seu, subset = scDblFinder.class == "singlet")
cat("\nSinglets kept:", ncol(seu), "\n")
print(table(seu$Dx))

writeLines(c(
  paste("Date:", Sys.Date()),
  paste("scDblFinder version:", packageVersion("scDblFinder")),
  "Run separately for each sample (samples = Sample_ID)",
  paste("Expected doublet rate (dbr):", dbr, "for every sample"),
  paste("Random seed:", seed),
  paste("Doublets removed:", nrow(calls) - ncol(seu)),
  paste("Singlets kept:", ncol(seu))
), file.path(out, "doublet_settings.txt"))

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
saveRDS(seu, "05_singlets.rds")
