############################################
# 14b_mast.R — nucleus-level DEG (ASD vs Control) per subtype with MAST (fixed effects only)
# Run in R from the folder containing 12_annotated.rds and gene_table_biotype.csv
# Packages: Seurat, Matrix, MAST, data.table, ggplot2, ggrepel
# Runtime: the largest subtypes (Ex-UL ~22k nuclei) can take an hour or more each on CPU
############################################

library(Seurat)
library(Matrix)
library(MAST)
library(data.table)
library(ggplot2)
library(ggrepel)

############################################
# Settings
############################################

out <- "14_DEG/MAST"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

group_by     <- "subtype"
min_pct      <- 0.10         # gene detected in >= 10% of nuclei in Control OR ASD
fdr          <- 0.05
lfc_report   <- 0.25         # reporting only (second DEG count column); not a test threshold
model        <- ~ cngeneson + Age_z + PMI_z + Dx
coef_name    <- "DxASD"
run_subtypes <- NULL         # NULL = all 12; or e.g. c("PVALB+", "SST+") to run a subset

############################################
# 1. Load
############################################

seu <- readRDS("12_annotated.rds")
counts <- LayerData(seu, assay = "RNA", layer = "counts")
md <- seu@meta.data
md$Dx <- factor(md$Dx, levels = c("Control", "ASD"))
if (anyNA(md[, c("Age", "PMI")])) stop("Age or PMI missing for some nuclei")

gt <- read.csv("gene_table_biotype.csv")
gt <- gt[match(rownames(counts), gt$gene_id), c("gene_id", "gene_name", "gene_biotype")]
stopifnot(identical(gt$gene_id, rownames(counts)))
cat("Genes in object (all biotypes):", nrow(counts), "\n")

subtypes <- levels(md[[group_by]])
if (is.null(subtypes)) subtypes <- sort(unique(md[[group_by]]))
if (!is.null(run_subtypes)) subtypes <- intersect(subtypes, run_subtypes)

volcano <- function(tab, title) {
  tab <- tab[!is.na(tab$padj), ]
  tab$sig <- ifelse(tab$padj < fdr & tab$log2FC > 0, "Up in ASD",
             ifelse(tab$padj < fdr & tab$log2FC < 0, "Down in ASD", "NS"))
  lab <- head(tab[tab$sig != "NS", ][order(tab$padj[tab$sig != "NS"]), ], 15)
  ggplot(tab, aes(log2FC, -log10(pvalue), color = sig)) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_text_repel(data = lab, aes(label = gene_name), size = 2.5, color = "black",
                    max.overlaps = 30, segment.size = 0.2) +
    scale_color_manual(values = c("Up in ASD" = "#E64B35", "Down in ASD" = "#4DBBD5", "NS" = "grey75")) +
    labs(title = title, x = "log2 fold change (ASD vs Control)", y = "-log10 p", color = NULL) +
    theme_classic(base_size = 9)
}

############################################
# 2. Per subtype: filter genes, log2(CPM+1), MAST hurdle model
############################################

summary_rows <- list()

for (st in subtypes) {
  t0 <- Sys.time()
  cat("\n====", st, "====\n")
  cells <- rownames(md)[md[[group_by]] == st]
  cnt <- counts[, cells]
  cd <- md[cells, c("Sample_ID", "Dx", "Age", "PMI")]
  cat("Nuclei: Control", sum(cd$Dx == "Control"), "| ASD", sum(cd$Dx == "ASD"),
      "| donors:", length(unique(cd$Sample_ID)), "\n")

  # Gene filter: detected in >= min_pct of nuclei in either group
  det <- cnt > 0
  pct_c <- Matrix::rowMeans(det[, cd$Dx == "Control", drop = FALSE])
  pct_a <- Matrix::rowMeans(det[, cd$Dx == "ASD", drop = FALSE])
  keep <- pct_c >= min_pct | pct_a >= min_pct
  cat("Genes passing", 100 * min_pct, "% filter:", sum(keep), "\n")

  # log2(CPM + 1); library size from all genes
  lib <- Matrix::colSums(cnt)
  x <- cnt[keep, ] %*% Diagonal(x = 1e6 / lib)        # CPM, still sparse
  x@x <- log2(x@x + 1)                                # log2(CPM + 1); zeros stay zero
  logcpm <- as.matrix(x)
  dimnames(logcpm) <- list(rownames(cnt)[keep], cells)
  rm(x)

  cdat <- data.frame(wellKey = cells, Dx = cd$Dx,
                     Age_z = as.numeric(scale(cd$Age)), PMI_z = as.numeric(scale(cd$PMI)))
  fdat <- data.frame(primerid = rownames(logcpm))
  sca <- FromMatrix(logcpm, cData = cdat, fData = fdat, check_sanity = FALSE)
  colData(sca)$cngeneson <- as.numeric(scale(colSums(assay(sca) > 0)))   # cellular detection rate
  rm(logcpm, det); gc()

  zfit <- zlm(model, sca, method = "bayesglm", ebayes = TRUE)
  sm <- summary(zfit, doLRT = coef_name)$datatable

  p  <- sm[contrast == coef_name & component == "H", .(primerid, pvalue = `Pr(>Chisq)`)]
  fc <- sm[contrast == coef_name & component == "logFC", .(primerid, log2FC = coef, ci.hi, ci.lo)]
  tab <- as.data.frame(merge(p, fc, by = "primerid"))
  tab$padj <- p.adjust(tab$pvalue, method = "BH")
  tab <- data.frame(gt[match(tab$primerid, gt$gene_id), ], tab[, -1],
                    pct_Control = round(100 * pct_c[tab$primerid], 1),
                    pct_ASD = round(100 * pct_a[tab$primerid], 1), row.names = NULL)
  tab <- tab[order(tab$padj, tab$pvalue), ]
  write.csv(tab, file.path(out, paste0("DEG_MAST_", st, ".csv")), row.names = FALSE)

  s <- !is.na(tab$padj) & tab$padj < fdr
  summary_rows[[st]] <- data.frame(
    subtype = st, nuclei_Control = sum(cd$Dx == "Control"), nuclei_ASD = sum(cd$Dx == "ASD"),
    genes_tested = nrow(tab),
    up = sum(s & tab$log2FC > 0, na.rm = TRUE), down = sum(s & tab$log2FC < 0, na.rm = TRUE),
    up_lfc = sum(s & tab$log2FC >= lfc_report, na.rm = TRUE),
    down_lfc = sum(s & tab$log2FC <= -lfc_report, na.rm = TRUE),
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  cat("FDR <", fdr, ": up", summary_rows[[st]]$up, "| down", summary_rows[[st]]$down,
      "| time", summary_rows[[st]]$minutes, "min\n")

  pdf(file.path(out, paste0("volcano_MAST_", st, ".pdf")), width = 7, height = 5.5)
  print(volcano(tab, paste0(st, " - MAST, fixed effects (FDR < ", fdr, ")")))
  print(ggplot(tab[!is.na(tab$pvalue), ], aes(pvalue)) + geom_histogram(bins = 50, fill = "grey50") +
          labs(title = paste(st, "- p-value distribution")) + theme_classic(base_size = 9))
  dev.off()

  rm(sca, zfit, sm, cnt); gc()
  # Save summary after every subtype so a crash doesn't lose finished results
  write.csv(do.call(rbind, summary_rows), file.path(out, "DEG_summary_MAST.csv"), row.names = FALSE)
}

############################################
# 3. Summary
############################################

summary_tab <- do.call(rbind, summary_rows)
cat("\nSummary (up/down at FDR <", fdr, "; *_lfc also |log2FC| >=", lfc_report, "):\n")
print(summary_tab, row.names = FALSE)

writeLines(c(
  paste("Grouping:", group_by),
  paste("Gene filter: all biotypes; detected in >=", 100 * min_pct, "% of nuclei in Control or ASD (per subtype)"),
  "Expression: log2(CPM + 1), library size = total counts per nucleus (all genes)",
  paste("Model: MAST zlm", deparse(model), "; method = bayesglm, ebayes = TRUE; fixed effects only (no donor term)"),
  "cngeneson: scaled number of genes detected per nucleus (after gene filter); Age, PMI centred and scaled within subtype",
  paste("Test: hurdle likelihood-ratio test for", coef_name, "; BH FDR <", fdr),
  "log2FC: MAST combined (discrete x continuous) coefficient for DxASD",
  "Caveat: nuclei treated as independent; p-values are not corrected for donor-level correlation",
  paste("Second DEG count uses |log2FC| >=", lfc_report, "; reporting only")),
  file.path(out, "MAST_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_14b.txt"))
cat("\nSaved", out, "/\n")
