############################################
# 14a_pseudobulk_deseq2.R — pseudobulk DEG (ASD vs Control) per subtype with DESeq2
# Run in R from the folder containing 12_annotated.rds and gene_table_biotype.csv
# Packages: Seurat, Matrix, DESeq2, edgeR, apeglm, ggplot2, ggrepel
############################################

library(Seurat)
library(Matrix)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)

############################################
# Settings
############################################

out <- "14_DEG/pseudobulk"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

group_by          <- "subtype"   # 12 subtypes
min_nuclei_donor  <- 10          # donor included for a subtype if it has >= 10 nuclei of that subtype
min_donors_per_dx <- 3           # subtype tested only if >= 3 donors per Dx group pass
fdr               <- 0.05
lfc_report        <- 0.25        # used only for the second DEG count column; not a test threshold
design_formula    <- ~ Age_z + PMI_z + Dx
coef_name         <- "Dx_ASD_vs_Control"

############################################
# 1. Load
############################################

seu <- readRDS("12_annotated.rds")
counts <- LayerData(seu, assay = "RNA", layer = "counts")
md <- seu@meta.data
md$Dx <- factor(md$Dx, levels = c("Control", "ASD"))

gt <- read.csv("gene_table_biotype.csv")
gt <- gt[match(rownames(counts), gt$gene_id), c("gene_id", "gene_name", "gene_biotype")]
stopifnot(identical(gt$gene_id, rownames(counts)))
cat("Genes in object (all biotypes):", nrow(counts), "\n")

donor_info <- unique(md[, c("Sample_ID", "Dx", "Age", "PMI")])
stopifnot(!any(duplicated(donor_info$Sample_ID)))
if (anyNA(donor_info[, c("Age", "PMI")])) stop("Age or PMI missing for some donors")
rownames(donor_info) <- donor_info$Sample_ID

subtypes <- levels(md[[group_by]])
if (is.null(subtypes)) subtypes <- sort(unique(md[[group_by]]))

############################################
# 2. Helpers
############################################

volcano <- function(tab, title, x = "log2FoldChange_shrunk") {
  tab <- tab[!is.na(tab$padj), ]
  tab$sig <- ifelse(tab$padj < fdr & tab[[x]] > 0, "Up in ASD",
             ifelse(tab$padj < fdr & tab[[x]] < 0, "Down in ASD", "NS"))
  lab <- head(tab[tab$sig != "NS", ][order(tab$padj[tab$sig != "NS"]), ], 15)
  ggplot(tab, aes(.data[[x]], -log10(pvalue), color = sig)) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_text_repel(data = lab, aes(label = gene_name), size = 2.5, color = "black",
                    max.overlaps = 30, segment.size = 0.2) +
    scale_color_manual(values = c("Up in ASD" = "#E64B35", "Down in ASD" = "#4DBBD5", "NS" = "grey75")) +
    labs(title = title, x = "log2 fold change (ASD vs Control)", y = "-log10 p", color = NULL) +
    theme_classic(base_size = 9)
}

############################################
# 3. Per subtype: aggregate, filter, DESeq2
############################################

summary_rows <- list()
pdf(file.path(out, "pseudobulk_figures.pdf"), width = 7, height = 5.5)

for (st in subtypes) {
  cat("\n====", st, "====\n")
  cells <- rownames(md)[md[[group_by]] == st]
  donor <- factor(md[cells, "Sample_ID"])
  n_per_donor <- table(donor)
  keep_donors <- names(n_per_donor)[n_per_donor >= min_nuclei_donor]

  cd <- donor_info[keep_donors, ]
  n_dx <- table(factor(cd$Dx, levels = c("Control", "ASD")))
  cat("Donors with >=", min_nuclei_donor, "nuclei: Control", n_dx["Control"], "| ASD", n_dx["ASD"], "\n")
  if (any(n_dx < min_donors_per_dx)) {
    cat("Skipped: fewer than", min_donors_per_dx, "donors in a Dx group\n")
    summary_rows[[st]] <- data.frame(subtype = st, tested = FALSE, donors_Control = n_dx[["Control"]],
                                     donors_ASD = n_dx[["ASD"]], nuclei = length(cells), genes_tested = NA,
                                     up = NA, down = NA, up_lfc = NA, down_lfc = NA)
    next
  }

  # Sum counts per donor (donors that pass only)
  sel <- donor %in% keep_donors
  d_sel <- droplevels(donor[sel])
  pb <- counts[, cells[sel]] %*% t(fac2sparse(d_sel))
  pb <- as.matrix(pb)[, keep_donors]

  cd$Age_z <- as.numeric(scale(cd$Age))
  cd$PMI_z <- as.numeric(scale(cd$PMI))
  cd$n_nuclei <- as.integer(n_per_donor[keep_donors])

  keep_genes <- filterByExpr(pb, group = cd$Dx)
  cat("Genes passing filterByExpr:", sum(keep_genes), "of", nrow(pb), "\n")

  dds <- DESeqDataSetFromMatrix(round(pb[keep_genes, ]), colData = cd, design = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  stopifnot(coef_name %in% resultsNames(dds))

  res <- results(dds, name = coef_name, alpha = fdr)
  shr <- tryCatch(lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE),
                  error = function(e) { cat("apeglm failed (", conditionMessage(e), "); using 'normal'\n");
                                        lfcShrink(dds, coef = coef_name, type = "normal", quiet = TRUE) })

  tab <- data.frame(gt[match(rownames(res), gt$gene_id), ],
                    baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                    log2FoldChange_shrunk = shr$log2FoldChange, stat = res$stat,
                    pvalue = res$pvalue, padj = res$padj, row.names = NULL)
  tab$donors_Control <- n_dx[["Control"]]; tab$donors_ASD <- n_dx[["ASD"]]
  tab <- tab[order(tab$padj, tab$pvalue), ]
  write.csv(tab, file.path(out, paste0("DEG_pseudobulk_", st, ".csv")), row.names = FALSE)
  write.csv(cbind(gene_id = rownames(pb), gene_name = gt$gene_name[match(rownames(pb), gt$gene_id)], pb),
            file.path(out, paste0("pseudobulk_counts_", st, ".csv")), row.names = FALSE)

  s <- !is.na(tab$padj) & tab$padj < fdr
  summary_rows[[st]] <- data.frame(
    subtype = st, tested = TRUE, donors_Control = n_dx[["Control"]], donors_ASD = n_dx[["ASD"]],
    nuclei = sum(sel), genes_tested = nrow(tab),
    up = sum(s & tab$log2FoldChange > 0), down = sum(s & tab$log2FoldChange < 0),
    up_lfc = sum(s & tab$log2FoldChange >= lfc_report), down_lfc = sum(s & tab$log2FoldChange <= -lfc_report))
  cat("FDR <", fdr, ": up", summary_rows[[st]]$up, "| down", summary_rows[[st]]$down, "\n")
  cat("Significant genes by biotype:\n"); print(table(tab$gene_biotype[s]))

  # QC: PCA of donors (variance-stabilized) and p-value histogram
  vs <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
  pc <- plotPCA(vs, intgroup = "Dx", returnData = TRUE, ntop = 500)
  pv <- round(100 * attr(pc, "percentVar"))
  print(ggplot(pc, aes(PC1, PC2, color = Dx)) + geom_point(size = 2) +
          geom_text_repel(aes(label = name), size = 2, show.legend = FALSE) +
          scale_color_manual(values = c(Control = "#4DBBD5", ASD = "#E64B35")) +
          labs(title = paste(st, "- pseudobulk PCA (top 500 genes)"),
               x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) +
          theme_classic(base_size = 9))
  print(ggplot(tab[!is.na(tab$pvalue), ], aes(pvalue)) + geom_histogram(bins = 50, fill = "grey50") +
          labs(title = paste(st, "- p-value distribution")) + theme_classic(base_size = 9))
  print(volcano(tab, paste0(st, " - pseudobulk DESeq2 (FDR < ", fdr, ")")))
}
dev.off()

############################################
# 4. Summary
############################################

summary_tab <- do.call(rbind, summary_rows)
write.csv(summary_tab, file.path(out, "DEG_summary_pseudobulk.csv"), row.names = FALSE)
cat("\nSummary (up/down at FDR <", fdr, "; *_lfc also |log2FC| >=", lfc_report, "):\n")
print(summary_tab, row.names = FALSE)

writeLines(c(
  paste("Grouping:", group_by),
  paste("Pseudobulk: raw counts summed per donor x", group_by),
  paste("Donor inclusion: >=", min_nuclei_donor, "nuclei; subtype tested if >=", min_donors_per_dx, "donors per Dx"),
  "Genes: all biotypes; per-subtype filter edgeR::filterByExpr(group = Dx), default settings",
  paste("Design:", deparse(design_formula), "(Age, PMI centred and scaled within each subtype)"),
  paste("Test: DESeq2 Wald,", coef_name, "; BH FDR <", fdr),
  "LFC shrinkage: apeglm (fallback 'normal'), used for plots only; tables report both",
  paste("Second DEG count uses |log2FC| >=", lfc_report, "(unshrunk); reporting only")),
  file.path(out, "pseudobulk_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_14a.txt"))
cat("\nSaved", out, "/\n")
