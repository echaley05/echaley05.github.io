############################################
# 18b_age_by_dx_pseudobulk.R - Age x Diagnosis interaction per subtype, donor level (pseudobulk DESeq2)
# 1. Main: all donors with >= 10 nuclei; ~ PMI + Age * Dx; LRT for Age:Dx
# 2. Sensitivity: donors restricted to the age range shared by both groups (drops the oldest Controls)
# 3. Comparison with the MAST interaction results from 18_age_by_dx.R
# 4. Donor-level plots of the top genes
# Run in R from the folder containing 12_annotated.rds, gene_table_biotype.csv and 18_age_by_dx/
# Packages: Seurat, Matrix, DESeq2, edgeR, ggplot2
############################################

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(DESeq2); library(edgeR); library(ggplot2)
})

############################################
# Settings
############################################

out          <- "18_age_by_dx/pseudobulk"
mast_dir     <- "18_age_by_dx"
group_by     <- "subtype"
run_subtypes <- c("Ex-UL", "Ex-DL", "PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+",
                  "Astrocyte", "Oligodendrocyte", "OPC", "Microglia")
min_nuclei_donor  <- 10
min_donors_per_dx <- 3
age_unit     <- 10          # slopes per decade
fdr          <- 0.05
n_mast_top   <- 200         # top MAST interaction genes compared against pseudobulk
n_plot       <- 8
dx_cols      <- c(Control = "#4DBBD5", ASD = "#E64B35")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

############################################
# Load
############################################

seu <- readRDS("12_annotated.rds")
counts <- LayerData(seu, assay = "RNA", layer = "counts")
md <- seu@meta.data
rm(seu); gc()
md$Dx <- factor(md$Dx, levels = c("Control", "ASD"))
gt <- read.csv("gene_table_biotype.csv")
gt <- gt[match(rownames(counts), gt$gene_id), c("gene_id", "gene_name", "gene_biotype")]
donors <- unique(md[, c("Sample_ID", "Dx", "Age", "PMI")])
rownames(donors) <- donors$Sample_ID

############################################
# Core function: one pseudobulk interaction test
############################################

run_pb <- function(pb, cd, tag) {
  cd$Age_c <- (cd$Age - mean(cd$Age)) / age_unit
  cd$PMI_z <- as.numeric(scale(cd$PMI))
  keep <- filterByExpr(pb, group = cd$Dx)
  dds <- DESeqDataSetFromMatrix(round(pb[keep, ]), colData = cd, design = ~ PMI_z + Age_c * Dx)
  dds <- DESeq(dds, test = "LRT", reduced = ~ PMI_z + Age_c + Dx, quiet = TRUE)
  res <- results(dds)
  cf <- coef(dds)                                       # log2 scale
  ic <- grep("Age_c.*Dx|Dx.*Age_c", colnames(cf), value = TRUE)
  dc <- grep("^Dx", colnames(cf), value = TRUE)
  tab <- data.frame(
    gt[match(rownames(res), gt$gene_id), ],
    baseMean = res$baseMean,
    slope_Control_per_decade = cf[, "Age_c"],
    slope_ASD_per_decade = cf[, "Age_c"] + cf[, ic],
    interaction_slope_diff = cf[, ic],
    Dx_log2FC_at_mean_age = cf[, dc],
    p_interaction = res$pvalue, padj_interaction = res$padj, row.names = NULL)
  tab <- tab[order(tab$padj_interaction, tab$p_interaction), ]
  tab$rank <- seq_len(nrow(tab))
  list(tab = tab, dds = dds, n_genes = sum(keep), tag = tag)
}

############################################
# Per subtype
############################################

summary_rows <- list(); cmp_rows <- list()
pdf_p <- file.path(out, "pvalue_histograms.pdf"); pdf(pdf_p, width = 7, height = 3)

for (st in run_subtypes) {
  cat("\n====", st, "====\n")
  fs <- gsub("[^A-Za-z0-9_-]", "", st)
  cells <- rownames(md)[md[[group_by]] == st]
  n_d <- table(as.character(md[cells, "Sample_ID"]))
  keep_d <- names(n_d)[n_d >= min_nuclei_donor]
  sel <- md[cells, "Sample_ID"] %in% keep_d
  d_f <- factor(as.character(md[cells[sel], "Sample_ID"]), levels = keep_d)
  pb <- as.matrix(counts[, cells[sel]] %*% t(fac2sparse(d_f)))[, keep_d]
  cd <- donors[keep_d, ]

  # Main and sensitivity donor sets
  rng <- c(max(tapply(cd$Age, cd$Dx, min)), min(tapply(cd$Age, cd$Dx, max)))
  in_rng <- cd$Age >= rng[1] & cd$Age <= rng[2]
  sets <- list(main = rep(TRUE, nrow(cd)), shared_age = in_rng)

  res_list <- list()
  for (nm in names(sets)) {
    s <- sets[[nm]]
    ndx <- table(factor(cd$Dx[s], levels = c("Control", "ASD")))
    if (any(ndx < min_donors_per_dx)) { cat(nm, ": too few donors, skipped\n"); next }
    r <- tryCatch(run_pb(pb[, s, drop = FALSE], cd[s, ], nm),
                  error = function(e) { cat(nm, "failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(r)) next
    r$donors <- ndx
    res_list[[nm]] <- r
    write.csv(r$tab, file.path(out, paste0("AgeDx_pseudobulk_", fs, "_", nm, ".csv")), row.names = FALSE)
    hist(r$tab$p_interaction, breaks = 50, col = "grey60", border = "white",
         main = paste(st, "-", nm, "- Age:Dx p-values"), xlab = "p")
    cat(sprintf("%-11s donors C %d / A %d | genes %d | FDR<%.2f: %d | p<0.01: %d\n", nm,
                ndx[["Control"]], ndx[["ASD"]], r$n_genes, fdr,
                sum(r$tab$padj_interaction < fdr, na.rm = TRUE), sum(r$tab$p_interaction < 0.01, na.rm = TRUE)))
  }
  if (!length(res_list)) next
  m <- res_list$main$tab; s2 <- res_list$shared_age$tab
  sig_main <- m$gene_id[!is.na(m$padj_interaction) & m$padj_interaction < fdr]
  sig_both <- if (!is.null(s2)) intersect(sig_main, s2$gene_id[!is.na(s2$padj_interaction) & s2$padj_interaction < fdr]) else character(0)

  # ---- Comparison with MAST ----
  mf <- file.path(mast_dir, paste0("AgeDx_MAST_", fs, ".csv"))
  cmp <- NULL
  if (file.exists(mf)) {
    mt <- read.csv(mf)
    mt <- mt[order(mt$padj_interaction, mt$p_interaction), ]
    common <- intersect(mt$gene_id, m$gene_id)
    top <- head(mt$gene_id[mt$gene_id %in% common], n_mast_top)
    pm <- m[match(top, m$gene_id), ]; mm <- mt[match(top, mt$gene_id), ]
    rho <- suppressWarnings(cor(mt$interaction_slope_diff[match(common, mt$gene_id)],
                                m$interaction_slope_diff[match(common, m$gene_id)], method = "spearman"))
    cmp <- data.frame(
      subtype = st, genes_common = length(common), spearman_slope_diff_all_genes = round(rho, 3),
      mast_top = length(top),
      top_with_pseudobulk_p_lt_0.05 = sum(pm$p_interaction < 0.05, na.rm = TRUE),
      top_same_direction = sum(sign(pm$interaction_slope_diff) == sign(mm$interaction_slope_diff), na.rm = TRUE),
      top_median_pseudobulk_percentile = round(100 * median(pm$rank / nrow(m)), 1))
    cmp_rows[[st]] <- cmp
    both <- data.frame(gene_id = top, gene_name = pm$gene_name,
                       MAST_padj = mm$padj_interaction, MAST_slope_diff = mm$interaction_slope_diff,
                       PB_p = pm$p_interaction, PB_padj = pm$padj_interaction, PB_slope_diff = pm$interaction_slope_diff,
                       PB_slope_Control = pm$slope_Control_per_decade, PB_slope_ASD = pm$slope_ASD_per_decade)
    write.csv(both, file.path(out, paste0("MAST_top_vs_pseudobulk_", fs, ".csv")), row.names = FALSE)
  }

  summary_rows[[st]] <- data.frame(
    subtype = st,
    donors_Control = res_list$main$donors[["Control"]], donors_ASD = res_list$main$donors[["ASD"]],
    genes_tested = res_list$main$n_genes,
    FDR_main = length(sig_main), p01_main = sum(m$p_interaction < 0.01, na.rm = TRUE),
    shared_age_range = paste(rng, collapse = "-"),
    donors_shared_range = if (!is.null(s2)) sum(res_list$shared_age$donors) else NA,
    FDR_shared_age = if (!is.null(s2)) sum(s2$padj_interaction < fdr, na.rm = TRUE) else NA,
    FDR_both = length(sig_both),
    top_genes_main = paste(head(m$gene_name[!grepl("^ENSG", m$gene_name)], 10), collapse = ", "))

  # ---- Donor-level plots: FDR genes in both sets; else top pseudobulk genes that are also MAST top ----
  pick <- sig_both
  if (!length(pick) && !is.null(cmp)) {
    cand <- both$gene_id[order(both$PB_p)]
    pick <- cand
  }
  if (!length(pick)) pick <- m$gene_id
  pick <- head(pick[!grepl("^ENSG", gt$gene_name[match(pick, gt$gene_id)])], n_plot)
  nc <- log2(counts(res_list$main$dds, normalized = TRUE)[pick, , drop = FALSE] + 1)
  long <- do.call(rbind, lapply(pick, function(g) data.frame(
    gene = gt$gene_name[gt$gene_id == g], Sample_ID = colnames(nc), expr = nc[g, ],
    Dx = cd[colnames(nc), "Dx"], Age = cd[colnames(nc), "Age"],
    lab = paste0(gt$gene_name[gt$gene_id == g], "\nPB FDR = ", signif(m$padj_interaction[m$gene_id == g], 2)))))
  long$lab <- factor(long$lab, levels = unique(long$lab))
  p <- ggplot(long, aes(Age, expr, color = Dx)) + geom_point(size = 1.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    geom_vline(xintercept = rng, linetype = "dotted", color = "grey50") +
    facet_wrap(~ lab, scales = "free_y", ncol = 4) +
    scale_color_manual(values = dx_cols) +
    labs(title = paste0(st, " - Age x Dx, donor level (pseudobulk)"),
         subtitle = if (length(sig_both)) "Genes with FDR < 0.05 in both the main and shared-age-range analyses"
                    else "No genes at FDR < 0.05 in both analyses; showing top pseudobulk genes among the MAST top hits",
         caption = "Points: donor pseudobulk, log2(DESeq2-normalised counts + 1); lines: per-group linear fit (display); dotted: shared age range",
         x = "Age (years)", y = "log2(normalised counts + 1)", color = NULL) +
    theme_classic(base_size = 8) + theme(legend.position = "top", strip.text = element_text(face = "italic"),
                                         plot.caption = element_text(size = 5.5, color = "grey40"))
  ggsave(file.path(out, paste0("AgeDx_pseudobulk_top_genes_", fs, ".pdf")), p, width = 8,
         height = 2 * ceiling(length(pick) / 4) + 1.2)
  write.csv(do.call(rbind, summary_rows), file.path(out, "AgeDx_pseudobulk_summary.csv"), row.names = FALSE)
}
dev.off()

############################################
# Summary
############################################

S <- do.call(rbind, summary_rows)
cat("\nPseudobulk Age x Dx summary:\n"); print(S[, setdiff(colnames(S), "top_genes_main")], row.names = FALSE)
if (length(cmp_rows)) {
  C <- do.call(rbind, cmp_rows)
  write.csv(C, file.path(out, "MAST_vs_pseudobulk_comparison.csv"), row.names = FALSE)
  cat("\nMAST top", n_mast_top, "interaction genes in pseudobulk:\n"); print(C, row.names = FALSE)
}

writeLines(c(
  "Pseudobulk: raw counts summed per donor x subtype; donors with >= 10 nuclei; subtype tested if >= 3 donors per Dx",
  "Genes: all biotypes; edgeR::filterByExpr(group = Dx) per subtype and donor set",
  paste0("Model: DESeq2 ~ PMI_z + Age_c * Dx; LRT vs ~ PMI_z + Age_c + Dx; BH FDR < ", fdr),
  paste0("Age_c = (Age - mean donor age) / ", age_unit, "; slopes are log2 fold change per decade (Control = Age_c coefficient; ASD = Age_c + interaction)"),
  "Dx_log2FC_at_mean_age: DxASD coefficient (ASD - Control at the mean donor age); no separate test",
  "Sensitivity: donors restricted to the age range covered by both groups in that subtype",
  paste0("MAST comparison: top ", n_mast_top, " MAST interaction genes (18_age_by_dx.R) located in the pseudobulk ranking; Spearman of slope differences over shared genes")),
  file.path(out, "AgeDx_pseudobulk_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_18b.txt"))
cat("\nSaved", out, "/\n")
