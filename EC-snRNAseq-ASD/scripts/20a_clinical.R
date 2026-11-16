############################################
# 20a_clinical.R - per-donor ASD profiles, clustering for patient subgroups (stop for review),
#                  signature scores vs ADI-R and co-occurring conditions
# Run in R from the folder containing 12_annotated.rds, gene_table_biotype.csv,
# clinical_rnaseq_merged.csv and 14_DEG/MAST/
# Packages: Seurat, Matrix, edgeR, pheatmap, cluster, ggplot2
#
# After this script: review 20_clinical/draft_heatmaps/ and 20_clinical/cluster_k_summary.csv,
# then edit 20_clinical/group_assignments.csv (columns `group` and `group_name`) before running 20b.
############################################

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(edgeR); library(pheatmap); library(cluster); library(ggplot2)
})

############################################
# Settings
############################################

out          <- "20_clinical"
clin_file    <- "clinical_rnaseq_merged.csv"
mast_dir     <- "14_DEG/MAST"
run_subtypes <- c("Ex-UL", "Ex-DL", "PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+",
                  "Astrocyte", "Oligodendrocyte", "OPC", "Microglia")      # Endothelial skipped (5 DEGs)
min_nuclei_donor <- 10
n_genes      <- 200          # genes per heatmap (variable genes and DEGs)
min_degs     <- 10           # DEG heatmap/score made only if >= 10 DEGs are in the pseudobulk data
fdr          <- 0.05
lfc_cut      <- 0.25
fc_cap       <- 2            # log2FC capped at +/-2 for clustering and display
k_range      <- 2:5
n_perm       <- 10000
seed         <- 42
adi_cols     <- c("ADI_A", "ADI_B", "ADI_C", "ADI_Combined")
conditions   <- c("ADHD", "Anxiety", "Nonverbal", "ID", "Epilepsy")
dx_cols      <- c(Control = "#4DBBD5", ASD = "#E64B35")

for (d in c(out, file.path(out, "draft_heatmaps"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

############################################
# Load
############################################

seu <- readRDS("12_annotated.rds")
counts <- LayerData(seu, assay = "RNA", layer = "counts")
md <- seu@meta.data
rm(seu); gc()
gt <- read.csv("gene_table_biotype.csv")
sym <- setNames(gt$gene_name, gt$gene_id)

clin <- read.csv(clin_file, check.names = FALSE, stringsAsFactors = FALSE)
for (cn in conditions) clin[[cn]] <- factor(clin[[cn]], levels = c("N", "Y"))
rownames(clin) <- clin$Sample_ID
donors <- unique(md[, c("Sample_ID", "Dx", "Age", "PMI")]); rownames(donors) <- donors$Sample_ID
asd_ids <- donors$Sample_ID[donors$Dx == "ASD"]
miss <- setdiff(asd_ids, clin$Sample_ID)
if (length(miss)) cat("ASD donors without clinical data:", paste(miss, collapse = ", "), "\n")
cat("ASD donors:", length(asd_ids), "| with ADI-R:", sum(clin$has_ADIR), "\n")

# Annotation table for heatmaps
ann_base <- clin[, c(adi_cols, conditions)]
ann_base$Age <- donors[rownames(ann_base), "Age"]; ann_base$PMI <- donors[rownames(ann_base), "PMI"]
ann_colors_base <- setNames(lapply(conditions, function(x) c(N = "grey92", Y = "grey25")), conditions)

############################################
# Per subtype
############################################

objs <- list(); scores <- list(); ksum <- list(); assign <- list()
heat_cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
brks <- seq(-fc_cap, fc_cap, length.out = 102)

for (st in run_subtypes) {
  cat("\n====", st, "====\n")
  fs <- gsub("[^A-Za-z0-9_-]", "", st)
  cells <- rownames(md)[md$subtype == st]
  n_d <- table(as.character(md[cells, "Sample_ID"]))
  keep_d <- names(n_d)[n_d >= min_nuclei_donor]
  sel <- md[cells, "Sample_ID"] %in% keep_d
  d_f <- factor(as.character(md[cells[sel], "Sample_ID"]), levels = keep_d)
  pb <- as.matrix(counts[, cells[sel]] %*% t(fac2sparse(d_f)))[, keep_d]
  dx <- donors[keep_d, "Dx"]
  keep_g <- filterByExpr(pb, group = dx)
  lc <- cpm(DGEList(pb[keep_g, ]), log = TRUE, prior.count = 1)

  ctrl <- keep_d[dx == "Control"]; asd <- intersect(keep_d[dx == "ASD"], clin$Sample_ID)
  ctrl_mean <- rowMeans(lc[, ctrl, drop = FALSE])
  fc_asd <- lc[, asd, drop = FALSE] - ctrl_mean
  fc_ctrl <- sapply(ctrl, function(cc) lc[, cc] - rowMeans(lc[, setdiff(ctrl, cc), drop = FALSE]))  # leave-one-out
  named <- rownames(lc)[!is.na(sym[rownames(lc)]) & !grepl("^ENSG", sym[rownames(lc)])]
  cat("Donors: Control", length(ctrl), "| ASD", length(asd), "| genes", nrow(lc), "\n")

  # ---- Gene sets ----
  sets <- list()
  v <- apply(fc_asd[named, , drop = FALSE], 1, sd)
  sets$variable <- names(sort(v, decreasing = TRUE))[seq_len(min(n_genes, length(v)))]
  degf <- file.path(mast_dir, paste0("DEG_MAST_", st, ".csv"))
  up <- dn <- character(0)
  if (file.exists(degf)) {
    dt <- read.csv(degf); dt$gene_id <- sub("\\..*$", "", dt$gene_id)
    dt <- dt[!is.na(dt$padj) & dt$padj < fdr & abs(dt$log2FC) >= lfc_cut & dt$gene_id %in% rownames(lc), ]
    dt <- dt[order(dt$padj, dt$pvalue), ]
    up <- dt$gene_id[dt$log2FC > 0]; dn <- dt$gene_id[dt$log2FC < 0]
    dsel <- head(dt$gene_id[dt$gene_id %in% named], n_genes)
    if (length(dsel) >= min_degs) sets$DEG <- dsel else cat("Fewer than", min_degs, "DEGs - DEG heatmap skipped\n")
    cat("MAST DEGs in pseudobulk data: up", length(up), "| down", length(dn), "\n")
  } else cat("No MAST table:", degf, "\n")

  # ---- Signature score: mean FC of up-DEGs minus mean FC of down-DEGs ----
  if (length(up) + length(dn) >= min_degs) {
    sc_fun <- function(fc) (if (length(up)) colMeans(fc[up, , drop = FALSE]) else 0) -
                           (if (length(dn)) colMeans(fc[dn, , drop = FALSE]) else 0)
    scores[[st]] <- rbind(
      data.frame(subtype = st, Sample_ID = asd, Dx = "ASD", score = sc_fun(fc_asd)),
      data.frame(subtype = st, Sample_ID = ctrl, Dx = "Control", score = sc_fun(fc_ctrl)))
  }

  # ---- Clustering of ASD donors per gene set ----
  for (gs in names(sets)) {
    M <- pmax(pmin(fc_asd[sets[[gs]], , drop = FALSE], fc_cap), -fc_cap)
    d <- dist(t(M)); hc <- hclust(d, method = "ward.D2")
    ks <- k_range[k_range < ncol(M)]
    sil <- sapply(ks, function(k) mean(silhouette(cutree(hc, k), d)[, 3]))
    best <- ks[which.max(sil)]
    memb <- sapply(ks, function(k) cutree(hc, k)); colnames(memb) <- paste0("k", ks)
    ksum[[paste(st, gs)]] <- data.frame(subtype = st, gene_set = gs, k = ks, silhouette = round(sil, 3),
      group_sizes = sapply(ks, function(k) paste0("n = ", paste(table(cutree(hc, k)), collapse = " | "))),   # text, so Excel won't read it as a date
      suggested = ks == best)
    assign[[paste(st, gs)]] <- data.frame(subtype = st, gene_set = gs, Sample_ID = colnames(M),
      memb, group = cutree(hc, best), group_name = "", row.names = NULL)
    cat(sprintf("  %-8s genes %3d | silhouette k=%s: %s | suggested k = %d\n", gs, nrow(M),
                paste(ks, collapse = ","), paste(round(sil, 2), collapse = ", "), best))

    ann <- ann_base[colnames(M), , drop = FALSE]
    ann$Draft_group <- factor(paste0("G", cutree(hc, best)))
    if (!is.null(scores[[st]])) ann$Signature <- scores[[st]]$score[match(colnames(M), scores[[st]]$Sample_ID)]
    ann <- ann[, rev(colnames(ann))]
    acol <- c(ann_colors_base, list(Draft_group = setNames(
      c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")[seq_len(best)], paste0("G", seq_len(best)))))
    title <- paste0(st, " - ", ifelse(gs == "DEG", paste0("top ", nrow(M), " MAST DEGs"),
                                       paste0("top ", nrow(M), " variable genes")),
                    " | per-donor log2FC vs Control mean | DRAFT (k = ", best, ", best silhouette)")
    pheatmap(M, color = heat_cols, breaks = brks, cluster_cols = hc, cutree_cols = best,
             clustering_method = "ward.D2", clustering_distance_rows = "euclidean",
             labels_row = sym[rownames(M)], fontsize_row = if (nrow(M) > 100) 3 else 5, fontsize_col = 6,
             annotation_col = ann, annotation_colors = acol, main = title, border_color = NA,
             filename = file.path(out, "draft_heatmaps", paste0("draft_", fs, "_", gs, ".pdf")),
             width = 9, height = max(6, 0.04 * nrow(M) + 4))
  }
  objs[[st]] <- list(fc_asd = fc_asd, fc_ctrl = fc_ctrl, sets = sets, up = up, down = dn)
}

############################################
# Save clustering outputs for review
############################################

K <- do.call(rbind, ksum); A <- do.call(rbind, assign)
write.csv(K, file.path(out, "cluster_k_summary.csv"), row.names = FALSE)
af <- file.path(out, "group_assignments.csv")
if (file.exists(af)) {
  af <- file.path(out, "group_assignments_NEW.csv")
  cat("\ngroup_assignments.csv already exists and was NOT overwritten; new draft saved as group_assignments_NEW.csv\n")
}
write.csv(A, af, row.names = FALSE)
saveRDS(objs, file.path(out, "20a_objects.rds"))

############################################
# Signature score vs ADI-R and conditions
############################################

S <- do.call(rbind, scores)
S <- merge(S, clin[, c("Sample_ID", adi_cols, conditions, "has_ADIR")], by = "Sample_ID", all.x = TRUE)
write.csv(S, file.path(out, "signature_scores.csv"), row.names = FALSE)

perm_cor <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  r <- cor(x, y, method = "spearman")
  null <- replicate(n_perm, cor(x, sample(y), method = "spearman"))
  c(n = length(x), rho = r, p_perm = (1 + sum(abs(null) >= abs(r))) / (n_perm + 1))
}
R <- do.call(rbind, lapply(names(scores), function(st) {
  s <- S[S$subtype == st & S$Dx == "ASD" & S$has_ADIR %in% TRUE, ]
  do.call(rbind, lapply(adi_cols, function(a) data.frame(subtype = st, ADI = a, t(perm_cor(s$score, s[[a]])))))
}))
R$fdr <- p.adjust(R$p_perm, "BH")
write.csv(R, file.path(out, "signature_vs_ADIR.csv"), row.names = FALSE)

Cn <- do.call(rbind, lapply(names(scores), function(st) {
  s <- S[S$subtype == st & S$Dx == "ASD", ]
  do.call(rbind, lapply(conditions, function(cnd) {
    y <- s$score[s[[cnd]] %in% "Y"]; n <- s$score[s[[cnd]] %in% "N"]
    if (length(y) < 3 || length(n) < 3) return(NULL)
    data.frame(subtype = st, condition = cnd, n_Y = length(y), n_N = length(n),
               median_Y = median(y), median_N = median(n),
               p = wilcox.test(y, n, exact = FALSE)$p.value)
  }))
}))
Cn$fdr <- p.adjust(Cn$p, "BH")
write.csv(Cn, file.path(out, "signature_vs_conditions.csv"), row.names = FALSE)

cat("\nSignature vs ADI-R (nominal p < 0.05):\n"); print(R[R$p_perm < 0.05, ], row.names = FALSE)
cat("\nSignature vs conditions (nominal p < 0.05):\n"); print(Cn[Cn$p < 0.05, ], row.names = FALSE)

############################################
# Figures for scores
############################################

S$subtype <- factor(S$subtype, levels = run_subtypes); S$Dx <- factor(S$Dx, levels = c("Control", "ASD"))
p1 <- ggplot(S, aes(Dx, score, color = Dx)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_jitter(width = 0.15, size = 1.2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.3) +
  facet_wrap(~ subtype, scales = "free_y", ncol = 6) + scale_color_manual(values = dx_cols) +
  labs(title = "ASD signature score per donor", y = "mean log2FC(up DEGs) - mean log2FC(down DEGs)", x = NULL,
       caption = "ASD donors vs Control mean; each Control vs the other Controls (leave-one-out)") +
  theme_classic(base_size = 8) + theme(legend.position = "none")

Sa <- S[S$Dx == "ASD" & S$has_ADIR %in% TRUE, ]
la <- do.call(rbind, lapply(adi_cols, function(a) data.frame(Sa[, c("subtype", "Sample_ID", "score")], ADI = a, value = Sa[[a]])))
lab <- R; lab$txt <- paste0("rho = ", round(lab$rho, 2), ", p = ", signif(lab$p_perm, 2))
p2 <- ggplot(la, aes(value, score)) + geom_point(size = 1) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, color = "grey40") +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.3, size = 2, inherit.aes = FALSE) +
  facet_grid(subtype ~ ADI, scales = "free") +
  labs(title = "Signature score vs ADI-R (Spearman, permutation p)", x = "ADI-R score", y = "Signature score") +
  theme_bw(base_size = 7)

Sc <- S[S$Dx == "ASD", ]
lc2 <- do.call(rbind, lapply(conditions, function(cnd) data.frame(Sc[, c("subtype", "score")], condition = cnd, status = Sc[[cnd]])))
lc2 <- lc2[!is.na(lc2$status), ]
p3 <- ggplot(lc2, aes(status, score)) + geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, size = 0.8) + facet_grid(subtype ~ condition, scales = "free_y") +
  labs(title = "Signature score by co-occurring condition (ASD donors; Wilcoxon in table)", x = NULL, y = "Signature score") +
  theme_bw(base_size = 7)

pdf(file.path(out, "signature_score_figures.pdf"), width = 10, height = 6); print(p1); dev.off()
pdf(file.path(out, "signature_vs_ADIR.pdf"), width = 7, height = 1.3 * length(scores) + 1); print(p2); dev.off()
pdf(file.path(out, "signature_vs_conditions.pdf"), width = 8, height = 1.3 * length(scores) + 1); print(p3); dev.off()

writeLines(c(
  paste0("Pseudobulk per donor x subtype (donors with >= ", min_nuclei_donor, " nuclei); filterByExpr(group = Dx); log2 CPM (prior.count = 1)"),
  "Per-donor log2FC: ASD donor - mean of Controls; Control donors - mean of the other Controls (leave-one-out)",
  paste0("Variable genes: top ", n_genes, " by SD of log2FC across ASD donors (genes with a symbol)"),
  paste0("DEG genes: top ", n_genes, " MAST DEGs (FDR < ", fdr, ", |log2FC| >= ", lfc_cut, ") by FDR, genes with a symbol; heatmap only if >= ", min_degs),
  paste0("Clustering: ASD donors, Euclidean distance on log2FC capped at +/-", fc_cap, ", Ward D2; suggested k = best mean silhouette among k = ", paste(range(k_range), collapse = "-")),
  "Signature score: mean log2FC of all up-DEGs minus mean log2FC of all down-DEGs (MAST DEGs present in the pseudobulk data)",
  paste0("ADI-R: Spearman, two-sided permutation p (", n_perm, "), ASD donors with ADI-R only; BH across subtypes x domains"),
  "Conditions: Wilcoxon (Y vs N, ASD donors), conditions with >= 3 donors per side; BH across subtypes x conditions",
  paste("Seed:", seed)), file.path(out, "clinical_20a_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_20a.txt"))

cat("\nDone. Review draft_heatmaps/ and cluster_k_summary.csv, then edit group_assignments.csv before 20b.\n")
