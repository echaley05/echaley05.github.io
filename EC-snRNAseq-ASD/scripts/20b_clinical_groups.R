############################################
# 20b_clinical_groups.R - final patient-subgroup heatmaps and clinical enrichment tests
# Uses 20_clinical/20a_objects.rds and 20_clinical/group_assignments.csv from 20a_clinical.R
# Groups: the `group` column (pre-filled with the best-silhouette k); optional `group_name`
# Run in R from the same folder as 20a
# Packages: pheatmap, ggplot2
############################################

suppressPackageStartupMessages({ library(pheatmap); library(ggplot2) })

############################################
# Settings
############################################

out        <- "20_clinical"
fd         <- file.path(out, "final")
clin_file  <- "clinical_rnaseq_merged.csv"
fc_cap     <- 2
n_sim      <- 10000          # Fisher simulated p-values
seed       <- 42
adi_cols   <- c("ADI_A", "ADI_B", "ADI_C", "ADI_Combined")
adi_labels <- c(ADI_A = "A: social", ADI_B = "B: communication", ADI_C = "C: restricted/repetitive",
                ADI_Combined = "Combined")
conditions <- c("ADHD", "Anxiety", "Nonverbal", "ID", "Epilepsy")
min_group_flag <- 3          # groups smaller than this are flagged
grp_cols   <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")
dpi        <- 600
dir.create(fd, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

############################################
# Load
############################################

objs <- readRDS(file.path(out, "20a_objects.rds"))
A <- read.csv(file.path(out, "group_assignments.csv"), stringsAsFactors = FALSE)
A$group_name[is.na(A$group_name)] <- ""
A$label <- ifelse(nzchar(A$group_name), A$group_name, paste0("G", A$group))
clin <- read.csv(clin_file, stringsAsFactors = FALSE); rownames(clin) <- clin$Sample_ID
for (cn in conditions) clin[[cn]] <- factor(clin[[cn]], levels = c("N", "Y"))
gt <- read.csv("gene_table_biotype.csv"); sym <- setNames(gt$gene_name, gt$gene_id)

# Dunn's test (rank-based post-hoc after Kruskal-Wallis, tie-corrected), BH within the comparison
dunn <- function(x, g) {
  ok <- !is.na(x); x <- x[ok]; g <- factor(g[ok]); N <- length(x); r <- rank(x)
  ties <- table(x); tc <- sum(ties^3 - ties) / (12 * (N - 1))
  rb <- tapply(r, g, mean); n <- table(g); prs <- combn(levels(g), 2)
  z <- apply(prs, 2, function(p) (rb[p[1]] - rb[p[2]]) /
               sqrt((N * (N + 1) / 12 - tc) * (1 / n[p[1]] + 1 / n[p[2]])))
  p <- 2 * pnorm(-abs(z))
  data.frame(comparison = paste(prs[1, ], "vs", prs[2, ]), z = round(z, 3), p = p, p_adj = p.adjust(p, "BH"))
}

heat_cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
brks <- seq(-fc_cap, fc_cap, length.out = 102)
adi_rows <- list(); cond_rows <- list(); dunn_rows <- list(); prof_rows <- list(); size_rows <- list()
plot_adi <- list(); plot_cond <- list()

############################################
# Per subtype x gene set
############################################

for (key in unique(paste(A$subtype, A$gene_set, sep = "||"))) {
  st <- sub("\\|\\|.*", "", key); gs <- sub(".*\\|\\|", "", key)
  if (is.null(objs[[st]]) || is.null(objs[[st]]$sets[[gs]])) next
  a <- A[A$subtype == st & A$gene_set == gs, ]
  ob <- objs[[st]]; fs <- gsub("[^A-Za-z0-9_-]", "", st)
  M <- pmax(pmin(ob$fc_asd[ob$sets[[gs]], a$Sample_ID, drop = FALSE], fc_cap), -fc_cap)
  hc <- hclust(dist(t(M)), method = "ward.D2")
  grp <- setNames(a$label, a$Sample_ID)
  lv <- unique(a$label[order(a$group)]); grp <- factor(grp, levels = lv)
  k <- length(lv)
  sizes <- table(grp)
  small <- names(sizes)[sizes < min_group_flag]
  size_rows[[key]] <- data.frame(subtype = st, gene_set = gs, k = k,
                                 group_sizes = paste0("n = ", paste(paste0(names(sizes), ":", sizes), collapse = " | ")),
                                 small_groups = paste(small, collapse = ", "))
  cat(sprintf("%-16s %-8s k = %d | %s%s\n", st, gs, k, paste(paste0(names(sizes), "=", sizes), collapse = " "),
              if (length(small)) paste0("  (small: ", paste(small, collapse = ", "), ")") else ""))

  # ---- Group profiles: mean log2FC of up- and down-DEGs per group ----
  for (g in lv) {
    ids <- names(grp)[grp == g]
    prof_rows[[paste(key, g)]] <- data.frame(subtype = st, gene_set = gs, group = g, n = length(ids),
      mean_log2FC_upDEGs = if (length(ob$up)) round(mean(ob$fc_asd[ob$up, ids]), 3) else NA,
      mean_log2FC_downDEGs = if (length(ob$down)) round(mean(ob$fc_asd[ob$down, ids]), 3) else NA,
      mean_log2FC_heatmap_genes = round(mean(ob$fc_asd[ob$sets[[gs]], ids]), 3),
      mean_abs_log2FC_heatmap_genes = round(mean(abs(ob$fc_asd[ob$sets[[gs]], ids])), 3))
  }

  # ---- ADI-R by group ----
  cl <- clin[names(grp), ]; cl$group <- grp
  for (adi in adi_cols) {
    s <- cl[!is.na(cl[[adi]]), ]
    ng <- table(droplevels(s$group))
    if (length(ng) < 2) next
    if (length(ng) == 2) {
      tst <- "Wilcoxon"; p <- wilcox.test(s[[adi]] ~ droplevels(s$group), exact = FALSE)$p.value
    } else {
      tst <- "Kruskal-Wallis"; p <- kruskal.test(s[[adi]] ~ droplevels(s$group))$p.value
      dunn_rows[[paste(key, adi)]] <- cbind(subtype = st, gene_set = gs, ADI = adi, dunn(s[[adi]], droplevels(s$group)))
    }
    adi_rows[[paste(key, adi)]] <- data.frame(subtype = st, gene_set = gs, ADI = adi, test = tst,
      n_by_group = paste(paste0(names(ng), ":", ng), collapse = " | "),
      medians = paste(paste0(names(ng), ":", tapply(s[[adi]], droplevels(s$group), median)), collapse = " | "), p = p)
  }

  # ---- Conditions by group (Fisher, simulated p) ----
  for (cnd in conditions) {
    tb <- table(cl$group, cl[[cnd]])
    if (ncol(tb) < 2 || any(colSums(tb) == 0)) next
    p <- fisher.test(tb, simulate.p.value = TRUE, B = n_sim)$p.value
    cond_rows[[paste(key, cnd)]] <- data.frame(subtype = st, gene_set = gs, condition = cnd,
      Y_by_group = paste(paste0(rownames(tb), ":", tb[, "Y"], "/", rowSums(tb)), collapse = " | "), p = p)
  }

  plot_adi[[key]] <- do.call(rbind, lapply(adi_cols, function(adi)
    data.frame(subtype = st, gene_set = gs, Sample_ID = rownames(cl), group = cl$group, ADI = adi, value = cl[[adi]])))
  plot_cond[[key]] <- do.call(rbind, lapply(conditions, function(cnd)
    data.frame(subtype = st, gene_set = gs, group = cl$group, condition = cnd, status = cl[[cnd]])))

  # ---- Final heatmap ----
  ann <- clin[colnames(M), c(rev(adi_cols), rev(conditions))]
  ann$Group <- grp[colnames(M)]
  ann <- ann[, c("Group", setdiff(colnames(ann), "Group"))]
  acol <- c(setNames(lapply(conditions, function(x) c(N = "grey92", Y = "grey25")), conditions),
            list(Group = setNames(grp_cols[seq_len(k)], lv)))
  # Draw the tree cut only when the groups match a cut of the tree (true for unedited assignments)
  same_cut <- nrow(unique(data.frame(cutree(hc, k), as.character(grp[hc$labels])))) == k
  title <- paste0(st, " - ", ifelse(gs == "DEG", paste0("top ", nrow(M), " MAST DEGs"), paste0("top ", nrow(M), " variable genes")),
                  "\nPer-donor log2FC vs Control mean (capped at +/-", fc_cap, "); ", k, " groups")
  for (ext in c("pdf", "png")) {
    pheatmap(M, color = heat_cols, breaks = brks, cluster_cols = hc, cutree_cols = if (same_cut) k else NA,
             clustering_method = "ward.D2", labels_row = sym[rownames(M)],
             fontsize_row = if (nrow(M) > 100) 3 else 5, fontsize_col = 6, fontsize = 7,
             annotation_col = ann, annotation_colors = acol, main = title, border_color = NA,
             filename = file.path(fd, paste0("heatmap_", fs, "_", gs, ".", ext)),
             width = 8, height = max(6, 0.04 * nrow(M) + 4))
  }
}

############################################
# Multiple testing and tables
############################################

R <- do.call(rbind, adi_rows); R$fdr <- p.adjust(R$p, "BH")
Cn <- do.call(rbind, cond_rows); Cn$fdr <- p.adjust(Cn$p, "BH")
write.csv(do.call(rbind, size_rows), file.path(fd, "group_sizes.csv"), row.names = FALSE)
write.csv(do.call(rbind, prof_rows), file.path(fd, "group_profiles.csv"), row.names = FALSE)
write.csv(R, file.path(fd, "ADIR_by_group.csv"), row.names = FALSE)
write.csv(Cn, file.path(fd, "conditions_by_group.csv"), row.names = FALSE)
if (length(dunn_rows)) write.csv(do.call(rbind, dunn_rows), file.path(fd, "ADIR_by_group_Dunn_posthoc.csv"), row.names = FALSE)

cat("\nADI-R by group (nominal p < 0.05):\n"); print(R[R$p < 0.05, c("subtype", "gene_set", "ADI", "test", "medians", "p", "fdr")], row.names = FALSE)
cat("\nConditions by group (nominal p < 0.05):\n"); print(Cn[Cn$p < 0.05, ], row.names = FALSE)
cat("\nTests: ADI-R", nrow(R), "| conditions", nrow(Cn), "| BH across all tests of each kind\n")

############################################
# Summary figures
############################################

PA <- do.call(rbind, plot_adi); PA <- PA[!is.na(PA$value), ]
PA$ADI <- factor(adi_labels[PA$ADI], levels = adi_labels)
R$ADI_lab <- factor(adi_labels[R$ADI], levels = adi_labels)
R$txt <- paste0("p = ", signif(R$p, 2))
PC <- do.call(rbind, plot_cond); PC <- PC[!is.na(PC$status), ]
Cn$txt <- paste0("p = ", signif(Cn$p, 2))

for (gs in unique(PA$gene_set)) {
  pa <- PA[PA$gene_set == gs, ]; ra <- R[R$gene_set == gs, ]
  p1 <- ggplot(pa, aes(group, value, color = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.12, size = 1) +
    geom_text(data = ra, aes(x = -Inf, y = Inf, label = txt), hjust = -0.1, vjust = 1.4, size = 2, inherit.aes = FALSE) +
    facet_grid(subtype ~ ADI, scales = "free") +
    scale_color_manual(values = grp_cols) +
    labs(title = paste0("ADI-R by patient subgroup (", gs, " heatmaps)"),
         caption = "Wilcoxon (2 groups) or Kruskal-Wallis (3+); FDR in ADIR_by_group.csv", x = NULL, y = "ADI-R score") +
    theme_bw(base_size = 7) + theme(legend.position = "none")
  ggsave(file.path(fd, paste0("ADIR_by_group_", gs, ".pdf")), p1, width = 7.5,
         height = 1.3 * length(unique(pa$subtype)) + 1, limitsize = FALSE)

  pc <- PC[PC$gene_set == gs, ]; rc <- Cn[Cn$gene_set == gs, ]
  p2 <- ggplot(pc, aes(group, fill = status)) + geom_bar(position = "fill", width = 0.7) +
    geom_text(data = rc, aes(x = -Inf, y = Inf, label = txt), hjust = -0.1, vjust = 1.3, size = 2, inherit.aes = FALSE) +
    facet_grid(subtype ~ condition) +
    scale_fill_manual(values = c(N = "grey85", Y = "grey25")) + scale_y_continuous(labels = scales::percent) +
    labs(title = paste0("Co-occurring conditions by patient subgroup (", gs, " heatmaps)"),
         caption = paste0("Fisher's exact test, ", n_sim, " simulated p-values; FDR in conditions_by_group.csv"),
         x = NULL, y = "Donors", fill = NULL) +
    theme_bw(base_size = 7)
  ggsave(file.path(fd, paste0("conditions_by_group_", gs, ".pdf")), p2, width = 8,
         height = 1.3 * length(unique(pc$subtype)) + 1, limitsize = FALSE)
}

writeLines(c(
  "Groups: group_assignments.csv `group` column (default: k with the best mean silhouette, 20a); labels from `group_name` if filled",
  paste0("Heatmaps: per-donor log2FC vs Control mean, capped at +/-", fc_cap, "; donors ordered by the 20a Ward D2 tree"),
  "ADI-R by group: Wilcoxon (2 groups) or Kruskal-Wallis + Dunn (tie-corrected, BH within comparison) for 3+; ASD donors with ADI-R only",
  paste0("Conditions by group: Fisher's exact test, simulated p (", n_sim, " replicates, seed ", seed, ")"),
  "BH FDR across all ADI-R tests, and separately across all condition tests (subtypes x gene sets)",
  paste0("Groups with fewer than ", min_group_flag, " donors are listed in group_sizes.csv")),
  file.path(fd, "clinical_20b_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(fd, "sessionInfo_20b.txt"))
cat("\nSaved", fd, "/\n")
