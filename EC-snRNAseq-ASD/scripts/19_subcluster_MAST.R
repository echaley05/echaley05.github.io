############################################
# 19_subcluster_MAST.R - MAST DEGs (ASD vs Control) within each neuronal subcluster
# Subclusters: sub_cluster_0.3 from 10_subcluster.py (excitatory 0-13, inhibitory 0-8)
# Same model and settings as 14b_mast.R; adds a per-subcluster breakdown with caution flags
# Run in R from the folder containing 12_annotated.rds and gene_table_biotype.csv
# Packages: Seurat, Matrix, MAST, data.table, ggplot2, ggrepel
# Runtime: roughly the same as 14b in total (the neurons are split into 23 smaller runs)
############################################

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(MAST); library(data.table); library(ggplot2); library(ggrepel)
})

############################################
# Settings
############################################

out          <- "19_subcluster_MAST"
min_pct      <- 0.10        # gene detected in >= 10% of nuclei in Control OR ASD
fdr          <- 0.05
lfc_report   <- 0.25        # DEG counts and top-gene lists use FDR < 0.05 and |log2FC| >= 0.25
model        <- ~ cngeneson + Age_z + PMI_z + Dx
coef_name    <- "DxASD"
n_top        <- 10
n_labels     <- 20          # volcano: label top 20 significant protein-coding/lncRNA genes with a symbol
label_biotypes <- c("protein_coding", "lncRNA")
dpi          <- 600
run_only     <- NULL        # NULL = all subclusters; or e.g. c("Inh_3", "Exc_0")

# Caution flags (reported only; every subcluster is still tested)
min_nuclei_donor <- 10      # a donor "counts" if it has >= 10 nuclei in the subcluster
min_donors_flag  <- 3       # flag if fewer than 3 such donors in either Dx group
skew_flag        <- 75      # flag if more than 75% of nuclei come from one Dx group
top_donor_flag   <- 50      # flag if one donor contributes more than 50% of nuclei

classes <- c(ExcitatoryNeuron = "Exc", InhibitoryNeuron = "Inh")
subtype_levels <- c("Ex-UL", "Ex-DL", "PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+")
dx_cols <- c(Control = "#4DBBD5", ASD = "#E64B35")

dir.create(file.path(out, "DEG"), recursive = TRUE, showWarnings = FALSE)

############################################
# 1. Load and define subclusters
############################################

seu <- readRDS("12_annotated.rds")
counts <- LayerData(seu, assay = "RNA", layer = "counts")
md <- seu@meta.data
rm(seu); gc()
md$Dx <- factor(md$Dx, levels = c("Control", "ASD"))
if (anyNA(md[, c("Age", "PMI")])) stop("Age or PMI missing")

gt <- read.csv("gene_table_biotype.csv")
gt <- gt[match(rownames(counts), gt$gene_id), c("gene_id", "gene_name", "gene_biotype")]
stopifnot(identical(gt$gene_id, rownames(counts)))

sub_raw <- as.character(md$sub_cluster_0.3)
neu <- md$cell_type %in% names(classes) & !is.na(sub_raw) & sub_raw != ""
md$subcluster <- NA_character_
md$subcluster[neu] <- paste0(classes[as.character(md$cell_type[neu])], "_", sub_raw[neu])

# One subtype per subcluster (subtypes were assigned by subcluster)
sc_tab <- unique(md[neu, c("subcluster", "cell_type", "subtype")])
stopifnot(!any(duplicated(sc_tab$subcluster)))
sc_tab$num <- as.integer(sub("^.*_", "", sc_tab$subcluster))
sc_tab <- sc_tab[order(match(sc_tab$cell_type, names(classes)), sc_tab$num), ]
sc_tab$label <- paste0(sc_tab$subcluster, " (", sc_tab$subtype, ")")
sc_tab$file <- gsub("[^A-Za-z0-9_-]", "", paste0(sc_tab$subcluster, "_", sc_tab$subtype))
if (!is.null(run_only)) sc_tab <- sc_tab[sc_tab$subcluster %in% run_only, ]
cat("Subclusters:", nrow(sc_tab), "\n"); print(sc_tab[, c("subcluster", "subtype")], row.names = FALSE)

############################################
# 2. Per subcluster: breakdown + MAST
############################################

rows <- list()

for (i in seq_len(nrow(sc_tab))) {
  sc <- sc_tab$subcluster[i]; lab <- sc_tab$label[i]
  t0 <- Sys.time()
  cat("\n====", lab, "====\n")
  cells <- rownames(md)[!is.na(md$subcluster) & md$subcluster == sc]
  cnt <- counts[, cells, drop = FALSE]
  cd <- md[cells, c("Sample_ID", "Dx", "Age", "PMI")]

  # ---- Breakdown ----
  n_dx <- table(cd$Dx)
  nd <- table(as.character(cd$Sample_ID))
  donor_dx <- tapply(as.character(cd$Dx), as.character(cd$Sample_ID), `[`, 1)
  big <- names(nd)[nd >= min_nuclei_donor]
  genes_per <- Matrix::colSums(cnt > 0); umis_per <- Matrix::colSums(cnt)
  b <- data.frame(
    subcluster = sc, label = lab, class = sc_tab$cell_type[i], subtype = sc_tab$subtype[i],
    nuclei = length(cells), nuclei_Control = n_dx[["Control"]], nuclei_ASD = n_dx[["ASD"]],
    pct_ASD = round(100 * n_dx[["ASD"]] / length(cells), 1),
    donors_Control = sum(donor_dx == "Control"), donors_ASD = sum(donor_dx == "ASD"),
    donors_ge10_Control = sum(donor_dx[big] == "Control"), donors_ge10_ASD = sum(donor_dx[big] == "ASD"),
    top_donor = names(nd)[which.max(nd)], top_donor_pct = round(100 * max(nd) / length(cells), 1),
    median_genes_Control = median(genes_per[cd$Dx == "Control"]), median_genes_ASD = median(genes_per[cd$Dx == "ASD"]),
    median_UMIs_Control = median(umis_per[cd$Dx == "Control"]), median_UMIs_ASD = median(umis_per[cd$Dx == "ASD"]))
  fl <- c(if (min(b$donors_ge10_Control, b$donors_ge10_ASD) < min_donors_flag) "few donors",
          if (max(b$pct_ASD, 100 - b$pct_ASD) > skew_flag) "Dx skew",
          if (b$top_donor_pct > top_donor_flag) "one-donor dominated")
  b$flags <- if (length(fl)) paste(fl, collapse = "; ") else ""
  cat("Nuclei: Control", b$nuclei_Control, "| ASD", b$nuclei_ASD, "(", b$pct_ASD, "% ASD )",
      "| donors >=10 nuclei: Control", b$donors_ge10_Control, "ASD", b$donors_ge10_ASD,
      if (nzchar(b$flags)) paste0("| FLAGS: ", b$flags) else "", "\n")

  res <- data.frame(genes_tested = NA, DEGs_up = NA, DEGs_down = NA, DEGs_up_protein_coding = NA,
                    DEGs_down_protein_coding = NA, DEGs_lncRNA = NA, top_up = "", top_down = "",
                    status = "", minutes = NA)

  if (any(n_dx == 0)) {
    res$status <- "not tested: one Dx group has no nuclei"
    cat(res$status, "\n")
  } else {
    # ---- Genes and expression ----
    det <- cnt > 0
    pct_c <- Matrix::rowMeans(det[, cd$Dx == "Control", drop = FALSE])
    pct_a <- Matrix::rowMeans(det[, cd$Dx == "ASD", drop = FALSE])
    keep <- pct_c >= min_pct | pct_a >= min_pct
    rm(det)
    x <- cnt[keep, , drop = FALSE] %*% Diagonal(x = 1e6 / umis_per)
    x@x <- log2(x@x + 1)
    logcpm <- as.matrix(x); dimnames(logcpm) <- list(rownames(cnt)[keep], cells); rm(x)

    cdat <- data.frame(wellKey = cells, Dx = cd$Dx,
                       Age_z = as.numeric(scale(cd$Age)), PMI_z = as.numeric(scale(cd$PMI)))
    cdat$Age_z[is.na(cdat$Age_z)] <- 0; cdat$PMI_z[is.na(cdat$PMI_z)] <- 0
    sca <- FromMatrix(logcpm, cData = cdat, fData = data.frame(primerid = rownames(logcpm)), check_sanity = FALSE)
    colData(sca)$cngeneson <- as.numeric(scale(colSums(assay(sca) > 0)))
    rm(logcpm); gc()

    tab <- tryCatch({
      zfit <- zlm(model, sca, method = "bayesglm", ebayes = TRUE)
      sm <- summary(zfit, doLRT = coef_name)$datatable
      p  <- sm[contrast == coef_name & component == "H", .(primerid, pvalue = `Pr(>Chisq)`)]
      fc <- sm[contrast == coef_name & component == "logFC", .(primerid, log2FC = coef, ci.hi, ci.lo)]
      as.data.frame(merge(p, fc, by = "primerid"))
    }, error = function(e) { cat("MAST failed:", conditionMessage(e), "\n"); NULL })

    if (is.null(tab)) {
      res$status <- "MAST failed (see console)"
    } else {
      tab$padj <- p.adjust(tab$pvalue, "BH")
      tab <- data.frame(gt[match(tab$primerid, gt$gene_id), ], tab[, -1],
                        pct_Control = round(100 * pct_c[tab$primerid], 1),
                        pct_ASD = round(100 * pct_a[tab$primerid], 1), row.names = NULL)
      tab$subcluster <- lab
      tab <- tab[order(tab$padj, tab$pvalue), ]
      write.csv(tab, file.path(out, "DEG", paste0("DEG_MAST_", sc_tab$file[i], ".csv")), row.names = FALSE)

      s <- !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$log2FC)
      up <- s & tab$log2FC >= lfc_report; dn <- s & tab$log2FC <= -lfc_report
      named <- !grepl("^ENSG", tab$gene_name)
      res$genes_tested <- nrow(tab)
      res$DEGs_up <- sum(up); res$DEGs_down <- sum(dn)
      res$DEGs_up_protein_coding <- sum(up & tab$gene_biotype == "protein_coding")
      res$DEGs_down_protein_coding <- sum(dn & tab$gene_biotype == "protein_coding")
      res$DEGs_lncRNA <- sum((up | dn) & tab$gene_biotype == "lncRNA")
      res$top_up   <- paste(head(tab$gene_name[up & named], n_top), collapse = ", ")
      res$top_down <- paste(head(tab$gene_name[dn & named], n_top), collapse = ", ")
      res$status <- "tested"
      cat("Genes tested:", nrow(tab), "| DEGs up", res$DEGs_up, "| down", res$DEGs_down, "\n")
      if (nzchar(res$top_up))   cat("  Top up:  ", res$top_up, "\n")
      if (nzchar(res$top_down)) cat("  Top down:", res$top_down, "\n")
    }
    rm(sca); gc()
  }
  res$minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  rows[[sc]] <- cbind(b, res)
  write.csv(do.call(rbind, rows), file.path(out, "subcluster_breakdown_and_DEGs.csv"), row.names = FALSE)
}

############################################
# 3. Summary figures
############################################

S <- do.call(rbind, rows)
S$subtype <- factor(S$subtype, levels = subtype_levels)
S$class <- factor(ifelse(S$class == "ExcitatoryNeuron", "Excitatory", "Inhibitory"), levels = c("Excitatory", "Inhibitory"))
lf <- paste0(S$label, ifelse(nzchar(S$flags), " !", ""))
S$label_flag <- factor(lf, levels = rev(lf))

nuc <- rbind(data.frame(label = S$label_flag, class = S$class, Dx = "Control", n = S$nuclei_Control),
             data.frame(label = S$label_flag, class = S$class, Dx = "ASD", n = S$nuclei_ASD))
nuc$Dx <- factor(nuc$Dx, levels = c("Control", "ASD"))
p1 <- ggplot(nuc, aes(n, label, fill = Dx)) + geom_col(width = 0.75) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = dx_cols) +
  labs(title = "Nuclei per neuronal subcluster by diagnosis", x = "Nuclei", y = NULL, fill = NULL,
       caption = paste0("! = flagged (few donors: < ", min_donors_flag, " donors with >= ", min_nuclei_donor,
                        " nuclei in a group; Dx skew: > ", skew_flag, "% from one group; one donor > ", top_donor_flag, "%)")) +
  theme_classic(base_size = 8) + theme(legend.position = "top", plot.caption = element_text(size = 6, color = "grey40"))

deg <- rbind(data.frame(label = S$label_flag, class = S$class, dir = "Up in ASD", n = S$DEGs_up),
             data.frame(label = S$label_flag, class = S$class, dir = "Down in ASD", n = -S$DEGs_down))
deg$dir <- factor(deg$dir, levels = c("Up in ASD", "Down in ASD"))
p2 <- ggplot(deg[!is.na(deg$n), ], aes(n, label, fill = dir)) + geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Up in ASD" = "#E64B35", "Down in ASD" = "#4DBBD5")) +
  scale_x_continuous(labels = abs) +
  labs(title = "MAST DEGs per neuronal subcluster", x = "DEGs (down | up)", y = NULL, fill = NULL,
       caption = paste0("DEGs: FDR < ", fdr, " and |log2FC| >= ", lfc_report, "; model ", deparse(model), "; ! = flagged")) +
  theme_classic(base_size = 8) + theme(legend.position = "top", plot.caption = element_text(size = 6, color = "grey40"))

h <- 0.22 * nrow(S) + 2
pdf(file.path(out, "subcluster_summary.pdf"), width = 6.5, height = h)
print(p1); print(p2)
dev.off()
ggsave(file.path(out, "subcluster_nuclei_by_Dx.png"), p1, width = 6.5, height = h, dpi = 600)
ggsave(file.path(out, "subcluster_DEG_counts.png"), p2, width = 6.5, height = h, dpi = 600)

cat("\nSummary:\n")
print(S[, c("label", "nuclei_Control", "nuclei_ASD", "donors_ge10_Control", "donors_ge10_ASD",
            "DEGs_up", "DEGs_down", "flags")], row.names = FALSE)

############################################
# 4. Volcano plots (same style as 14c_volcano.R: square-root y axis, p = 0 as triangles)
############################################

lfc_cut <- lfc_report
cols <- c("Up in ASD" = "#E64B35", "Down in ASD" = "#4DBBD5", "NS" = "grey80")

# Tick positions for a square-root y axis (real -log10 p values)
sqrt_breaks <- function(lim) {
  top <- lim[2]
  cand <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300, 500)
  b <- cand[cand <= top]
  if (length(b) > 7) b <- b[c(1, seq(length(b) - 5, length(b)))]
  b
}

volcano <- function(tab, lfc_col, title) {
  tab <- tab[!is.na(tab$pvalue) & !is.na(tab$padj) & !is.na(tab[[lfc_col]]), ]
  tab$lfc <- tab[[lfc_col]]

  # p = 0 (numerical underflow): plot at the smallest non-zero p in this table
  tab$zero <- tab$pvalue == 0
  n_zero <- sum(tab$zero)
  if (n_zero > 0) tab$pvalue[tab$pvalue == 0] <- min(tab$pvalue[tab$pvalue > 0])
  tab$y <- -log10(tab$pvalue)

  tab$sig <- "NS"
  tab$sig[tab$padj < fdr & tab$lfc >=  lfc_cut] <- "Up in ASD"
  tab$sig[tab$padj < fdr & tab$lfc <= -lfc_cut] <- "Down in ASD"
  tab$sig <- factor(tab$sig, levels = names(cols))

  # FDR line: p-value of the least significant gene with padj < fdr
  fdr_y <- if (any(tab$padj < fdr)) -log10(max(tab$pvalue[tab$padj < fdr])) else NA

  # Labels: significant protein-coding/lncRNA genes with a real symbol (no Ensembl IDs)
  has_symbol <- !is.na(tab$gene_name) & tab$gene_name != "" & !grepl("^ENSG", tab$gene_name)
  lab <- tab[tab$sig != "NS" & tab$gene_biotype %in% label_biotypes & has_symbol, ]
  lab <- head(lab[order(lab$padj, lab$pvalue, -abs(lab$lfc)), ], n_labels)

  n_up <- sum(tab$sig == "Up in ASD"); n_dn <- sum(tab$sig == "Down in ASD")
  xr <- max(abs(tab$lfc)) * 1.05
  ytop <- max(tab$y) * 1.45                     # headroom so labels can sit above the top genes

  tab <- tab[order(tab$sig != "NS"), ]          # significant points drawn on top

  p <- ggplot(tab, aes(lfc, y)) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    { if (!is.na(fdr_y)) geom_hline(yintercept = fdr_y, linetype = "dashed", linewidth = 0.3, color = "grey40") } +
    geom_point(data = tab[!tab$zero, ], aes(color = sig), size = 0.7, alpha = 0.8, stroke = 0) +
    geom_point(data = tab[tab$zero, ], aes(color = sig), size = 1.4, shape = 17) +
    geom_text_repel(data = lab, aes(label = gene_name), size = 2.3, fontface = "italic",
                    max.overlaps = Inf, segment.size = 0.2, segment.color = "grey50",
                    min.segment.length = 0, box.padding = 0.4, point.padding = 0.15,
                    force = 3, force_pull = 0.5, max.time = 3, max.iter = 1e5, seed = 42) +
    scale_color_manual(values = cols, drop = FALSE) +
    scale_x_continuous(limits = c(-xr, xr)) +
    scale_y_continuous(trans = "sqrt", limits = c(0, ytop), breaks = sqrt_breaks,
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, subtitle = paste0("Down: ", n_dn, "   |   Up: ", n_up),
         x = expression(log[2]~"fold change (ASD vs Control)"),
         y = expression(-log[10]~italic(p)~"(square-root scale)"), color = NULL,
         caption = paste0("Significant: FDR < ", fdr, " and |log2FC| >= ", lfc_cut,
                          "; labels: top ", n_labels, " protein-coding/lncRNA genes with a symbol",
                          if (n_zero > 0) paste0("\nTriangles: ", n_zero, " genes with p below the numerical limit (reported as 0), plotted at the smallest non-zero p") else "")) +
    theme_classic(base_size = 8) +
    theme(legend.position = "top", legend.key.size = unit(3, "mm"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 7),
          plot.caption = element_text(size = 5.5, color = "grey40"),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3)) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1, shape = 16)))

  list(plot = p, n_up = n_up, n_down = n_dn, n_zero = n_zero, max_y = max(tab$y))
}

vdir <- file.path(out, "volcano")
dir.create(vdir, showWarnings = FALSE)
vrows <- list()
pdf(file.path(vdir, "volcano_all_subclusters.pdf"), width = 4.5, height = 4.5)
for (i in seq_len(nrow(sc_tab))) {
  f <- file.path(out, "DEG", paste0("DEG_MAST_", sc_tab$file[i], ".csv"))
  if (!file.exists(f)) next
  tab <- read.csv(f, check.names = FALSE)
  if (sum(!is.na(tab$padj)) < 10) next
  v <- volcano(tab, "log2FC", paste0(sc_tab$label[i], " - MAST"))
  print(v$plot)
  ggsave(file.path(vdir, paste0("volcano_", sc_tab$file[i], ".pdf")), v$plot, width = 4.5, height = 4.5)
  ggsave(file.path(vdir, paste0("volcano_", sc_tab$file[i], ".png")), v$plot, width = 4.5, height = 4.5, dpi = dpi)
  vrows[[i]] <- data.frame(subcluster = sc_tab$label[i], up = v$n_up, down = v$n_down,
                           p_zero = v$n_zero, max_neglog10p = round(v$max_y, 1))
}
dev.off()
if (length(vrows)) write.csv(do.call(rbind, vrows), file.path(vdir, "volcano_counts.csv"), row.names = FALSE)
cat("Volcano plots saved in", vdir, "\n")

writeLines(c(
  "Subclusters: sub_cluster_0.3 within ExcitatoryNeuron (Exc_*) and InhibitoryNeuron (Inh_*); parent subtype in brackets",
  paste0("Gene filter: all biotypes; detected in >= ", 100 * min_pct, "% of nuclei in Control or ASD (per subcluster)"),
  "Expression: log2(CPM + 1), library size = total counts per nucleus",
  paste("Model: MAST zlm", deparse(model), "; bayesglm, ebayes = TRUE; fixed effects only"),
  "cngeneson: scaled number of genes detected per nucleus after the gene filter; Age and PMI scaled within subcluster",
  paste0("Test: hurdle LRT for ", coef_name, "; BH within subcluster; DEGs counted at FDR < ", fdr, " and |log2FC| >= ", lfc_report),
  paste0("Flags (reported only): < ", min_donors_flag, " donors with >= ", min_nuclei_donor, " nuclei in either group; > ",
         skew_flag, "% nuclei from one Dx group; top donor > ", top_donor_flag, "% of nuclei"),
  "Top genes: first 10 by FDR among DEGs, genes without a symbol skipped",
  "Caveat: nuclei treated as independent (no donor term), as in 14b"),
  file.path(out, "subcluster_MAST_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_19.txt"))
cat("\nSaved", out, "/\n")
