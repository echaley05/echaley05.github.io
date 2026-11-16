############################################
# 18_age_by_dx.R - Age x Diagnosis interaction per subtype (MAST, fixed effects)
# 1. Age coverage check by Dx and subtype
# 2. MAST hurdle model ~ cngeneson + PMI + Age * Dx (ebayes = FALSE; required for interaction models)
#    LRT on Age:Dx; age slopes in Control and ASD (log2FC per decade); Dx effect at the mean age
# 3. Donor-level plots of the top interaction genes
# Run in R from the folder containing 12_annotated.rds and gene_table_biotype.csv
# Packages: Seurat, Matrix, MAST, data.table, ggplot2
############################################

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(MAST); library(data.table); library(ggplot2)
})

############################################
# Settings
############################################

out          <- "18_age_by_dx"
group_by     <- "subtype"
run_subtypes <- c("Ex-UL", "Ex-DL", "PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+",
                  "Astrocyte", "Oligodendrocyte", "OPC", "Microglia")   # Endothelial skipped (too few nuclei)
min_pct      <- 0.10       # gene detected in >= 10% of nuclei in Control OR ASD (same as 14b)
fdr          <- 0.05       # interaction genes: FDR < 0.05 on Age:Dx (no slope cutoff)
age_unit     <- 10         # Age scaled per decade
min_nuclei_donor <- 10     # for coverage summary and donor-level plots
n_plot       <- 12         # top interaction genes plotted per subtype
dir.create(out, showWarnings = FALSE)

dx_cols <- c(Control = "#4DBBD5", ASD = "#E64B35")

############################################
# Load
############################################

seu <- readRDS("12_annotated.rds")
counts <- LayerData(seu, assay = "RNA", layer = "counts")
md <- seu@meta.data
md$Dx <- factor(md$Dx, levels = c("Control", "ASD"))
gt <- read.csv("gene_table_biotype.csv")
gt <- gt[match(rownames(counts), gt$gene_id), c("gene_id", "gene_name", "gene_biotype")]
rm(seu); gc()

############################################
# 1. Age coverage
############################################

donors <- unique(md[, c("Sample_ID", "Dx", "Age", "PMI")])
stopifnot(!any(duplicated(donors$Sample_ID)), !anyNA(donors$Age), !anyNA(donors$PMI))
wt <- wilcox.test(Age ~ Dx, data = donors, exact = FALSE)
cat("Donor age by Dx:\n")
print(aggregate(Age ~ Dx, donors, function(x) c(n = length(x), min = min(x), median = median(x), max = max(x))))
cat("Wilcoxon p (age differs by Dx):", signif(wt$p.value, 3), "\n")

cov_rows <- list()
for (st in run_subtypes) {
  n <- table(md$Sample_ID[md[[group_by]] == st])
  d <- donors[donors$Sample_ID %in% names(n)[n >= min_nuclei_donor], ]
  rc <- range(d$Age[d$Dx == "Control"]); ra <- range(d$Age[d$Dx == "ASD"])
  cov_rows[[st]] <- data.frame(
    subtype = st, donors_Control = sum(d$Dx == "Control"), donors_ASD = sum(d$Dx == "ASD"),
    age_Control = paste(rc, collapse = "-"), age_ASD = paste(ra, collapse = "-"),
    overlap = paste(max(rc[1], ra[1]), min(rc[2], ra[2]), sep = "-"),
    overlap_fraction = round(max(0, min(rc[2], ra[2]) - max(rc[1], ra[1])) / (max(rc[2], ra[2]) - min(rc[1], ra[1])), 2))
}
coverage <- do.call(rbind, cov_rows)
write.csv(donors, file.path(out, "donor_age.csv"), row.names = FALSE)
write.csv(coverage, file.path(out, "age_coverage_by_subtype.csv"), row.names = FALSE)
cat("\nAge coverage by subtype (donors with >=", min_nuclei_donor, "nuclei):\n"); print(coverage, row.names = FALSE)

p_age <- ggplot(donors, aes(Dx, Age, color = Dx)) +
  geom_jitter(width = 0.12, height = 0, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.3) +
  scale_color_manual(values = dx_cols) +
  labs(title = "Donor age by diagnosis", subtitle = paste0("Wilcoxon p = ", signif(wt$p.value, 2)),
       x = NULL, y = "Age (years)") +
  theme_classic(base_size = 9) + theme(legend.position = "none")
ggsave(file.path(out, "donor_age_by_Dx.pdf"), p_age, width = 3, height = 3.5)
ggsave(file.path(out, "donor_age_by_Dx.png"), p_age, width = 3, height = 3.5, dpi = 600)

############################################
# 2. MAST interaction model per subtype
############################################

summary_rows <- list()

for (st in run_subtypes) {
  t0 <- Sys.time()
  cat("\n====", st, "====\n")
  cells <- rownames(md)[md[[group_by]] == st]
  cnt <- counts[, cells]
  cd <- md[cells, c("Sample_ID", "Dx", "Age", "PMI")]

  det <- cnt > 0
  pct_c <- Matrix::rowMeans(det[, cd$Dx == "Control", drop = FALSE])
  pct_a <- Matrix::rowMeans(det[, cd$Dx == "ASD", drop = FALSE])
  keep <- pct_c >= min_pct | pct_a >= min_pct
  rm(det)
  cat("Nuclei: Control", sum(cd$Dx == "Control"), "| ASD", sum(cd$Dx == "ASD"),
      "| genes passing filter:", sum(keep), "\n")

  lib <- Matrix::colSums(cnt)
  x <- cnt[keep, ] %*% Diagonal(x = 1e6 / lib)
  x@x <- log2(x@x + 1)
  logcpm <- as.matrix(x); dimnames(logcpm) <- list(rownames(cnt)[keep], cells)
  rm(x, cnt); gc()

  # Age centred on the mean DONOR age in this subtype, in decades
  donor_age_mean <- mean(unique(cd[, c("Sample_ID", "Age")])$Age)
  cdat <- data.frame(wellKey = cells, Dx = cd$Dx,
                     Age_c = (cd$Age - donor_age_mean) / age_unit,
                     PMI_z = as.numeric(scale(cd$PMI)))
  sca <- FromMatrix(logcpm, cData = cdat, fData = data.frame(primerid = rownames(logcpm)), check_sanity = FALSE)
  colData(sca)$cngeneson <- as.numeric(scale(colSums(assay(sca) > 0)))

  zfit <- zlm(~ cngeneson + PMI_z + Age_c * Dx, sca, method = "bayesglm", ebayes = FALSE)

  # Likelihood-ratio tests (hurdle = combined discrete + continuous)
  # CoefficientHypothesis tests the single coefficient (a term name is refused because of marginality)
  lrt_int <- lrTest(zfit, CoefficientHypothesis("Age_c:DxASD"))
  lrt_dx  <- lrTest(zfit, CoefficientHypothesis("DxASD"))

  # Combined log2FC from the fitted coefficients, computed directly (MAST::logFC fails with these contrasts).
  # Expected log2(CPM+1) = P(detected) * mean when detected = plogis(D . x) * (C . x),
  # evaluated at mean CDR and PMI (both scaled, so 0).
  Cc <- coef(zfit, "C"); Dc <- coef(zfit, "D")
  xvec <- function(age, asd) {
    v <- setNames(rep(0, ncol(Cc)), colnames(Cc))
    v["(Intercept)"] <- 1; v["Age_c"] <- age; v["DxASD"] <- asd; v["Age_c:DxASD"] <- age * asd
    v
  }
  mu <- function(v) plogis(as.vector(Dc[, names(v)] %*% v)) * as.vector(Cc[, names(v)] %*% v)
  m00 <- mu(xvec(0, 0)); m10 <- mu(xvec(1, 0)); m01 <- mu(xvec(0, 1)); m11 <- mu(xvec(1, 1))
  slope_ctrl <- m10 - m00          # change per decade, Control
  slope_asd  <- m11 - m01          # change per decade, ASD
  dx_mean    <- m01 - m00          # ASD - Control at the mean age

  tab <- data.frame(
    gene_id = rownames(logcpm),
    gene_name = gt$gene_name[match(rownames(logcpm), gt$gene_id)],
    gene_biotype = gt$gene_biotype[match(rownames(logcpm), gt$gene_id)],
    slope_Control_per_decade = slope_ctrl,
    slope_ASD_per_decade = slope_asd,
    interaction_slope_diff = slope_asd - slope_ctrl,
    p_interaction = lrt_int[, "hurdle", "Pr(>Chisq)"],
    Dx_log2FC_at_mean_age = dx_mean,
    p_Dx_at_mean_age = lrt_dx[, "hurdle", "Pr(>Chisq)"],
    pct_Control = round(100 * pct_c[rownames(logcpm)], 1),
    pct_ASD = round(100 * pct_a[rownames(logcpm)], 1))
  tab$padj_interaction <- p.adjust(tab$p_interaction, "BH")
  tab$padj_Dx_at_mean_age <- p.adjust(tab$p_Dx_at_mean_age, "BH")

  sig <- !is.na(tab$padj_interaction) & tab$padj_interaction < fdr
  tab$pattern <- ""
  s_c <- tab$slope_Control_per_decade; s_a <- tab$slope_ASD_per_decade
  tab$pattern[sig] <- ifelse(sign(s_c[sig]) != sign(s_a[sig]),
                             ifelse(s_a[sig] > 0, "Opposite: rises in ASD, falls in Control",
                                    "Opposite: falls in ASD, rises in Control"),
                             ifelse(abs(s_a[sig]) > abs(s_c[sig]),
                                    ifelse(s_a[sig] > 0, "Same direction: rises faster in ASD", "Same direction: falls faster in ASD"),
                                    ifelse(s_c[sig] > 0, "Same direction: rises faster in Control", "Same direction: falls faster in Control")))
  tab <- tab[order(tab$padj_interaction, tab$p_interaction), ]
  fs <- gsub("[^A-Za-z0-9_-]", "", st)
  write.csv(tab, file.path(out, paste0("AgeDx_MAST_", fs, ".csv")), row.names = FALSE)

  pat <- table(factor(tab$pattern[tab$pattern != ""], levels = c(
    "Opposite: rises in ASD, falls in Control", "Opposite: falls in ASD, rises in Control",
    "Same direction: rises faster in ASD", "Same direction: falls faster in ASD",
    "Same direction: rises faster in Control", "Same direction: falls faster in Control")))
  summary_rows[[st]] <- data.frame(
    subtype = st, nuclei = length(cells), donors = length(unique(cd$Sample_ID)),
    genes_tested = nrow(tab), interaction_genes = sum(sig),
    protein_coding = sum(sig & tab$gene_biotype == "protein_coding"),
    lncRNA = sum(sig & tab$gene_biotype == "lncRNA"),
    Dx_at_mean_age_genes = sum(tab$padj_Dx_at_mean_age < fdr, na.rm = TRUE),
    t(setNames(as.vector(pat), names(pat))), check.names = FALSE,
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  cat("Interaction genes (FDR <", fdr, "):", sum(sig), "| time", summary_rows[[st]]$minutes, "min\n")
  if (sum(sig)) print(pat)

  # ---- 3. Donor-level plots of top interaction genes (symbol-labelled) ----
  top <- head(tab$gene_id[sig & !grepl("^ENSG", tab$gene_name)], n_plot)
  if (length(top)) {
    n_d <- table(cd$Sample_ID); ok_d <- names(n_d)[n_d >= min_nuclei_donor]
    sel <- cd$Sample_ID %in% ok_d
    dm <- rowsum(t(logcpm[top, sel, drop = FALSE]), as.character(cd$Sample_ID[sel]))
    dm <- dm / as.vector(table(as.character(cd$Sample_ID[sel]))[rownames(dm)])
    dd <- donors[match(rownames(dm), donors$Sample_ID), ]
    long <- do.call(rbind, lapply(top, function(g) data.frame(
      gene = tab$gene_name[tab$gene_id == g], Sample_ID = rownames(dm), Dx = dd$Dx, Age = dd$Age, expr = dm[, g])))
    long$gene <- factor(long$gene, levels = unique(long$gene))
    p <- ggplot(long, aes(Age, expr, color = Dx)) +
      geom_point(size = 1.3) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
      facet_wrap(~ gene, scales = "free_y", ncol = 4) +
      scale_color_manual(values = dx_cols) +
      labs(title = paste0(st, " - top Age x Dx interaction genes (FDR < ", fdr, ")"),
           subtitle = paste0("Points: donor mean log2(CPM+1) (donors with >= ", min_nuclei_donor,
                             " nuclei); lines: per-group linear fit for display"),
           x = "Age (years)", y = "Mean log2(CPM + 1)", color = NULL) +
      theme_classic(base_size = 8) +
      theme(strip.text = element_text(face = "italic"), legend.position = "top")
    ggsave(file.path(out, paste0("AgeDx_top_genes_", fs, ".pdf")), p, width = 8,
           height = 1.9 * ceiling(length(top) / 4) + 1)
  }

  rm(sca, zfit, lrt_int, lrt_dx, logcpm); gc()
  write.csv(do.call(rbind, summary_rows), file.path(out, "AgeDx_summary.csv"), row.names = FALSE)
}

############################################
# Summary and settings
############################################

S <- do.call(rbind, summary_rows)
cat("\nSummary:\n"); print(S, row.names = FALSE)

writeLines(c(
  paste("Subtypes:", paste(run_subtypes, collapse = ", ")),
  paste0("Gene filter: all biotypes; detected in >= ", 100 * min_pct, "% of nuclei in Control or ASD (per subtype)"),
  "Expression: log2(CPM + 1); library size = total counts per nucleus",
  "Model: MAST zlm ~ cngeneson + PMI_z + Age_c * Dx; method = bayesglm, ebayes = FALSE; fixed effects only",
  paste0("Age_c: (Age - mean donor age in the subtype) / ", age_unit, " (per decade); PMI_z: scaled within subtype"),
  "Tests: hurdle likelihood-ratio tests for Age_c:DxASD (interaction) and DxASD (Dx effect at the mean age); BH within subtype",
  "Slopes: change in expected log2(CPM+1) per decade = plogis(D.x)*(C.x) from MAST hurdle coefficients at mean CDR and PMI, in Control and ASD; interaction = ASD slope - Control slope",
  paste0("Interaction genes: padj < ", fdr, " (no slope cutoff); patterns from the signs and sizes of the two slopes"),
  paste0("Plots: donor mean log2(CPM+1) vs age (donors with >= ", min_nuclei_donor, " nuclei), per-group lm lines for display only"),
  paste("Donor age by Dx, Wilcoxon p =", signif(wt$p.value, 3))),
  file.path(out, "AgeDx_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_18.txt"))
cat("\nSaved", out, "/\n")
