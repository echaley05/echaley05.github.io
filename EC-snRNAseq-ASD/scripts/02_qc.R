############################################
# 02_qc.R — QC metrics, tables, and figures on ALL genes
# All 78,899 genes kept (biotypes annotated for reporting only)
# No nuclei, samples, or genes are removed
############################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

out <- "02_QC"
dir.create(out, showWarnings = FALSE)
save_tbl <- function(x, f) write.csv(x, file.path(out, f), row.names = FALSE)

seu <- readRDS("01_raw_loaded.rds")
gene_table <- read.csv("gene_table.csv")

############################################
# 1. Annotate biotypes (GENCODE v50 GTF) — annotation only, no filtering
############################################

gtf_file <- "gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"

# Read only the "gene" lines, 500k lines at a time
con <- gzfile(gtf_file, "r")
gene_lines <- character()
repeat {
  chunk <- readLines(con, n = 500000)
  if (length(chunk) == 0) break
  gene_lines <- c(gene_lines, chunk[grepl("\tgene\t", chunk, fixed = TRUE)])
}
close(con)

get_attr <- function(x, key) sub(paste0('.*', key, ' "([^"]*)".*'), "\\1", x)
gtf_genes <- data.frame(
  id   = sub("\\.[0-9]+", "", get_attr(gene_lines, "gene_id")),  # drop version, keep _PAR_Y
  name = get_attr(gene_lines, "gene_name"),
  type = get_attr(gene_lines, "gene_type"))
rm(gene_lines)
gtf_genes <- gtf_genes[!duplicated(gtf_genes$id), ]

i <- match(gene_table$gene_id, gtf_genes$id)
gene_table$gene_biotype <- ifelse(is.na(i), "not_in_GTF", gtf_genes$type[i])
cat("Genes not found in GTF:", sum(is.na(i)), "of", nrow(gene_table), "(kept, labelled not_in_GTF)\n")

write.csv(gene_table, "gene_table_biotype.csv", row.names = FALSE)

############################################
# 2. Biotype summary (genes and UMIs per biotype)
############################################

counts <- LayerData(seu, layer = "counts")
stopifnot(identical(rownames(counts), gene_table$gene_id))

biotype_summary <- gene_table %>%
  mutate(umis = Matrix::rowSums(counts)) %>%
  group_by(gene_biotype) %>%
  summarise(n_genes = n(), total_umis = sum(umis), .groups = "drop") %>%
  mutate(pct_umis = round(100 * total_umis / sum(total_umis), 2)) %>%
  arrange(desc(total_umis))
print(biotype_summary, n = 15)
save_tbl(biotype_summary, "biotype_summary.csv")
rm(counts); invisible(gc())

sym <- gene_table$gene_name

############################################
# 3. QC metrics (all genes)
############################################

is_mt <- grepl("^MT-", sym)
is_rb <- grepl("^RPS|^RPL", sym)
cat(sum(is_mt), "MT features\n")
cat(sum(is_rb), "RPS/RPL features\n")

seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = rownames(seu)[is_mt])
seu[["percent.rb"]] <- PercentageFeatureSet(seu, features = rownames(seu)[is_rb])
seu$log10GenesPerUMI <- ifelse(seu$nCount_RNA > 1,
                               log10(seu$nFeature_RNA) / log10(seu$nCount_RNA), NA)

saveRDS(seu, "02_qc_metrics.rds")

md <- seu@meta.data %>% mutate(Dx = factor(Dx, levels = c("Control", "ASD")))

############################################
# 4. Nucleus-level threshold tables
############################################

f_nfeat <- c(200, 300, 500, 750, 1000)
f_ncount <- c(500, 1000, 2000)
c_nfeat <- c(5000, 7500, 10000, 15000)
c_ncount <- c(50000, 100000, 200000)
c_mt <- c(1, 5, 10, 20)
c_rb <- c(1, 2, 5, 10)

rules <- c(
  setNames(lapply(f_nfeat,  function(t) md$nFeature_RNA <= t), paste0("nFeature<=", f_nfeat)),
  setNames(lapply(f_ncount, function(t) md$nCount_RNA <= t),   paste0("nCount<=", f_ncount)),
  setNames(lapply(c_nfeat,  function(t) md$nFeature_RNA > t),  paste0("nFeature>", c_nfeat)),
  setNames(lapply(c_ncount, function(t) md$nCount_RNA > t),    paste0("nCount>", c_ncount)),
  setNames(lapply(c_mt,     function(t) md$percent.mt > t),    paste0("percent.mt>", c_mt)),
  setNames(lapply(c_rb,     function(t) md$percent.rb > t),    paste0("percent.rb>", c_rb)))

fail_df <- data.frame(Sample_ID = md$Sample_ID, Dx = md$Dx,
                      do.call(cbind, rules), check.names = FALSE)

pct_lost <- ~ round(100 * mean(.x), 1)
sens_by_sample <- fail_df %>% group_by(Sample_ID, Dx) %>%
  summarise(n_nuclei = n(), across(all_of(names(rules)), pct_lost), .groups = "drop")
sens_by_dx <- fail_df %>% group_by(Dx) %>%
  summarise(n_nuclei = n(), across(all_of(names(rules)), pct_lost))

adaptive <- md %>% group_by(Sample_ID, Dx) %>% summarise(
  nCount_low    = 10^(median(log10(nCount_RNA))   - 3 * mad(log10(nCount_RNA))),
  nCount_high   = 10^(median(log10(nCount_RNA))   + 3 * mad(log10(nCount_RNA))),
  nFeature_low  = 10^(median(log10(nFeature_RNA)) - 3 * mad(log10(nFeature_RNA))),
  nFeature_high = 10^(median(log10(nFeature_RNA)) + 3 * mad(log10(nFeature_RNA))),
  percent.mt_high = median(percent.mt) + 3 * mad(percent.mt),
  percent.rb_high = median(percent.rb) + 3 * mad(percent.rb),
  .groups = "drop")

upper_tail <- md %>% group_by(Sample_ID, Dx) %>%
  summarise(across(c(nFeature_RNA, nCount_RNA),
                   list(q99 = ~ unname(quantile(.x, .99)),
                        q999 = ~ unname(quantile(.x, .999)),
                        max = max),
                   .names = "{.col}__{.fn}"),
            .groups = "drop")

save_tbl(sens_by_sample, "threshold_sensitivity_by_sample.csv")
save_tbl(sens_by_dx, "threshold_sensitivity_by_Dx.csv")
save_tbl(adaptive, "adaptive_bounds_by_sample.csv")
save_tbl(upper_tail, "upper_tail_by_sample.csv")

############################################
# 5. Sample-level table
############################################

sample_tbl <- md %>% group_by(Sample_ID, Dx, Age, PMI, RIN) %>% summarise(
  n_nuclei = n(),
  med_nCount = median(nCount_RNA),       iqr_nCount = IQR(nCount_RNA),
  med_nFeature = median(nFeature_RNA),   iqr_nFeature = IQR(nFeature_RNA),
  med_percent.mt = median(percent.mt),   iqr_percent.mt = IQR(percent.mt),
  med_percent.rb = median(percent.rb),   iqr_percent.rb = IQR(percent.rb),
  med_log10GenesPerUMI = median(log10GenesPerUMI, na.rm = TRUE),
  reads_per_UMI = sum(mread_count) / sum(tscp_count),
  .groups = "drop") %>%
  left_join(sens_by_sample %>%
              select(Sample_ID, starts_with("nFeature"), starts_with("nCount"),
                     starts_with("percent.mt")) %>%
              rename_with(~ paste0("pct_", .x), -Sample_ID),
            by = "Sample_ID")

fit <- MASS::rlm(log10(med_nFeature) ~ log10(med_nCount), data = sample_tbl)
sample_tbl$complexity_resid <- unname(resid(fit))
r <- sample_tbl$complexity_resid
sample_tbl$complexity_resid_mads <- round((r - median(r)) / mad(r), 1)

sample_levels <- sample_tbl %>% arrange(Dx, med_nCount) %>% pull(Sample_ID)
md$Sample_f <- factor(md$Sample_ID, levels = sample_levels)

sm_metrics <- c("med_nCount", "med_nFeature", "med_percent.mt", "med_percent.rb",
                "n_nuclei", "reads_per_UMI", "complexity_resid")

cov_tests <- do.call(rbind, lapply(sm_metrics, function(m) {
  x <- sample_tbl[[m]]
  o <- data.frame(metric = m, p_Dx = suppressWarnings(
    wilcox.test(x[sample_tbl$Dx == "ASD"], x[sample_tbl$Dx == "Control"])$p.value))
  for (v in c("Age", "PMI", "RIN")) {
    t <- suppressWarnings(cor.test(x, sample_tbl[[v]], method = "spearman"))
    o[[paste0("rho_", v)]] <- unname(t$estimate)
    o[[paste0("p_", v)]] <- t$p.value
  }
  o
}))

save_tbl(sample_tbl, "sample_qc_summary.csv")
save_tbl(cov_tests, "qc_vs_covariates_tests.csv")

############################################
# 6. Figures (one PDF)
############################################

dx_cols <- c(Control = "#4DBBD5", ASD = "#E64B35")
dx_col <- scale_color_manual(values = dx_cols)
dx_fill <- scale_fill_manual(values = dx_cols)
th <- theme_bw(base_size = 10)
th90 <- th + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
mt_x <- scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5),
                           breaks = c(0, 1, 5, 10, 20, 40, 80))

pdf(file.path(out, "QC_figures.pdf"), width = 16, height = 9)

# Per-sample distributions
for (m in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "log10GenesPerUMI")) {
  p <- ggplot(md, aes(Sample_f, .data[[m]], fill = Dx)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    stat_summary(fun = median, geom = "point", size = 0.7) +
    dx_fill + labs(title = m, x = NULL) + th90
  if (m %in% c("nCount_RNA", "nFeature_RNA")) p <- p + scale_y_log10()
  print(p)
}

# Retention curves
idx <- split(seq_len(nrow(md)), md$Sample_ID)
curves <- function(metric, grid, floor) {
  do.call(rbind, lapply(names(idx), function(s) {
    v <- md[[metric]][idx[[s]]]
    data.frame(Sample_ID = s, Dx = md$Dx[idx[[s]][1]], x = grid,
               retained = sapply(grid, function(g) if (floor) mean(v > g) else mean(v <= g)))
  }))
}
curve_plot <- function(d, xlab, trans) {
  ggplot(d, aes(x, retained, group = Sample_ID, color = Dx)) +
    geom_line(alpha = 0.7) + scale_x_continuous(trans = trans) +
    dx_col + labs(x = xlab, y = "Fraction of nuclei retained") + th
}
print(curve_plot(curves("nFeature_RNA", 10^seq(2, 4.3, length.out = 150), TRUE),
                 "nFeature_RNA floor", "log10"))
print(curve_plot(curves("nCount_RNA", 10^seq(2, 5.5, length.out = 150), TRUE),
                 "nCount_RNA floor", "log10"))
print(curve_plot(curves("nFeature_RNA", 10^seq(3, 4.5, length.out = 150), FALSE),
                 "nFeature_RNA ceiling", "log10"))
print(curve_plot(curves("nCount_RNA", 10^seq(4, 6, length.out = 150), FALSE),
                 "nCount_RNA ceiling", "log10"))
print(curve_plot(curves("percent.mt", seq(0, 100, length.out = 300), FALSE),
                 "percent.mt ceiling", "sqrt"))

# Threshold sensitivity heatmap
sens_long <- sens_by_sample %>%
  pivot_longer(-c(Sample_ID, Dx, n_nuclei), names_to = "rule", values_to = "pct_lost") %>%
  mutate(Sample_ID = factor(Sample_ID, sample_levels), rule = factor(rule, names(rules)))
print(ggplot(sens_long, aes(rule, Sample_ID, fill = pct_lost)) +
        geom_tile() + geom_text(aes(label = round(pct_lost)), size = 2) +
        scale_fill_gradient(low = "white", high = "#B2182B", name = "% nuclei lost") +
        labs(x = NULL, y = NULL) + th90)

# Upper tail
ut <- upper_tail %>%
  pivot_longer(-c(Sample_ID, Dx), names_to = c("metric", "stat"), names_sep = "__") %>%
  mutate(Sample_ID = factor(Sample_ID, sample_levels))
print(ggplot(ut, aes(Sample_ID, value, color = stat)) +
        geom_point() + scale_y_log10() +
        facet_grid(metric ~ Dx, scales = "free", space = "free_x") + labs(x = NULL) + th90)

# Barcode-rank plots
rank_df <- md %>% group_by(Sample_f) %>%
  arrange(desc(nCount_RNA), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>% ungroup()
print(ggplot(rank_df, aes(rank, nCount_RNA, color = Dx)) +
        geom_line(linewidth = 0.4) + scale_x_log10() + scale_y_log10() +
        facet_wrap(~ Sample_f, ncol = 8) + dx_col + th + theme(text = element_text(size = 7)))

# Relationship plots
print(ggplot(md, aes(nCount_RNA, nFeature_RNA)) +
        geom_bin2d(aes(fill = after_stat(log10(count))), bins = 40) +
        scale_x_log10() + scale_y_log10() + scale_fill_viridis_c(name = "log10 nuclei") +
        facet_wrap(~ Sample_f, ncol = 8) + th + theme(text = element_text(size = 7)))
print(ggplot(md, aes(nCount_RNA, percent.mt)) +
        geom_bin2d(aes(fill = after_stat(log10(count))), bins = 40) +
        scale_x_log10() + scale_fill_viridis_c(name = "log10 nuclei") +
        facet_wrap(~ Sample_f, ncol = 8) + th + theme(text = element_text(size = 7)))

# MT histograms
print(ggplot(md, aes(percent.mt, fill = Dx)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = c(5, 10, 20), linetype = "dashed", color = "grey40", linewidth = 0.3) +
        mt_x + dx_fill +
        facet_wrap(~ Sample_f, ncol = 8, scales = "free_y") +
        labs(title = "percent.mt per sample (dashed = 5, 10, 20%)", y = "Nuclei") +
        th + theme(text = element_text(size = 7)))

# Sample-level ranked dot plots (lines = cohort median ± 3 MAD, reference only)
rank_metrics <- c(med_nCount = TRUE, med_nFeature = TRUE, n_nuclei = TRUE,
                  med_percent.mt = FALSE, med_percent.rb = FALSE,
                  reads_per_UMI = FALSE, complexity_resid = FALSE)
rank_plot <- function(m, logscale) {
  v <- sample_tbl[[m]]
  tv <- if (logscale) log10(v) else v
  ref <- median(tv) + c(-3, 0, 3) * mad(tv)
  if (logscale) ref <- 10^ref
  d <- sample_tbl %>% mutate(val = v) %>% arrange(val) %>%
    mutate(Sample_ID = factor(Sample_ID, Sample_ID))
  p <- ggplot(d, aes(Sample_ID, val, color = Dx)) + geom_point(size = 2) +
    geom_hline(yintercept = ref, linetype = c("dashed", "solid", "dashed"), color = "grey50") +
    dx_col + labs(title = m, y = m, x = NULL) + th90
  if (logscale) p <- p + scale_y_log10()
  p
}
for (m in names(rank_metrics)) print(rank_plot(m, rank_metrics[[m]]))

# Complexity view
print(ggplot(sample_tbl, aes(med_nCount, med_nFeature, color = Dx)) +
        geom_line(aes(y = 10^fitted(fit)), color = "grey40") +
        geom_point(size = 2) +
        geom_text(aes(label = Sample_ID), size = 2.5, vjust = -0.8,
                  check_overlap = TRUE, show.legend = FALSE) +
        scale_x_log10() + scale_y_log10() + dx_col + th)

# Nuclei count vs depth
print(ggplot(sample_tbl, aes(med_nCount, n_nuclei, color = Dx)) +
        geom_point(size = 2) +
        geom_text(aes(label = Sample_ID), size = 2.5, vjust = -0.8,
                  check_overlap = TRUE, show.legend = FALSE) +
        scale_x_log10() + scale_y_log10() + dx_col + th +
        labs(subtitle = paste("Spearman rho =",
                              round(cor(sample_tbl$n_nuclei, sample_tbl$med_nCount,
                                        method = "spearman"), 2))))

# QC vs covariates
sl <- sample_tbl %>%
  pivot_longer(all_of(sm_metrics), names_to = "metric", values_to = "value")
print(sl %>%
        pivot_longer(c(Age, PMI, RIN), names_to = "covariate", values_to = "cov_value") %>%
        ggplot(aes(cov_value, value, color = Dx)) + geom_point(size = 1.2) +
        facet_grid(metric ~ covariate, scales = "free") + dx_col + th +
        theme(text = element_text(size = 8)))
print(ggplot(sl, aes(Dx, value, fill = Dx)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, size = 1) +
        facet_wrap(~ metric, scales = "free_y") + dx_fill + th)

dev.off()
