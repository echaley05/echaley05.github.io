############################################
# 03_filter_preview.R — preview the final filters
# Nothing is saved as a filtered object
############################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

out <- "03_filter_preview"
dir.create(out, showWarnings = FALSE)

seu <- readRDS("02_qc_metrics.rds")

############################################
# Filter settings
############################################

exclude_samples <- c("77_4KGA_ERC", "57_RMGG_ERC", "95_K47B_ERC")
nFeature_min <- 500
nFeature_max <- 15000
nCount_min   <- 1000
nCount_max   <- 200000
mt_max       <- 10
rb_max       <- 1

############################################
# 1. Flag nuclei
############################################

md <- seu@meta.data %>%
  mutate(Dx = factor(Dx, levels = c("Control", "ASD")),
         excluded_sample = Sample_ID %in% exclude_samples)

md$f_nFeature_low  <- md$nFeature_RNA < nFeature_min
md$f_nCount_low    <- md$nCount_RNA   < nCount_min
md$f_percent.mt    <- md$percent.mt   > mt_max
md$f_nFeature_high <- md$nFeature_RNA > nFeature_max
md$f_nCount_high   <- md$nCount_RNA   > nCount_max
md$f_percent.rb    <- md$percent.rb  >= rb_max

f_cols <- grep("^f_", names(md), value = TRUE)
md$retained <- !md$excluded_sample & !Reduce(`|`, md[f_cols])

# Cascade: each nucleus is counted under the first filter it fails
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
reasons <- c("Excluded sample",
             paste0("nFeature < ", fmt(nFeature_min)),
             paste0("nCount < ", fmt(nCount_min)),
             paste0("percent.mt > ", mt_max),
             paste0("nFeature > ", fmt(nFeature_max)),
             paste0("nCount > ", fmt(nCount_max)),
             paste0("percent.rb >= ", rb_max),
             "Retained")

md$reason <- factor(case_when(
  md$excluded_sample ~ reasons[1],
  md$f_nFeature_low  ~ reasons[2],
  md$f_nCount_low    ~ reasons[3],
  md$f_percent.mt    ~ reasons[4],
  md$f_nFeature_high ~ reasons[5],
  md$f_nCount_high   ~ reasons[6],
  md$f_percent.rb    ~ reasons[7],
  TRUE               ~ reasons[8]), levels = reasons)

sample_levels <- md %>% group_by(Sample_ID, Dx) %>%
  summarise(d = median(nCount_RNA), .groups = "drop") %>%
  arrange(Dx, d) %>% pull(Sample_ID)
md$Sample_f <- factor(md$Sample_ID, levels = sample_levels)

############################################
# 2. Tables
############################################

by_sample <- md %>% dplyr::count(Sample_ID, Dx, reason) %>%
  pivot_wider(names_from = reason, values_from = n, values_fill = 0)
for (r in setdiff(reasons, names(by_sample))) by_sample[[r]] <- 0
by_sample <- by_sample %>%
  mutate(total = rowSums(across(all_of(reasons))),
         pct_retained = round(100 * .data[["Retained"]] / total, 1))

# Each filter on its own (nuclei can fail more than one)
independent <- md %>% filter(!excluded_sample) %>% group_by(Sample_ID, Dx) %>%
  summarise(n_nuclei = n(),
            across(all_of(f_cols), ~ round(100 * mean(.x), 1)),
            pct_retained = round(100 * mean(retained), 1), .groups = "drop")

dx_summary <- md %>% group_by(Dx) %>% summarise(
  samples_total = n_distinct(Sample_ID),
  samples_kept = n_distinct(Sample_ID[!excluded_sample]),
  nuclei_all = n(),
  nuclei_in_kept_samples = sum(!excluded_sample),
  nuclei_retained = sum(retained),
  pct_retained_of_kept_samples = round(100 * sum(retained) / sum(!excluded_sample), 1))

write.csv(by_sample, file.path(out, "removal_by_filter_cascade.csv"), row.names = FALSE)
write.csv(independent, file.path(out, "removal_by_filter_independent.csv"), row.names = FALSE)
write.csv(dx_summary, file.path(out, "summary_by_Dx.csv"), row.names = FALSE)

print(as.data.frame(dx_summary))
cat("\nLowest retention among kept samples:\n")
print(independent %>% arrange(pct_retained) %>% select(Sample_ID, Dx, n_nuclei, pct_retained) %>% head(12))

############################################
# 3. Figures
############################################

dx_cols <- c(Control = "#4DBBD5", ASD = "#E64B35")
th <- theme_bw(base_size = 10)
th90 <- th + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf(file.path(out, "filter_preview.pdf"), width = 16, height = 9)

# Share of nuclei removed by each filter (cascade order = legend order)
print(ggplot(md, aes(Sample_f, fill = reason)) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c("grey35", "#D95F02", "#E6AB02", "#1B9E77",
                                     "#7570B3", "#E7298A", "#66A61E", "grey85")) +
        facet_grid(~ Dx, scales = "free_x", space = "free_x") +
        labs(title = "Nuclei removed by each filter (counted under the first filter failed)",
             x = NULL, y = "Fraction of nuclei", fill = NULL) + th90)

# Nuclei before and after
counts_df <- by_sample %>%
  transmute(Sample_ID, Dx, Before = total, After = Retained) %>%
  pivot_longer(c(Before, After), names_to = "stage", values_to = "n") %>%
  mutate(stage = factor(stage, levels = c("Before", "After")),
         Sample_ID = factor(Sample_ID, sample_levels))
print(ggplot(counts_df, aes(Sample_ID, n, fill = stage)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c(Before = "grey75", After = "#2C3E50")) +
        facet_grid(~ Dx, scales = "free_x", space = "free_x") +
        labs(title = "Nuclei per sample before and after filtering", x = NULL, y = "Nuclei", fill = NULL) +
        th90)

# Each filter on its own, kept samples only
ind_long <- independent %>%
  pivot_longer(all_of(f_cols), names_to = "filter", values_to = "pct") %>%
  mutate(Sample_ID = factor(Sample_ID, sample_levels))
print(ggplot(ind_long, aes(filter, Sample_ID, fill = pct)) +
        geom_tile() + geom_text(aes(label = round(pct)), size = 2.5) +
        scale_fill_gradient(low = "white", high = "#B2182B", name = "% of nuclei") +
        labs(title = "% of nuclei failing each filter on its own", x = NULL, y = NULL) + th90)

dev.off()
