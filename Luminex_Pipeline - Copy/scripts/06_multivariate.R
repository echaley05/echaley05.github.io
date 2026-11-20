# =============================================================================
# STEP 6 (optional)  Multivariate analyses: co-expression modules, PERMANOVA,
# PLS-DA, and within-case region contrasts.
#
# Input : output/Luminex_normalized.xlsx, output/Luminex_QC.xlsx (Detection, Std_Curves)
# Output: output/Luminex_Multivariate.xlsx + output/Luminex_Multivariate_plots.pdf
#
# How covariates are handled, by analysis:
#   Modules      clustering runs on values with the covariates, Plate and Region
#                removed, so modules reflect co-regulation, not batch or anatomy.
#                Group is NEVER removed and never used to build modules.
#   Eigengenes   tested per region with the step-4 model (lmPerm, Model A).
#   PERMANOVA    covariates enter before Group (by = "margin"), within each region,
#                so each case appears once.
#   PLS-DA       run on the covariate-adjusted matrix; judged by cross-validation
#                and label permutation, not by how the score plot looks.
#   Region       within-case differences (e.g. ERC - STG), so case-level covariates
#   contrasts    cancel. Only when the study has 2+ regions.
#
# Imputation (THIS STEP ONLY, IMPUTE_MV): "OOR <" -> LLOQ/2 for that plate and analyte,
# converted with that sample's protein; "OOR >" -> ULOQ. "***"/"---" stay missing and
# samples with any gap are dropped from the matrix. Region contrasts use observed values only.
# =============================================================================

if (!file.exists("R/pipeline_helpers.R"))
  stop("Run this from the Luminex_Pipeline folder: open Luminex_Pipeline.Rproj in RStudio, ",
       "or setwd(\"path/to/Luminex_Pipeline\") first.", call. = FALSE)
source("R/pipeline_helpers.R")
need_packages(c("vegan", "pls", "cluster"), attach = FALSE)
need_file(FILES$norm, "01_normalize.R"); need_file(FILES$qc, "02_qc.R")
load_lmperm()
set.seed(SEED)

# ---- Load, with below/above-range codes kept for imputation ---------------------------------
long <- load_long(min_detect = min(MIN_DETECT, MV_DETECT), observed_only = FALSE)
limits <- read_excel(FILES$qc, sheet = "Std_Curves") |>
  transmute(Plate = factor(as.integer(Plate)), Analyte, LLOQ = as.numeric(LLOQ), ULOQ = as.numeric(ULOQ))
detect <- attr(long, "detect")
multi_region <- nlevels(long$Region) > 1
lev <- levels(long$Group); cols <- group_colors(lev); porder <- plot_order(lev)

long <- long |> left_join(limits, by = c("Plate", "Analyte")) |>
  mutate(lo_sub = if (IMPUTE_LOW == "lloq_sqrt2") LLOQ / sqrt(2) else LLOQ / 2,
         Value_imp = case_when(
           !is.na(Value) ~ Value,
           !IMPUTE_MV ~ NA_real_,
           coalesce(Code, "") == "OOR <" & !is.na(LLOQ) ~ lo_sub / Protein_pgmL * SCALE,
           coalesce(Code, "") == "OOR >" & IMPUTE_HIGH & !is.na(ULOQ) ~ ULOQ / Protein_pgmL * SCALE,
           TRUE ~ NA_real_),
         Imputed = is.na(Value) & !is.na(Value_imp),
         y_imp = log2(Value_imp))

keep_mv <- detect |> filter(Pct_Detected >= MV_DETECT) |> pull(Analyte)
keep_rc <- detect |> filter(Pct_Detected >= MIN_DETECT) |> pull(Analyte)
meta <- long |> select(Sample_ID, Case_ID, Region, Plate, Group, all_of(COVARIATES)) |> distinct()

wide <- long |> filter(Analyte %in% keep_mv) |> select(Sample_ID, Analyte, y_imp) |>
  pivot_wider(names_from = Analyte, values_from = y_imp)
complete <- wide |> filter(if_all(-Sample_ID, ~ !is.na(.x))) |> arrange(Sample_ID)
X <- as.matrix(complete |> select(-Sample_ID)); rownames(X) <- complete$Sample_ID
mv_meta <- meta |> filter(Sample_ID %in% rownames(X)) |> arrange(match(Sample_ID, rownames(X)))
stopifnot(identical(rownames(X), mv_meta$Sample_ID))
message(ncol(X), " analytes x ", nrow(X), " samples in the multivariate matrix (",
        n_distinct(mv_meta$Case_ID), " cases); ", nrow(wide) - nrow(X), " samples dropped for gaps")
if (ncol(X) < 3) stop("Fewer than 3 analytes reach MV_DETECT = ", MV_DETECT, "%. Lower MV_DETECT.", call. = FALSE)

# ---- Covariate adjustment (Group deliberately excluded) -------------------------------------------
adjust <- function(X, meta, drop_region) {
  terms <- c(usable_terms(meta), if (drop_region && multi_region) "Region")
  f <- as.formula(paste("v ~", if (length(terms)) paste(terms, collapse = " + ") else "1"))
  out <- apply(X, 2, function(v) as.numeric(residuals(lm(f, data = mutate(meta, v = v)))))
  rownames(out) <- rownames(X); out
}
X_adj_noregion <- adjust(X, mv_meta, drop_region = TRUE)    # for clustering
X_adj          <- adjust(X, mv_meta, drop_region = FALSE)   # for PLS-DA

# =============================================================================
# 1. MODULES AND EIGENGENES
# =============================================================================
cor_mat  <- cor(X_adj_noregion, method = "spearman")
dist_mat <- as.dist(1 - abs(cor_mat))
hc <- hclust(dist_mat, method = "average")
ks <- 2:max(2, min(6, ncol(X) - 1))
sil <- map_dbl(ks, function(k) mean(cluster::silhouette(cutree(hc, k), dist_mat)[, 3]))
k_best <- ks[which.max(sil)]
module_tbl <- tibble(Analyte = names(cutree(hc, k_best)), Module = paste0("M", cutree(hc, k_best))) |>
  arrange(Module, Analyte)
message(k_best, " modules (silhouette ", round(max(sil), 3), ")")

# eigengene = PC1 of the module's standardized adjusted values, signed to match the module mean
eigen_tbl <- map_dfr(unique(module_tbl$Module), function(m) {
  a <- module_tbl$Analyte[module_tbl$Module == m]
  sub <- scale(X_adj[, a, drop = FALSE])
  p <- prcomp(sub, center = TRUE, scale. = FALSE)
  pc1 <- p$x[, 1]; if (length(a) > 1 && cor(pc1, rowMeans(sub)) < 0) pc1 <- -pc1
  tibble(Sample_ID = rownames(X), Module = m, Eigengene = as.numeric(pc1),
         Var_Explained = summary(p)$importance[2, 1], N_Analytes = length(a))
})
eig <- eigen_tbl |> left_join(mv_meta, by = "Sample_ID")

eig_fits <- eig |> group_split(Module, Region) |> map(function(d) {
  set.seed(fit_seed("eigengene", d$Module[1], as.character(d$Region[1])))
  r <- fit_group_model(droplevels(mutate(d, y = Eigengene)), interaction = FALSE)
  if (is.null(r)) return(NULL)
  base <- tibble(Module = d$Module[1], Region = as.character(d$Region[1]))
  list(terms = bind_cols(base, r$terms), pairs = bind_cols(base, r$pairs |> rename(Difference = Log2FC) |>
                                                              select(-Fold_Change)))
}) |> compact()
eigen_terms <- map_dfr(eig_fits, "terms") |>
  group_by(Region, Term) |> mutate(FDR = p.adjust(Perm_P, "BH")) |> ungroup() |> arrange(Term, Region, Perm_P)
eigen_contrasts <- map_dfr(eig_fits, "pairs") |>
  group_by(Region, Contrast) |> mutate(FDR = p.adjust(Perm_P, "BH")) |> ungroup() |>
  mutate(Signif = FDR < FDR_ALPHA) |> arrange(FDR)

# =============================================================================
# 2. PERMANOVA (within region: one sample per case, so rows are independent)
# =============================================================================
by_region <- split(seq_len(nrow(mv_meta)), mv_meta$Region, drop = TRUE)
permanova <- imap_dfr(by_region, function(idx, rg) {
  d <- droplevels(mv_meta[idx, ]); if (length(idx) < 15 || nlevels(d$Group) < 2) return(NULL)
  dm <- vegan::vegdist(X[idx, , drop = FALSE], method = "euclidean")
  set.seed(fit_seed("permanova", rg))
  f <- as.formula(paste("dm ~", paste(c(usable_terms(d), "Group"), collapse = " + ")))
  a <- vegan::adonis2(f, data = d, permutations = MV_N_PERM, by = "margin")
  bp <- vegan::permutest(vegan::betadisper(dm, d$Group), permutations = MV_N_PERM)
  as.data.frame(a) |> tibble::rownames_to_column("Term") |>
    transmute(Region = rg, Term, Df, SumOfSqs, R2, F = `F`, P = `Pr(>F)`) |>
    mutate(Group_Dispersion_P = bp$tab$`Pr(>F)`[1])
})
pcoa <- imap_dfr(by_region, function(idx, rg) {
  if (length(idx) < 15) return(NULL)
  p <- cmdscale(vegan::vegdist(X[idx, , drop = FALSE], method = "euclidean"), k = 2, eig = TRUE)
  tibble(Region = rg, Sample_ID = mv_meta$Sample_ID[idx], Group = as.character(mv_meta$Group[idx]),
         Axis1 = p$points[, 1], Axis2 = p$points[, 2])
})

# =============================================================================
# 3. PLS-DA (covariate-adjusted data, cross-validated and permutation-tested)
# =============================================================================
plsda_fit <- function(Xm, grp, ncomp) {
  Y <- model.matrix(~ 0 + grp); colnames(Y) <- levels(grp)
  pls::plsr(Y ~ X, ncomp = ncomp, data = data.frame(Y = I(Y), X = I(Xm)), scale = TRUE)
}
plsda_cv_acc <- function(Xm, grp, ncomp, folds = 5) {
  n <- length(grp); f <- sample(rep_len(1:folds, n)); correct <- rep(NA, n)
  for (k in 1:folds) {
    tr <- f != k; te <- !tr
    if (n_distinct(grp[tr]) < 2) next
    m <- plsda_fit(Xm[tr, , drop = FALSE], droplevels(grp[tr]), ncomp)
    pr <- predict(m, newdata = data.frame(X = I(Xm[te, , drop = FALSE])), ncomp = ncomp)
    pm <- matrix(pr[, , 1], nrow = sum(te))                  # samples x groups (kept 2-D for one test sample)
    correct[te] <- dimnames(pr)[[2]][max.col(pm, ties.method = "first")] == as.character(grp[te])
  }
  mean(correct, na.rm = TRUE)
}
vip_scores <- function(m, ncomp) {
  W <- m$loading.weights[, 1:ncomp, drop = FALSE]; Q <- m$Yloadings[, 1:ncomp, drop = FALSE]
  Tm <- m$scores[, 1:ncomp, drop = FALSE]; SS <- colSums(Q^2) * colSums(Tm^2)
  Wn <- sweep(W^2, 2, colSums(W^2), "/")
  sqrt(nrow(W) * (Wn %*% SS) / sum(SS))[, 1]
}
plsda_perf <- list(); plsda_scores <- list(); plsda_vip <- list()
for (rg in names(by_region)) {
  idx <- by_region[[rg]]; if (length(idx) < 20) next
  Xm <- X_adj[idx, , drop = FALSE]; grp <- droplevels(mv_meta$Group[idx])
  if (nlevels(grp) < 2) next
  set.seed(fit_seed("plsda", rg))
  ncomp <- min(MAX_COMP, ncol(Xm) - 1, length(grp) - 2)
  acc <- map_dbl(1:ncomp, ~ plsda_cv_acc(Xm, grp, .x)); best <- which.max(acc)
  n_show <- min(max(2, best), ncomp)
  m <- plsda_fit(Xm, grp, n_show)
  null_acc <- map_dbl(1:MV_N_PERM, function(i) plsda_cv_acc(Xm, sample(grp), best))
  p_perm <- (1 + sum(null_acc >= acc[best])) / (1 + MV_N_PERM)
  plsda_perf[[rg]] <- tibble(Region = rg, N = length(grp), N_Comp = best, CV_Accuracy = acc[best],
    Null_Mean_Accuracy = mean(null_acc), Chance_Accuracy = max(prop.table(table(grp))), Perm_P = p_perm,
    Verdict = if_else(p_perm < 0.05, "separates better than chance", "no separation beyond chance"))
  sc <- as.data.frame(unclass(m$scores)[, 1:n_show, drop = FALSE]); names(sc) <- paste0("Comp", seq_len(ncol(sc)))
  plsda_scores[[rg]] <- bind_cols(tibble(Region = rg, Group = as.character(grp)), sc)
  plsda_vip[[rg]] <- tibble(Region = rg, Analyte = colnames(Xm), VIP = vip_scores(m, best)) |> arrange(desc(VIP))
}
plsda_perf <- bind_rows(plsda_perf); plsda_scores <- bind_rows(plsda_scores); plsda_vip <- bind_rows(plsda_vip)

# =============================================================================
# 4. WITHIN-CASE REGION CONTRASTS (2+ regions only)
# =============================================================================
dunn_test <- function(x, g) {
  g <- droplevels(factor(g)); r <- rank(x); N <- length(x)
  ties <- table(x); sigma2 <- (N * (N + 1) / 12) - sum(ties^3 - ties) / (12 * (N - 1))
  map_dfr(combn(levels(g), 2, simplify = FALSE), function(pr) {
    i <- pr[1]; j <- pr[2]; ni <- sum(g == i); nj <- sum(g == j)
    z <- (mean(r[g == j]) - mean(r[g == i])) / sqrt(sigma2 * (1 / ni + 1 / nj))
    tibble(Contrast = paste(j, "vs", i), Z = z, P = 2 * pnorm(-abs(z)))
  })
}
region_contrasts <- tibble(); region_posthoc <- tibble(); rc_long <- NULL
if (multi_region) {
  rc_long <- long |> filter(Analyte %in% keep_rc, !is.na(y)) |>
    select(Case_ID, Group, Region, Analyte, y) |> pivot_wider(names_from = Region, values_from = y)
  rg_pairs <- combn(levels(long$Region), 2, simplify = FALSE)
  region_contrasts <- map_dfr(rg_pairs, function(pr) {
    a <- pr[2]; b <- pr[1]
    d <- rc_long |> mutate(Diff = .data[[a]] - .data[[b]]) |> filter(!is.na(Diff))
    d |> group_split(Analyte) |> map_dfr(function(dd) {
      ng <- table(droplevels(dd$Group))
      if (length(ng) < 2 || any(ng < MIN_GROUP)) return(NULL)
      kw <- kruskal.test(Diff ~ Group, data = dd)
      med <- dd |> group_by(Group) |> summarise(M = median(Diff), N = n(), .groups = "drop")
      tibble(Comparison = paste(a, "-", b), Analyte = dd$Analyte[1], N = nrow(dd),
             Overall_Median_Diff = median(dd$Diff), Overall_P = wilcox.test(dd$Diff)$p.value,
             KW_Chi = unname(kw$statistic), KW_P = kw$p.value,
             !!!set_names(as.list(round(med$M, 3)), paste0("Median_", med$Group)),
             !!!set_names(as.list(med$N), paste0("N_", med$Group)))
    })
  })
  if (nrow(region_contrasts)) {
    region_contrasts <- region_contrasts |> group_by(Comparison) |>
      mutate(Overall_FDR = p.adjust(Overall_P, "BH"), KW_FDR = p.adjust(KW_P, "BH")) |> ungroup() |> arrange(KW_FDR)
    sig_rc <- region_contrasts |> filter(KW_FDR < FDR_ALPHA)
    if (nrow(sig_rc)) region_posthoc <- map_dfr(seq_len(nrow(sig_rc)), function(i) {
      r <- sig_rc[i, ]; pr <- str_split(r$Comparison, " - ")[[1]]
      d <- rc_long |> filter(Analyte == r$Analyte) |> mutate(Diff = .data[[pr[1]]] - .data[[pr[2]]]) |> filter(!is.na(Diff))
      dunn_test(d$Diff, d$Group) |> mutate(Comparison = r$Comparison, Analyte = r$Analyte, .before = 1)
    }) |> mutate(FDR = p.adjust(P, "BH"))
  }
}

# =============================================================================
# WRITE
# =============================================================================
imp_tbl <- long |> group_by(Analyte) |>
  summarise(N = n(), Observed = sum(!is.na(Value)),
            Imputed_low = sum(Imputed & Code %in% "OOR <"), Imputed_high = sum(Imputed & Code %in% "OOR >"),
            Still_missing = sum(is.na(y_imp)), Pct_Imputed = round(100 * mean(Imputed), 1),
            In_Multivariate = Analyte[1] %in% keep_mv, .groups = "drop") |> arrange(desc(Pct_Imputed))
write_workbook(list(
  Modules = module_tbl |> left_join(eigen_tbl |> distinct(Module, Var_Explained, N_Analytes), by = "Module"),
  Eigengene_Terms = eigen_terms, Eigengene_Contrasts = eigen_contrasts,
  PERMANOVA = permanova |> mutate(Signif = Term == "Group" & P < FDR_ALPHA),
  PLSDA_Performance = plsda_perf, PLSDA_VIP = plsda_vip,
  Region_Contrasts = if (multi_region) region_contrasts |> mutate(Signif = KW_FDR < FDR_ALPHA),
  Region_PostHoc = if (multi_region) region_posthoc,
  Imputation = imp_tbl,
  Analyte_Correlations = as.data.frame(round(cor_mat, 3)) |> tibble::rownames_to_column("Analyte")),
  FILES$mv)

# =============================================================================
# PLOTS
# =============================================================================
pdf(FILES$mv_pdf, width = 11, height = 8.5)
plot(hc, main = paste0("Analyte clustering (", k_best, " modules)"), xlab = "", sub = "", cex = .7)
rect.hclust(hc, k = k_best, border = "red")
hm <- as.data.frame(cor_mat) |> tibble::rownames_to_column("A") |> pivot_longer(-A, names_to = "B", values_to = "rho") |>
  mutate(A = factor(A, levels = hc$labels[hc$order]), B = factor(B, levels = hc$labels[hc$order]))
print(ggplot(hm, aes(A, B, fill = rho)) + geom_tile() +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
        labs(title = "Analyte correlation (covariate-adjusted, Spearman)", x = NULL, y = NULL) +
        theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)))
print(ggplot(eig |> mutate(Group = factor(as.character(Group), levels = porder)), aes(Group, Eigengene, fill = Group)) +
        geom_boxplot(outlier.shape = NA, alpha = .85, width = .6) + geom_jitter(width = .16, size = .8, alpha = .5) +
        scale_fill_manual(values = cols) + facet_grid(Module ~ Region, scales = "free_y") +
        labs(title = "Module eigengenes by group and region", subtitle = "values are covariate-adjusted", x = NULL) +
        theme_bw(base_size = 10) + theme(legend.position = "none"))
if (nrow(pcoa)) print(ggplot(pcoa |> mutate(Group = factor(Group, levels = porder)), aes(Axis1, Axis2, colour = Group)) +
        geom_point(size = 2, alpha = .85) + stat_ellipse(level = .68) + scale_colour_manual(values = cols) +
        facet_wrap(~ Region, scales = "free") +
        labs(title = "PCoA of analyte profiles (Euclidean, log2 values)",
             subtitle = "ellipses are 68% contours; see the PERMANOVA sheet for tests") + theme_bw(base_size = 10))
if (nrow(plsda_scores)) {
  print(ggplot(plsda_scores |> mutate(Group = factor(Group, levels = porder)), aes(Comp1, Comp2, colour = Group)) +
          geom_point(size = 2, alpha = .85) + stat_ellipse(level = .68) + scale_colour_manual(values = cols) +
          facet_wrap(~ Region, scales = "free") +
          labs(title = "PLS-DA scores (covariate-adjusted data)",
               subtitle = "supervised projection: read only alongside the permutation p") + theme_bw(base_size = 10))
  print(ggplot(plsda_vip |> group_by(Region) |> slice_max(VIP, n = 12) |> ungroup(), aes(reorder(Analyte, VIP), VIP)) +
          geom_col(fill = "grey40") + coord_flip() + geom_hline(yintercept = 1, linetype = 2, colour = "red") +
          facet_wrap(~ Region, scales = "free_y") +
          labs(title = "PLS-DA variable importance (VIP > 1 conventionally notable)", x = NULL) + theme_bw(base_size = 9))
}
if (multi_region && nrow(region_contrasts)) {
  top_rc <- region_contrasts |> slice_head(n = 8)
  d <- map_dfr(seq_len(nrow(top_rc)), function(i) {
    r <- top_rc[i, ]; pr <- str_split(r$Comparison, " - ")[[1]]
    rc_long |> filter(Analyte == r$Analyte) |>
      transmute(Analyte, Group = as.character(Group), Comparison = r$Comparison, Diff = .data[[pr[1]]] - .data[[pr[2]]]) |>
      filter(!is.na(Diff))
  })
  print(ggplot(d |> mutate(Group = factor(Group, levels = porder)), aes(Group, Diff, fill = Group)) +
          geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
          geom_boxplot(outlier.shape = NA, alpha = .85, width = .6) + geom_jitter(width = .16, size = .7, alpha = .5) +
          scale_fill_manual(values = cols) + facet_wrap(Analyte ~ Comparison, scales = "free_y", ncol = 4) +
          labs(title = "Within-case regional differences (log2)",
               subtitle = "each point is one case; 0 = no difference between regions", x = NULL, y = "log2 difference") +
          theme_bw(base_size = 9) + theme(legend.position = "none"))
}
invisible(dev.off())
message("\nWrote ", FILES$mv, " and ", FILES$mv_pdf)
print(plsda_perf)
print(permanova |> filter(Term == "Group") |> select(Region, R2, P, Group_Dispersion_P))
