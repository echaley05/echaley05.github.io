# =============================================================================
# STEP 3  Distribution check: which transform, and what do the analytes look like?
#
# Input : output/Luminex_normalized.xlsx, output/Luminex_QC.xlsx (Detection)
# Output: output/Distribution_Check.xlsx   normality under each transform
#         output/Distribution_QQ.pdf       QQ plots, raw vs log2 vs ln(x+1) vs sqrt
#         output/Analyte_Boxplots.pdf      box plots per analyte (by group, by plate)
#
# Transforms: none, log2, ln(x+1), sqrt and a fitted Box-Cox. (log10, log2 and ln
# differ only by a constant multiplier, so they give identical skewness,
# Shapiro-Wilk W and model p-values; only log2 is shown.)
# Normality is judged three ways, because they answer different questions:
#   1. pooled values per analyte
#   2. within each group x region cell
#   3. residuals of the step-4 model, fitted per region: what actually has to be normal
# =============================================================================

if (!file.exists("R/pipeline_helpers.R"))
  stop("Run this from the Luminex_Pipeline folder: open Luminex_Pipeline.Rproj in RStudio, ",
       "or setwd(\"path/to/Luminex_Pipeline\") first.", call. = FALSE)
source("R/pipeline_helpers.R")
need_packages("MASS", attach = FALSE)   # used as MASS::boxcox (attaching it would mask dplyr::select)

DIST_XLSX <- out_path("Distribution_Check.xlsx")
QQ_PDF    <- out_path("Distribution_QQ.pdf")
BOX_PDF   <- out_path("Analyte_Boxplots.pdf")
BOX_SCALE <- "log2"       # y axis for the box plots: "log2", "log10" or "identity"

long <- load_long(complete_covariates = FALSE)
lev  <- levels(long$Group); cols <- group_colors(lev)
message(n_distinct(long$Analyte), " analytes assessed, ", nrow(long), " values")

TRANSFORMS <- list(none = function(x) x, log2 = log2, `ln(x+1)` = log1p, sqrt = sqrt)

skewness <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  if (n < 3 || sd(x) == 0) return(NA_real_); sum((x - mean(x))^3) / (n * sd(x)^3) }
kurtosis <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  if (n < 4 || sd(x) == 0) return(NA_real_); sum((x - mean(x))^4) / (n * sd(x)^4) - 3 }
shapiro <- function(x, what = "W") {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(x) > 5000 || sd(x) == 0) return(NA_real_)
  s <- try(shapiro.test(x), silent = TRUE)
  if (inherits(s, "try-error")) return(NA_real_)
  if (what == "W") unname(s$statistic) else s$p.value
}
boxcox_lambda <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 10) return(NA_real_)
  bc <- try(MASS::boxcox(v ~ 1, data = data.frame(v = x), lambda = seq(-2, 2, 0.05), plotit = FALSE), silent = TRUE)
  if (inherits(bc, "try-error")) NA_real_ else bc$x[which.max(bc$y)]
}
bc_apply <- function(x, l) if (is.na(l)) rep(NA_real_, length(x)) else if (abs(l) < 1e-8) log(x) else (x^l - 1) / l

# ---- 1. Pooled values per analyte ------------------------------------------------
pooled <- long |> group_split(Analyte) |> map_dfr(function(d) {
  lam <- boxcox_lambda(d$Value)
  tr <- c(TRANSFORMS, setNames(list(function(x) bc_apply(x, lam)), paste0("Box-Cox (lambda=", round(lam, 2), ")")))
  imap_dfr(tr, function(f, nm) { v <- f(d$Value)
    tibble(Transform = nm, N = sum(is.finite(v)), Skewness = skewness(v), Kurtosis = kurtosis(v),
           Shapiro_W = shapiro(v), Shapiro_p = shapiro(v, "p")) }) |>
    mutate(Analyte = d$Analyte[1], BoxCox_Lambda = lam, .before = 1)
})
best_pooled <- pooled |> group_by(Analyte) |> slice_max(Shapiro_W, n = 1, with_ties = FALSE) |>
  transmute(Analyte, Best_Transform = Transform, Best_W = Shapiro_W, BoxCox_Lambda = round(BoxCox_Lambda, 2)) |>
  ungroup()

# ---- 2. Within group x region cells ------------------------------------------------
by_cell <- long |> group_by(Analyte, Group, Region) |> filter(n() >= 8) |> group_split() |>
  map_dfr(function(d) imap_dfr(TRANSFORMS, function(f, nm) { v <- f(d$Value)
    tibble(Analyte = d$Analyte[1], Group = as.character(d$Group[1]), Region = as.character(d$Region[1]),
           Transform = nm, N = nrow(d), Shapiro_W = shapiro(v), Shapiro_p = shapiro(v, "p")) }))
cell_summary <- by_cell |> group_by(Transform) |>
  summarise(Cells = n(), Median_W = median(Shapiro_W, na.rm = TRUE),
            Pct_Normal_p05 = 100 * mean(Shapiro_p > 0.05, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(Median_W))

# ---- 3. Residuals of the step-4 model (Model A, per region) ----------------------------
resid_check <- long |> filter(if_all(all_of(COVARIATES), ~ !is.na(.x))) |>
  group_split(Analyte) |> map_dfr(function(d) imap_dfr(TRANSFORMS, function(f, nm) {
    r <- map(split(d, d$Region, drop = TRUE), function(dr) {
      dr <- dr |> mutate(yt = f(Value)) |> filter(is.finite(yt)) |> droplevels()
      if (nrow(dr) < MIN_N || nlevels(dr$Group) < 2) return(NULL)
      covs <- usable_terms(dr)
      m <- try(lm(as.formula(paste("yt ~ Group", if (length(covs)) paste("+", paste(covs, collapse = " + ")))),
                  data = dr), silent = TRUE)
      if (inherits(m, "try-error")) NULL else m
    }) |> compact()
    if (!length(r)) return(NULL)
    res <- unlist(map(r, residuals)); fit <- unlist(map(r, fitted))
    tibble(Analyte = d$Analyte[1], Transform = nm, N = length(res), Resid_Skew = skewness(res),
           Shapiro_W = shapiro(res), Shapiro_p = shapiro(res, "p"),
           Het_p = suppressWarnings(summary(lm(res^2 ~ fit))$coefficients[2, 4]))   # variance vs fitted
  }))
resid_summary <- resid_check |> group_by(Transform) |>
  summarise(Analytes = n(), Median_W = median(Shapiro_W, na.rm = TRUE),
            Pct_Resid_Normal = 100 * mean(Shapiro_p > 0.05, na.rm = TRUE),
            Pct_Homoscedastic = 100 * mean(Het_p > 0.05, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(Median_W))
if (nrow(resid_summary))
  message("Best transform by residual normality: ", resid_summary$Transform[1],
          "  (steps 4-6 use log2)")

write_workbook(list(Summary_Residuals = resid_summary, Summary_Cells = cell_summary,
                    Best_Per_Analyte = best_pooled, Pooled_All = pooled, Cells_All = by_cell,
                    Residuals_All = resid_check), DIST_XLSX, width = 15)

# ---- QQ plots ---------------------------------------------------------------------------
an <- sort(unique(long$Analyte))
qq <- long |> select(Analyte, Value) |>
  mutate(none = Value, log2 = log2(Value), `ln(x+1)` = log1p(Value), sqrt = sqrt(Value)) |>
  select(-Value) |> pivot_longer(-Analyte, names_to = "Transform", values_to = "v") |>
  mutate(Transform = factor(Transform, levels = names(TRANSFORMS)))
pdf(QQ_PDF, width = 11, height = 8.5)
for (ch in chunks(an, 6)) suppressWarnings(print(
  ggplot(qq |> filter(Analyte %in% ch), aes(sample = v)) +
    stat_qq(size = .5) + stat_qq_line(colour = "red") +
    facet_grid(Analyte ~ Transform, scales = "free") +
    labs(title = "Normal QQ plots by transform (pooled values)", x = "Theoretical", y = "Observed") +
    theme_bw(base_size = 8)))
invisible(dev.off())

# ---- Box plots -----------------------------------------------------------------------------
sc <- switch(BOX_SCALE, log2 = scale_y_continuous(trans = "log2"), log10 = scale_y_log10(), scale_y_continuous())
lp <- long |> mutate(Group = factor(as.character(Group), levels = plot_order(lev)))
pdf(BOX_PDF, width = 11, height = 8.5)
suppressWarnings(print(ggplot(long, aes(reorder(Analyte, Value, median), Value)) +
  geom_boxplot(outlier.size = .5, fill = "grey90") + sc + coord_flip() +
  labs(title = "All analytes (bridges excluded)", x = NULL, y = paste0("Normalized value (", BOX_SCALE, " scale)")) +
  theme_bw(base_size = 8)))
for (ch in chunks(an, 6)) suppressWarnings(print(
  ggplot(lp |> filter(Analyte %in% ch), aes(Group, Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = .8) + geom_jitter(width = .18, size = .4, alpha = .5) +
    scale_fill_manual(values = cols) + facet_grid(Analyte ~ Region, scales = "free_y") + sc +
    labs(title = "By group and region", x = NULL, y = "Normalized value") +
    theme_bw(base_size = 8) + theme(legend.position = "none")))
if (nlevels(long$Plate) > 1) for (ch in chunks(an, 6)) suppressWarnings(print(
  ggplot(long |> filter(Analyte %in% ch), aes(Plate, Value)) +
    geom_boxplot(outlier.size = .4, fill = "grey90") + facet_wrap(~ Analyte, scales = "free_y", ncol = 3) + sc +
    labs(title = "By plate", x = "Plate", y = "Normalized value") + theme_bw(base_size = 8)))
invisible(dev.off())

message("Wrote ", DIST_XLSX, ", ", QQ_PDF, ", ", BOX_PDF)
print(resid_summary)
