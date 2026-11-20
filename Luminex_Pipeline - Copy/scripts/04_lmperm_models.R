# =============================================================================
# STEP 4  Main analysis: permutation linear models per region (lmPerm)
#
# Input : output/Luminex_normalized.xlsx, output/Luminex_QC.xlsx (Detection)
# Output: output/Luminex_lmPerm_Results.xlsx
#
# Fitted separately in each region (one row per case) for every analyte detected
# in >= MIN_DETECT % of samples:
#
#   Model A:  log2(value) ~ Group + COVARIATES [+ Plate]
#   Model B:  Model A + INTERACTION_COVARIATE:Group       (skipped if NULL)
#
#   Group      factor with any number of levels; the first in GROUP_LEVELS is the reference
#   Covariates numeric ones are centred within the region, so Model B's group
#              comparisons are differences at the region's mean (e.g. mean age).
#              A covariate with one value in a region (all male) is dropped there.
#   Plate      fixed factor, when ADJUST_FOR_PLATE and the region spans > 1 plate
#
# Permutation: lmPerm::lmp(), perm = "Prob", maxIter = N_PERM, with early stopping
# disabled so every term runs the full N_PERM iterations. lmPerm's one built-in
# exception: a term whose first 50 iterations ALL beat the observed statistic stops
# at 51 with p = 1 (clearly non-significant; see the Iter column).
#
# Per analyte x region x model:
#   1. factor coding    -> term p-values (Group, each covariate, Plate, [Age:Group])
#   2. indicator coding -> every pairwise group comparison as its own 1-df term
#      (one fit per reference group), plus slope comparisons in Model B.
# Groups with fewer than MIN_GROUP cases in a region are left out of that region.
#
# Multiple testing: Benjamini-Hochberg across analytes, separately within each
#   region x term and each region x comparison.
#
# lmPerm 2.1.0 is required: 2.1.6 (CRAN, Nov 2025) returns badly miscalibrated
# p-values when a factor has 3+ levels. A calibration self-test on your own design
# runs before the analysis and stops if p-values are not trustworthy.
# =============================================================================

if (!file.exists("R/pipeline_helpers.R"))
  stop("Run this from the Luminex_Pipeline folder: open Luminex_Pipeline.Rproj in RStudio, ",
       "or setwd(\"path/to/Luminex_Pipeline\") first.", call. = FALSE)
source("R/pipeline_helpers.R")
need_file(FILES$norm, "01_normalize.R"); need_file(FILES$qc, "02_qc.R")
load_lmperm()
message("lmPerm version: ", packageVersion("lmPerm"))

MODELS <- list(A = list(interaction = FALSE))
if (!is.null(INTERACTION_COVARIATE)) MODELS$B <- list(interaction = TRUE)
model_label <- function(m) {     # terms usable somewhere in the data; per-fit terms are in the Covariates column
  rhs <- paste(c("Group", LABEL_TERMS), collapse = " + ")
  if (MODELS[[m]]$interaction) rhs <- paste0(rhs, " + ", INTERACTION_COVARIATE, ":Group")
  paste("log2(value) ~", rhs)
}

# ---- Load -----------------------------------------------------------------------------
long <- load_long()
detect <- attr(long, "detect")
if (!is.null(INTERACTION_COVARIATE) && !is.numeric(long[[INTERACTION_COVARIATE]]))
  stop("INTERACTION_COVARIATE '", INTERACTION_COVARIATE, "' must be numeric (e.g. Age).", call. = FALSE)
lev <- levels(long$Group)
LABEL_TERMS <- unique(unlist(map(split(long, long$Region, drop = TRUE), usable_terms)))
LABEL_TERMS <- c(intersect(COVARIATES, LABEL_TERMS), intersect("Plate", LABEL_TERMS))
message(n_distinct(long$Analyte), " analytes (>= ", MIN_DETECT, "% detected), ",
        n_distinct(long$Case_ID), " cases, ", nrow(long), " values")
message("Groups (first = reference): ", paste(lev, collapse = ", "))
print(as.data.frame(design_table(long)))
if (length(lev) < 2) stop("Only one group in the data; nothing to compare.", call. = FALSE)

# ---- Calibration self-test ------------------------------------------------------------
# Pure-noise outcomes on your real design (covariates, groups, plates): no term has a
# real effect, so a working permutation test gives p < 0.05 about 5% of the time.
# Version 2.1.6 gave 30-85% in testing.
if (RUN_SELFTEST) {
  set.seed(SEED)
  rg0 <- names(which.max(table(distinct(long, Region, Case_ID)$Region)))   # the region with most cases
  d0 <- long |> filter(Region == rg0) |> distinct(Case_ID, .keep_all = TRUE) |> droplevels()
  covs0 <- usable_terms(d0)
  for (v in covs0) if (is.numeric(d0[[v]])) d0[[v]] <- d0[[v]] - mean(d0[[v]])
  f0 <- as.formula(paste("y ~ Group", if (length(covs0)) paste("+", paste(covs0, collapse = " + "))))
  check_terms <- c("Group", head(covs0[vapply(covs0, function(v) is.numeric(d0[[v]]), logical(1))], 1))
  n_sim <- 100
  message("Self-test: ", n_sim, " null datasets in ", rg0, " (n = ", nrow(d0), ") ...")
  sim_p <- replicate(n_sim, {
    d0$y <- rnorm(nrow(d0))
    t <- lmp_terms(run_lmp(f0, d0))
    t$Perm_P[match(check_terms, t$Term)]
  })
  fp <- rowMeans(matrix(sim_p, nrow = length(check_terms)) < 0.05)
  selftest_txt <- paste0(paste(sprintf("%s %.2f", check_terms, fp), collapse = ", "),
                         " false-positive rate on ", n_sim, " null datasets (expected ~0.05)")
  message("Self-test: ", selftest_txt)
  if (any(fp > 0.12))
    stop("lmPerm failed the calibration self-test: p-values are far too small on data with ",
         "no real effects. Install lmPerm 2.1.0 (see README) and rerun.", call. = FALSE)
  selftest_txt <- paste("Passed:", selftest_txt)
} else selftest_txt <- "Not run"

# ---- Run ---------------------------------------------------------------------------------
# Each fit is seeded from SEED + model + region + analyte (see fit_seed), so the same
# data and settings give identical p-values on every run, on every computer with the
# same R and lmPerm versions.
jobs <- long |> group_by(Region, Analyte) |> group_split()
res <- list()
for (m in names(MODELS)) {
  message("Model ", m, ": ", model_label(m))
  t0 <- Sys.time()
  res[[m]] <- map(jobs, function(d) {
    set.seed(fit_seed(m, as.character(d$Region[1]), d$Analyte[1]))
    r <- fit_group_model(droplevels(d), interaction = MODELS[[m]]$interaction)
    if (is.null(r)) return(NULL)
    base <- tibble(Model = m, Region = as.character(d$Region[1]), Analyte = d$Analyte[1])
    list(terms = bind_cols(base, r$terms), pairs = bind_cols(base, r$pairs),
         slopes = if (!is.null(r$slopes)) bind_cols(base, r$slopes))
  }, .progress = interactive()) |> compact()
  message(sprintf("  %d region x analyte fits in %.1f min", length(res[[m]]),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (!length(res$A)) stop("No region x analyte could be fitted. Check MIN_N / MIN_GROUP against the design table above.",
                         call. = FALSE)

det_lookup <- detect |> select(Analyte, Pct_Detected)
region_order <- function(x) x |> mutate(Region = factor(Region, levels = levels(long$Region)))
contrast_order <- group_pairs(lev)$Contrast

finish_terms <- function(x) x |>
  group_by(Region, Term) |> mutate(FDR = p.adjust(Perm_P, "BH")) |> ungroup() |>
  mutate(Signif = FDR < FDR_ALPHA) |>
  relocate(Df, SumSq, .after = Term) |> relocate(FDR, Signif, .after = Iter) |>
  region_order() |> arrange(Term, Region, Perm_P)

finish_pairs <- function(x) x |>
  group_by(Region, Contrast) |> mutate(FDR = p.adjust(Perm_P, "BH")) |> ungroup() |>
  mutate(Signif = FDR < FDR_ALPHA, Contrast = factor(Contrast, levels = contrast_order)) |>
  left_join(det_lookup, by = "Analyte") |>
  relocate(FDR, Signif, .after = Perm_P) |>
  region_order() |> arrange(Region, Contrast, Perm_P) |> mutate(Contrast = as.character(Contrast))

out <- list()
for (m in names(MODELS)) {
  out[[paste0("Terms_", m)]]    <- finish_terms(map_dfr(res[[m]], "terms"))
  out[[paste0("Pairwise_", m)]] <- finish_pairs(map_dfr(res[[m]], "pairs"))
  if (MODELS[[m]]$interaction)
    out[[paste0("Slopes_", m)]] <- finish_pairs(map_dfr(res[[m]], "slopes")) |> select(-Pct_Detected)
}

# Iteration check: every test should run N_PERM iterations, except lmPerm's early stop at 51 when p = 1.
all_iter <- bind_rows(lapply(out, \(x) select(x, Perm_P, Iter)))
odd <- all_iter |> filter(Iter != N_PERM, !(Iter == 51 & Perm_P == 1))
if (nrow(odd)) warning(nrow(odd), " tests ran an unexpected number of iterations; see the Iter columns.", call. = FALSE)
n_early <- sum(all_iter$Iter == 51 & all_iter$Perm_P == 1)

dropped <- map_dfr(res$A, "terms") |> filter(Groups_Dropped != "") |> distinct(Region, Groups_Dropped)
if (nrow(dropped)) message("Groups below MIN_GROUP left out: ",
                           paste(paste0(dropped$Groups_Dropped, " in ", dropped$Region), collapse = "; "))

# ---- Summary + settings --------------------------------------------------------------------
summ <- bind_rows(lapply(names(MODELS), function(m) out[[paste0("Pairwise_", m)]] |>
  group_by(Region, Contrast) |>
  summarise(Model = m, Tested = n(), P_lt_05 = sum(Perm_P < 0.05), FDR_lt_05 = sum(FDR < FDR_ALPHA),
            Significant = paste(Analyte[FDR < FDR_ALPHA], collapse = ", "), .groups = "drop"))) |>
  relocate(Model) |> mutate(Contrast = factor(Contrast, levels = contrast_order)) |>
  arrange(Model, Region, Contrast) |> mutate(Contrast = as.character(Contrast))

settings <- tibble(
  Setting = c("Model A", "Model B", "Groups (reference first)", "Group plot order", "Interaction covariate",
              "Fitted", "Analytes", "Permutation method", "Iterations per term", "Early stopping", "Seed",
              "lmPerm version", "Self-test", "Pairwise comparisons", "Multiple testing", "Perm_P = 0",
              "LM_P column", "Tests stopped at 51 (p = 1)", "Covariates dropped automatically",
              "Random seeds", "R version"),
  Value = c(model_label("A"), if (!is.null(MODELS$B)) model_label("B") else "not run",
            paste(lev, collapse = ", "), paste(plot_order(lev), collapse = ", "),
            INTERACTION_COVARIATE %||% "none",
            paste0("Separately in each region (", paste(levels(long$Region), collapse = ", "),
                   "); one row per case; observed values only (no imputation)"),
            paste0(n_distinct(long$Analyte), " with >= ", MIN_DETECT, "% detection"),
            "lmPerm::lmp, perm = 'Prob' (random pairwise exchanges of the response), unique SS",
            N_PERM, paste("Disabled (Ca =", CA_TINY, ")"), SEED,
            as.character(packageVersion("lmPerm")), selftest_txt,
            "Indicator terms: every group against every earlier group in GROUP_LEVELS (labelled 'later vs earlier')",
            "Benjamini-Hochberg across analytes within each region x term or region x comparison",
            paste0("No permuted statistic exceeded the observed one: p < 1/", N_PERM),
            "Ordinary F-test p for the same term, for reference only", n_early,
            "A covariate with a single value in a region, or Plate in a single-plate region, is left out there (Covariates column)",
            paste0("One per fit, from Seed + model + region + analyte; rows sorted by Case_ID. ",
                   "Same data + settings = identical p-values on every run"),
            paste0(R.version$major, ".", R.version$minor, "; RNG ", paste(RNGkind(), collapse = " / "))))

data_used <- long |> select(Region, Plate, Sample_ID, Case_ID, Group, all_of(COVARIATES), Analyte, Value, y) |>
  arrange(Analyte, Region, Group, Case_ID)

write_workbook(c(list(Settings = settings, Summary = summ), out, list(Data_Used = data_used)),
               FILES$lmperm, width = 14)
message("Wrote ", FILES$lmperm)
print(summ |> filter(FDR_lt_05 > 0) |> select(Model, Region, Contrast, FDR_lt_05, Significant))
