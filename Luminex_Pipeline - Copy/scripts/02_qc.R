# =============================================================================
# STEP 2  Luminex QC: run this before any group comparison
#
# Input : Bio-Plex exports + BCA files in DATA_DIR (same pairing as step 1)
# Output: output/Luminex_QC.xlsx + output/Luminex_QC_plots.pdf
#   Std_Curves      one row per plate x analyte: 5PL fit stats, LLOQ/ULOQ,
#                   number of standards recovering within tolerance
#   Std_Recovery    every standard: expected vs observed, % recovery, pass/fail
#   Detection       per analyte: % of samples in range, overall and per plate
#                   (later steps use this to decide which analytes to analyse)
#   Curve_vs_Detect Spearman correlation of curve quality with detection rate
#   Bead_Counts     per plate x analyte bead statistics, low-bead wells
#   Kit_Controls    C1/C2 across plates, inter-plate CV
#   Bridges         bridge samples across plates: CV, plate factors
#                   (skipped if BRIDGE_PATTERN is NULL or there is one plate)
#   Plate_Effects   per-plate shift from bridges; technical vs biological spread
#   Sample_Outliers per-sample median z-score across analytes
# =============================================================================

if (!file.exists("R/pipeline_helpers.R"))
  stop("Run this from the Luminex_Pipeline folder: open Luminex_Pipeline.Rproj in RStudio, ",
       "or setwd(\"path/to/Luminex_Pipeline\") first.", call. = FALSE)
source("R/pipeline_helpers.R")
to_pgmL <- switch(BCA_UNITS, "ug/mL" = 1e6, "mg/mL" = 1e9)

# 5PL / 4PL equations from the Standard Curve sheet. Bio-Plex puts the analyte
# names only inside the plot images and omits curves it could not fit, so each
# equation is matched back to its analyte by comparing predicted with observed FI.
fpl <- function(x, p) p$a + (p$d - p$a) / (1 + (x / p$c)^p$b)^p$g

read_curves <- function(path, std_obs) {
  raw <- suppressMessages(read_excel(path, sheet = "Standard Curve",
                                     col_names = FALSE, col_types = "text"))
  mat <- as.matrix(raw)
  idx <- which(!is.na(mat) & str_detect(mat, "^(Std\\. Curve:|FitProb)"), arr.ind = TRUE)
  if (!nrow(idx)) return(NULL)
  idx <- idx[order(idx[, "col"], idx[, "row"]), , drop = FALSE]
  recs <- list()
  for (k in seq_len(nrow(idx))) {
    txt <- mat[idx[k, "row"], idx[k, "col"]]
    if (str_detect(txt, "^Std")) {
      recs[[length(recs) + 1]] <- list(row = idx[k, "row"], col = idx[k, "col"], eq = txt, fit = NA)
    } else if (length(recs) && recs[[length(recs)]]$col == idx[k, "col"]) {
      recs[[length(recs)]]$fit <- txt
    }
  }
  if (!length(recs)) return(NULL)
  cv <- map_dfr(recs, as_tibble) |> arrange(row, col)
  n <- "(-?[0-9.eE+-]+)"
  g <- function(pat) as.numeric(str_match(cv$eq, pat)[, 2])
  cv <- cv |> transmute(
    a = g(paste0("FI = ", n)),
    d = g(paste0("FI = -?[0-9.eE+-]+ \\+ \\(", n)),
    c = g(paste0("Conc / ", n)),
    b = g(paste0("Conc / -?[0-9.eE+-]+\\)\\^", n)),
    g = coalesce(g(paste0("\\)\\)\\^", n, "\\s*$")), 1),
    Model    = if_else(str_detect(eq, "\\)\\)\\^-?[0-9.]+\\s*$"), "5PL", "4PL"),
    FitProb  = as.numeric(str_match(fit, "FitProb\\. = ([0-9.]+)")[, 2]),
    ResVar   = as.numeric(str_match(fit, "ResVar\\. = ([0-9.]+)")[, 2]))
  an <- unique(std_obs$Analyte)
  an <- an[order(as.integer(str_match(an, "\\((\\d+)\\)\\s*$")[, 2]))]   # bead-region order
  nC <- nrow(cv); nA <- length(an)
  if (nC > nA) return(NULL)
  cost <- matrix(Inf, nC, nA)
  for (i in seq_len(nC)) for (j in seq_len(nA)) {
    o <- std_obs[std_obs$Analyte == an[j] & !is.na(std_obs$FI) & std_obs$FI > 0, ]
    pr <- fpl(o$Exp_Conc, cv[i, ]); ok <- is.finite(pr) & pr > 0
    if (sum(ok) >= 3) cost[i, j] <- median(abs(log(pr[ok]) - log(o$FI[ok])))
  }
  if (all(is.infinite(cost))) {
    warning("could not match standard curves to analytes; curve stats skipped")
    return(NULL)
  }
  D <- matrix(Inf, nC + 1, nA + 1); D[1, ] <- 0
  for (i in 1:nC) for (j in i:nA) D[i + 1, j + 1] <- min(D[i + 1, j], D[i, j] + cost[i, j])
  asg <- integer(nC); j <- nA
  for (i in nC:1) {
    while (j > i && D[i + 1, j + 1] == D[i + 1, j]) j <- j - 1   # j >= i keeps the index valid
    asg[i] <- j; j <- j - 1
  }
  err  <- cost[cbind(seq_len(nC), asg)]
  poor <- !is.finite(err) | err > 0.5          # equation does not reproduce that analyte's FI
  if (any(poor)) warning("curve could not be matched to an analyte (", sum(poor),
                         " of ", nC, "); those curve stats are left blank")
  miss <- setdiff(an, an[asg[!poor]])
  if (length(miss)) message("   no curve equation in export for: ", paste(miss, collapse = ", "))
  cv |> mutate(Analyte = if_else(poor, NA_character_, an[asg]), Match_Err = err) |>
    filter(!is.na(Analyte))
}

# ---- Read every plate --------------------------------------------------------
pairs <- find_plate_files()

plates <- map2(pairs$Luminex, pairs$BCA, function(f, bf) {
  p <- plate_num(f)
  message("Plate ", p, ": ", basename(f))
  oc <- read_bioplex_sheet(f, "Obs Conc")
  fi <- read_bioplex_sheet(f, "FI")
  bd <- read_bioplex_sheet(f, "Bead Count")$wells
  ec <- read_bioplex_sheet(f, "Exp Conc")$summary
  cr <- read_bioplex_sheet(f, "Conc in Range")$summary
  if (is.null(oc) || is.null(ec) || is.null(cr) || is.null(fi))
    stop(basename(f), " is missing one of the sheets QC needs: Obs Conc, FI, Exp Conc, Conc in Range.", call. = FALSE)
  if (is.null(oc$wells)) stop(basename(f), ": the Obs Conc sheet has no per-well block.", call. = FALSE)
  if (is.null(bd)) {
    message("   no 'Bead Count' sheet; bead checks skipped for this plate")
    bd <- tibble(Well = character(), Analyte = character(), Raw = character())
  }

  ids <- oc$summary |> filter(str_detect(Type, "^X")) |> distinct(Type, Sample_ID)

  std <- oc$summary |> filter(str_detect(Type, "^S\\d")) |>
    transmute(Plate = p, Type, Analyte, Obs_Conc = num_of(Raw),
              Flag = case_when(is.na(Raw) ~ "blank", str_detect(Raw, "OOR") ~ Raw,
                               str_detect(Raw, "^\\*\\*\\*|^---") ~ "dropped",
                               str_detect(Raw, "^\\*") ~ "extrapolated", TRUE ~ "reported")) |>
    left_join(ec |> transmute(Type, Analyte, Exp_Conc = num_of(Raw)), by = c("Type", "Analyte")) |>
    left_join(cr |> transmute(Type, Analyte,
                              In_Range = !str_detect(coalesce(Raw, "OOR"), "OOR|\\*")),
              by = c("Type", "Analyte")) |>
    mutate(Recovery = 100 * Obs_Conc / Exp_Conc,
           Pass = Flag == "reported" & In_Range & !is.na(Recovery) &
             Recovery >= RECOVERY_MIN & Recovery <= RECOVERY_MAX)

  std_fi <- fi$wells |> filter(str_detect(Type, "^S\\d")) |>
    transmute(Type, Analyte, FI = num_of(Raw)) |>
    left_join(std |> distinct(Type, Analyte, Exp_Conc), by = c("Type", "Analyte")) |>
    group_by(Type, Analyte, Exp_Conc) |> summarise(FI = mean(FI, na.rm = TRUE), .groups = "drop")
  curves <- tryCatch(read_curves(f, std_fi), error = function(e) {
    message("   could not read the Standard Curve sheet (", conditionMessage(e), "); curve stats skipped")
    NULL })

  wells <- oc$wells |> select(-Sample_ID) |>
    transmute(Plate = p, Type, Well, Analyte, Raw,
              Obs_Conc = num_of(Raw),
              Code = if_else(is.na(num_of(Raw)), coalesce(Raw, ""), NA_character_),
              Extrapolated = str_detect(coalesce(Raw, ""), "^\\*\\d")) |>
    left_join(fi$wells |> transmute(Well, Analyte, FI = num_of(Raw)), by = c("Well", "Analyte")) |>
    left_join(bd |> transmute(Well, Analyte, Beads = num_of(Raw)), by = c("Well", "Analyte")) |>
    left_join(ids, by = "Type") |>
    left_join(read_bca(bf) |> transmute(Sample_ID, Protein = BCA), by = "Sample_ID") |>
    mutate(Class = case_when(Type == "B" ~ "blank",
                             str_detect(Type, "^S\\d") ~ "standard",
                             str_detect(Type, "^C\\d") ~ "kit_control",
                             is_bridge(Sample_ID) ~ "bridge",
                             TRUE ~ "sample"),
           Norm = Obs_Conc / (Protein * to_pgmL) * SCALE)
  list(std = std, wells = wells,
       curves = if (is.null(curves)) NULL else curves |> mutate(Plate = p))
})

std    <- map_dfr(plates, "std")
wells  <- map_dfr(plates, "wells")
curves <- map_dfr(plates, "curves")

# ---- 1. Standard curves: limits of quantification ---------------------------
limits <- std |> group_by(Plate, Analyte) |>
  summarise(N_Std_Pass = sum(Pass),
            LLOQ = if (any(Pass)) min(Exp_Conc[Pass]) else NA_real_,
            ULOQ = if (any(Pass)) max(Exp_Conc[Pass]) else NA_real_,
            Median_Recovery = median(Recovery[Pass], na.rm = TRUE), .groups = "drop")
if (nrow(curves) && "FitProb" %in% names(curves))
  limits <- left_join(limits, curves |> select(Plate, Analyte, Model, FitProb, ResVar),
                      by = c("Plate", "Analyte"))
if (!"FitProb" %in% names(limits)) limits <- limits |> mutate(Model = NA_character_, FitProb = NA_real_, ResVar = NA_real_)
limits <- limits |>
  mutate(Curve_Flag = case_when(
    N_Std_Pass < 4                      ~ "fewer than 4 usable standards",
    !is.na(FitProb) & FitProb < 0.05    ~ "poor fit (FitProb < 0.05)",
    TRUE                                ~ ""))

# ---- 2. Detection ------------------------------------------------------------
samples <- wells |> filter(Class %in% c("sample", "bridge")) |>
  left_join(limits |> select(Plate, Analyte, LLOQ, ULOQ), by = c("Plate", "Analyte")) |>
  mutate(Status = case_when(
    !is.na(Beads) & Beads < BEAD_MIN                  ~ "low_beads",
    coalesce(Code, "") %in% c("***", "---")           ~ "missing",
    is.na(LLOQ)                                       ~ "no_curve",
    coalesce(Code, "") == "OOR <"                     ~ "below_LLOQ",
    coalesce(Code, "") == "OOR >"                     ~ "above_ULOQ",
    !is.na(Obs_Conc) & Obs_Conc < LLOQ                ~ "below_LLOQ",
    !is.na(Obs_Conc) & Obs_Conc > ULOQ                ~ "above_ULOQ",
    TRUE                                              ~ "ok"))

detect <- samples |> filter(Class == "sample") |>
  group_by(Analyte) |>
  summarise(N = n(),
            Pct_Detected = 100 * mean(Status == "ok"),
            Pct_Below    = 100 * mean(Status == "below_LLOQ"),
            Pct_Above    = 100 * mean(Status == "above_ULOQ"),
            Pct_LowBead  = 100 * mean(Status == "low_beads"), .groups = "drop") |>
  mutate(Tier = cut(Pct_Detected / 100, c(-Inf, DETECT_TIERS, Inf),
                    labels = c("rarely detected", "partly detected", "well detected"))) |>
  arrange(desc(Pct_Detected))

detect_by_plate <- samples |> filter(Class == "sample") |>
  group_by(Analyte, Plate) |>
  summarise(Pct_Detected = 100 * mean(Status == "ok"), .groups = "drop") |>
  pivot_wider(names_from = Plate, values_from = Pct_Detected, names_prefix = "Plate_")

# ---- 3. Does curve quality predict detection? --------------------------------
curve_detect <- limits |> group_by(Analyte) |>
  summarise(Median_LLOQ = median(LLOQ, na.rm = TRUE),
            Min_FitProb = suppressWarnings(min(FitProb, na.rm = TRUE)),
            Mean_N_Std_Pass = mean(N_Std_Pass),
            Plates_Flagged = sum(Curve_Flag != ""), .groups = "drop") |>
  mutate(across(everything(), ~ replace(.x, is.infinite(.x), NA))) |>
  left_join(detect |> select(Analyte, Pct_Detected, Tier), by = "Analyte") |>
  left_join(samples |> filter(Class == "sample") |> group_by(Analyte) |>
              summarise(Median_FI = median(FI, na.rm = TRUE), .groups = "drop"), by = "Analyte")

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(c(NA, NA))
  t <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(unname(t$estimate), t$p.value)
}
curve_vs_detect <- tibble(
  Comparison = c("Median LLOQ vs % detected",
                 "Min FitProb vs % detected",
                 "Usable standards vs % detected",
                 "Median sample FI vs % detected"),
  Spearman_rho = c(sp(curve_detect$Median_LLOQ, curve_detect$Pct_Detected)[1],
                   sp(curve_detect$Min_FitProb, curve_detect$Pct_Detected)[1],
                   sp(curve_detect$Mean_N_Std_Pass, curve_detect$Pct_Detected)[1],
                   sp(curve_detect$Median_FI, curve_detect$Pct_Detected)[1]),
  P_value = c(sp(curve_detect$Median_LLOQ, curve_detect$Pct_Detected)[2],
              sp(curve_detect$Min_FitProb, curve_detect$Pct_Detected)[2],
              sp(curve_detect$Mean_N_Std_Pass, curve_detect$Pct_Detected)[2],
              sp(curve_detect$Median_FI, curve_detect$Pct_Detected)[2]))

# ---- 4. Bead counts ----------------------------------------------------------
beads <- wells |> filter(!is.na(Beads)) |>
  group_by(Plate, Analyte) |>
  summarise(Min_Beads = min(Beads), Median_Beads = median(Beads),
            N_Wells_Low = sum(Beads < BEAD_MIN), .groups = "drop") |>
  arrange(Min_Beads)
low_bead_wells <- wells |> filter(!is.na(Beads), Beads < BEAD_MIN) |>
  select(Plate, Well, Sample_ID, Analyte, Beads)

# ---- 5. Kit controls ---------------------------------------------------------
kit <- wells |> filter(Class == "kit_control") |>
  group_by(Type, Analyte) |>
  summarise(N_Plates = sum(!is.na(Obs_Conc)),
            Mean = mean(Obs_Conc, na.rm = TRUE),
            CV_pct = 100 * sd(Obs_Conc, na.rm = TRUE) / mean(Obs_Conc, na.rm = TRUE),
            .groups = "drop") |> arrange(desc(CV_pct))

# ---- 6. Bridges and plate effects --------------------------------------------
n_plates <- n_distinct(wells$Plate)
bridge <- samples |> filter(Class == "bridge", Status == "ok", is.finite(Norm))
if (n_plates < 2 && nrow(bridge)) {
  message("One plate only: bridge / plate-effect sections skipped")
  bridge <- bridge[0, ]
}
if (!nrow(bridge)) message("No bridge samples found (BRIDGE_PATTERN = ",
                           if (is.null(BRIDGE_PATTERN)) "NULL" else BRIDGE_PATTERN, "): bridge sections left empty")

suppressWarnings({   # summaries of empty tables (no bridges) warn harmlessly
bridge_cv <- bridge |> group_by(Analyte, Sample_ID) |>
  summarise(N_Plates = n_distinct(Plate),
            Mean_Norm = mean(Norm),
            CV_pct = 100 * sd(Norm) / mean(Norm),
            Max_Fold_Range = max(Norm) / min(Norm), .groups = "drop") |>
  arrange(desc(CV_pct))

# Plate factor: how far each plate sits from the across-plate geometric mean
bridge_ref <- bridge |> group_by(Analyte, Sample_ID) |>
  filter(n_distinct(Plate) == n_plates) |>
  summarise(Ref = exp(mean(log(Norm))), .groups = "drop")

plate_factors <- bridge |> inner_join(bridge_ref, by = c("Analyte", "Sample_ID")) |>
  group_by(Plate, Analyte) |>
  summarise(N_Bridges = n(),
            Plate_Factor = exp(median(log(Ref / Norm))), .groups = "drop") |>
  mutate(Plate_Factor = if_else(N_Bridges >= 2, Plate_Factor, NA_real_))

# Technical (bridge) vs biological (sample) spread, per analyte, on the log scale
spread <- bridge_cv |> group_by(Analyte) |>
  summarise(Technical_CV = median(CV_pct), .groups = "drop") |>
  left_join(samples |> filter(Class == "sample", Status == "ok", is.finite(Norm)) |>
              group_by(Analyte) |>
              summarise(Biological_CV = 100 * sd(Norm) / mean(Norm),
                        N_Samples = n(), .groups = "drop"), by = "Analyte") |>
  mutate(Tech_Over_Bio = Technical_CV / Biological_CV,
         Correction = case_when(
           is.na(Technical_CV)         ~ "no bridge data - cannot assess",
           Technical_CV > 25 & Tech_Over_Bio > 0.5 ~ "plate correction recommended",
           Technical_CV > 25           ~ "plate effect present but small vs biology",
           TRUE                        ~ "no correction needed")) |>
  arrange(desc(Technical_CV))
})

# Whole-plate shift: median plate factor across analytes (a systematic plate offset)
plate_shift <- plate_factors |> filter(!is.na(Plate_Factor)) |>
  group_by(Plate) |>
  summarise(N_Analytes = n(),
            Median_Plate_Factor = median(Plate_Factor),
            IQR_Plate_Factor = IQR(Plate_Factor), .groups = "drop")

# ---- 7. Sample-level outliers -------------------------------------------------
sample_z <- samples |> filter(Class == "sample", Status == "ok", is.finite(Norm), Norm > 0) |>
  group_by(Analyte) |>
  mutate(z = (log10(Norm) - mean(log10(Norm))) / sd(log10(Norm))) |>
  group_by(Plate, Well, Sample_ID) |>
  summarise(N_Analytes = n(), Median_z = median(z, na.rm = TRUE), .groups = "drop") |>
  left_join(wells |> distinct(Plate, Well, Protein), by = c("Plate", "Well")) |>
  mutate(Flag = if_else(abs(Median_z) > 2, "check", "")) |>
  arrange(desc(abs(Median_z)))

# ---- Write workbook -----------------------------------------------------------
wb <- createWorkbook()
hs <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                  fgFill = "#D9D9D9", border = "Bottom", wrapText = TRUE)
add <- function(name, df) {
  if (!nrow(df)) df <- tibble(Note = "nothing to report for this study")
  addWorksheet(wb, name); writeData(wb, name, df, headerStyle = hs)
  addStyle(wb, name, createStyle(fontName = "Arial", fontSize = 10),
           rows = 2:(nrow(df) + 1), cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
  freezePane(wb, name, firstRow = TRUE); setColWidths(wb, name, seq_len(ncol(df)), 14)
}
add("Std_Curves", limits)
add("Std_Recovery", std)
add("Detection", detect |> left_join(detect_by_plate, by = "Analyte"))
add("Curve_vs_Detect", curve_vs_detect)
add("Curve_Detect_Data", curve_detect)
add("Bead_Counts", beads)
add("Low_Bead_Wells", low_bead_wells)
add("Kit_Controls", kit)
add("Bridges", bridge_cv)
add("Plate_Factors", plate_factors)
add("Plate_Effects", spread)
add("Plate_Shift", plate_shift)
add("Sample_Outliers", sample_z)
saveWorkbook(wb, FILES$qc, overwrite = TRUE)

# ---- Plots ---------------------------------------------------------------------
pdf(FILES$qc_pdf, width = 11, height = 8.5)
print(ggplot(curve_detect, aes(Median_LLOQ, Pct_Detected)) +
        geom_point(aes(colour = Tier), size = 2) + scale_x_log10() +
        labs(title = "Detection rate vs assay sensitivity",
             x = "Median LLOQ across plates (pg/mL)", y = "% of samples in range") +
        theme_bw())
print(ggplot(curve_detect |> filter(is.finite(Min_FitProb)),
             aes(Min_FitProb, Pct_Detected)) +
        geom_point(aes(colour = Tier), size = 2) +
        geom_vline(xintercept = 0.05, linetype = 2) +
        labs(title = "Detection rate vs curve fit quality",
             x = "Worst FitProb across plates", y = "% of samples in range") +
        theme_bw())
print(ggplot(detect_by_plate |> pivot_longer(-Analyte, names_to = "Plate",
                                             values_to = "Pct"),
             aes(Plate, reorder(Analyte, Pct), fill = Pct)) +
        geom_tile() + scale_fill_viridis_c(option = "magma") +
        labs(title = "% detected per analyte per plate", y = NULL, fill = "%") +
        theme_bw(base_size = 7))
if (nrow(bridge)) {
  top <- bridge_cv |> slice_head(n = 24) |> pull(Analyte) |> unique()
  print(ggplot(bridge |> filter(Analyte %in% top),
               aes(factor(Plate), Norm, group = Sample_ID, colour = Sample_ID)) +
          geom_line() + geom_point(size = 1) +
          facet_wrap(~ Analyte, scales = "free_y") +
          labs(title = "Bridge samples across plates (most variable analytes)",
               x = "Plate", y = "Normalized") +
          theme_bw(base_size = 7) + theme(legend.position = "bottom"))
  print(ggplot(spread |> filter(is.finite(Technical_CV), is.finite(Biological_CV)),
               aes(Biological_CV, Technical_CV)) +
          geom_abline(slope = 1, linetype = 2) + geom_point(size = 2) +
          labs(title = "Technical (bridge) vs biological spread",
               subtitle = "Points above the line: plate variation exceeds between-sample variation",
               x = "Between-sample CV %", y = "Bridge CV across plates %") +
          theme_bw())
}
print(ggplot(sample_z, aes(reorder(paste(Plate, Sample_ID), Median_z), Median_z,
                           colour = factor(Plate))) +
        geom_point() + geom_hline(yintercept = c(-2, 2), linetype = 2) +
        coord_flip() + labs(title = "Sample-level shift across all analytes",
                            x = NULL, y = "Median z (log10 scale)", colour = "Plate") +
        theme_bw(base_size = 6))
invisible(dev.off())

message("\nWrote ", FILES$qc, " and ", FILES$qc_pdf)
message("Analytes well detected (>", 100 * DETECT_TIERS[2], "%): ",
        sum(detect$Tier == "well detected"),
        " | partly: ", sum(detect$Tier == "partly detected"),
        " | rarely: ", sum(detect$Tier == "rarely detected"))
