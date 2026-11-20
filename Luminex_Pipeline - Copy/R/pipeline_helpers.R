# =============================================================================
# Shared code for the Luminex pipeline. You should not need to edit this file;
# study-specific settings live in study_config.R.
#
# Every step script starts by sourcing this file, which
#   1. checks the working directory is the pipeline folder,
#   2. loads study_config.R and checks it for common mistakes,
#   3. defines the readers, group handling and model code all steps share.
# =============================================================================

if (!file.exists("study_config.R") || !dir.exists("R"))
  stop("Run the pipeline from its own folder. In RStudio open Luminex_Pipeline.Rproj, ",
       "or run setwd(\"path/to/Luminex_Pipeline\") first.", call. = FALSE)

# ---- Packages -------------------------------------------------------------------
need_packages <- function(pkgs, attach = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing R package(s): ", paste(miss, collapse = ", "), "\nInstall with:\n  install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE)
  if (attach) suppressPackageStartupMessages(for (p in pkgs) library(p, character.only = TRUE))
}
need_packages(c("readxl", "dplyr", "tidyr", "stringr", "purrr", "openxlsx", "ggplot2"))

# If the conflicted package is active (e.g. from .Rprofile or alongside Seurat /
# Bioconductor), let base, stats and dplyr win every name clash.
if ("conflicted" %in% loadedNamespaces()) {
  conflicted::conflict_prefer_all("stats", quiet = TRUE)
  conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
  conflicted::conflict_prefer_all("base",  quiet = TRUE)
}

# ---- Load and check the settings -------------------------------------------------
source("study_config.R", local = FALSE)

check_config <- function() {
  err <- character(0)
  if (!dir.exists(DATA_DIR)) err <- c(err, paste0("DATA_DIR folder '", DATA_DIR, "' does not exist."))
  if (!BCA_UNITS %in% c("ug/mL", "mg/mL")) err <- c(err, "BCA_UNITS must be \"ug/mL\" or \"mg/mL\".")
  out_names <- names(CASE_INFO_COLUMNS)
  if (is.null(out_names) || any(out_names == ""))
    err <- c(err, "CASE_INFO_COLUMNS needs names, e.g. c(Diagnosis = \"Dx\", Age = \"Age\").")
  bad <- out_names[out_names != make.names(out_names)]
  if (length(bad)) err <- c(err, paste0("CASE_INFO_COLUMNS output names must have no spaces or symbols: ",
                                         paste(bad, collapse = ", ")))
  reserved <- c("Plate", "Well", "Sample_ID", "Case_ID", "Region", "Group", "Value", "y", "Analyte")
  if (any(out_names %in% reserved))
    err <- c(err, paste0("These names are used internally; rename them in CASE_INFO_COLUMNS: ",
                         paste(intersect(out_names, reserved), collapse = ", ")))
  if (!GROUP_COLUMN %in% out_names)
    err <- c(err, paste0("GROUP_COLUMN '", GROUP_COLUMN, "' is not one of the CASE_INFO_COLUMNS output names (",
                         paste(out_names, collapse = ", "), ")."))
  miss_cov <- setdiff(COVARIATES, out_names)
  if (length(miss_cov)) err <- c(err, paste0("COVARIATES not found in CASE_INFO_COLUMNS: ", paste(miss_cov, collapse = ", ")))
  if (GROUP_COLUMN %in% COVARIATES) err <- c(err, "GROUP_COLUMN cannot also be a covariate.")
  if (!is.null(INTERACTION_COVARIATE) && !INTERACTION_COVARIATE %in% COVARIATES)
    err <- c(err, "INTERACTION_COVARIATE must also be listed in COVARIATES.")
  if (!is.null(GROUP_LEVELS) && anyDuplicated(GROUP_LEVELS)) err <- c(err, "GROUP_LEVELS has duplicates.")
  if (!is.null(GROUP_RULES) && is.null(names(GROUP_RULES)))
    err <- c(err, "GROUP_RULES needs group names, e.g. c(CTRL = \"(?i)control\", CASE = \".\").")
  if (length(err)) stop("Problem(s) in study_config.R:\n  - ", paste(err, collapse = "\n  - "), call. = FALSE)
}
check_config()
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

out_path <- function(f) file.path(OUTPUT_DIR, f)
FILES <- list(
  norm   = out_path("Luminex_normalized.xlsx"),
  qc     = out_path("Luminex_QC.xlsx"),
  qc_pdf = out_path("Luminex_QC_plots.pdf"),
  dist   = out_path("Distribution_Check.xlsx"),
  lmperm = out_path("Luminex_lmPerm_Results.xlsx"),
  mv     = out_path("Luminex_Multivariate.xlsx"),
  mv_pdf = out_path("Luminex_Multivariate_plots.pdf"))

need_file <- function(f, step) if (!file.exists(f))
  stop("Can't find ", f, ". Run ", step, " first.", call. = FALSE)

has_regions <- function() length(REGIONS) > 0
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Bio-Plex and BCA readers ------------------------------------------------------
plate_num <- function(f) as.integer(str_match(basename(f), "(?i)plate[ _-]?(\\d+)")[, 2])
num_of    <- function(x) suppressWarnings(as.numeric(str_remove(str_trim(x), "^\\*")))

# Pair every Bio-Plex export with its BCA file by plate number.
find_plate_files <- function() {
  files <- list.files(DATA_DIR, full.names = TRUE)
  files <- files[!str_detect(basename(files), "^~\\$")]            # Excel lock files
  lx  <- grep(LUMINEX_PATTERN, files, value = TRUE, perl = TRUE)
  bca <- grep(BCA_PATTERN, files, value = TRUE, perl = TRUE)
  lx  <- setdiff(lx, bca)
  if (!length(lx)) stop("No Bio-Plex files found in '", DATA_DIR, "'. File names need 'plex' and ",
                        "'Plate <number>' (see LUMINEX_PATTERN in study_config.R). Files there:\n  ",
                        paste(basename(files), collapse = "\n  "), call. = FALSE)
  if (anyNA(plate_num(lx)) || anyNA(plate_num(bca)))
    stop("Could not read a plate number from: ",
         paste(basename(c(lx[is.na(plate_num(lx))], bca[is.na(plate_num(bca))])), collapse = ", "), call. = FALSE)
  for (v in list(lx, bca)) if (anyDuplicated(plate_num(v)))
    stop("Two files share a plate number: ",
         paste(basename(v[plate_num(v) %in% plate_num(v)[duplicated(plate_num(v))]]), collapse = ", "), call. = FALSE)
  lx <- lx[order(plate_num(lx))]
  b  <- bca[match(plate_num(lx), plate_num(bca))]
  if (anyNA(b)) stop("No BCA file for plate(s) ", paste(plate_num(lx)[is.na(b)], collapse = ", "), call. = FALSE)
  tibble(Plate = plate_num(lx), Luminex = lx, BCA = b)
}

# A Bio-Plex sheet holds two blocks: replicate-averaged rows (with Sample IDs),
# then one row per well. Returns both as long tables, or NULL if the sheet is absent.
read_bioplex_sheet <- function(path, sheet) {
  if (!sheet %in% excel_sheets(path)) return(NULL)
  raw <- as.data.frame(suppressMessages(read_excel(path, sheet = sheet, col_names = FALSE, col_types = "text")))
  hdrs <- which(raw[[1]] %in% "Type" & raw[[2]] %in% "Well")
  if (!length(hdrs)) stop("Sheet '", sheet, "' in ", basename(path),
                          " has no 'Type | Well' header row. Is this a Bio-Plex Manager export?", call. = FALSE)
  block <- function(h, stop_at) {
    arow   <- as.character(raw[h - 1, ])
    a_cols <- which(str_detect(arow, "\\(\\d+\\)\\s*$"))
    if (!length(a_cols)) stop("No analyte names like 'IL-6 (25)' above row ", h, " of sheet '", sheet,
                              "' in ", basename(path), call. = FALSE)
    rows <- which(str_detect(coalesce(raw[[1]], ""), "^(B|S\\d+|C\\d+|X\\d+)$"))
    rows <- rows[rows > h & rows < stop_at]
    out <- raw[rows, c(1, 2, if (ncol(raw) >= 4) 4 else 2)]
    names(out) <- c("Type", "Well", "Sample_ID")
    vals <- raw[rows, a_cols, drop = FALSE]
    names(vals) <- str_trim(arow[a_cols])
    bind_cols(as_tibble(out), as_tibble(vals)) |>
      mutate(Sample_ID = str_trim(Sample_ID)) |>
      pivot_longer(-c(Type, Well, Sample_ID), names_to = "Analyte", values_to = "Raw") |>
      mutate(Raw = str_trim(Raw))
  }
  list(summary = block(hdrs[1], if (length(hdrs) > 1) hdrs[2] else Inf),
       wells   = if (length(hdrs) > 1) block(hdrs[2], Inf) else NULL)
}

read_bca <- function(path) {
  raw <- suppressMessages(read_excel(path, col_names = FALSE, col_types = "text"))
  hdr <- which(apply(raw, 1, function(r) any(str_trim(r) %in% "Sample ID")))[1]
  if (is.na(hdr)) stop("No 'Sample ID' header in BCA file ", basename(path), call. = FALSE)
  nm <- str_trim(as.character(unlist(raw[hdr, ])))
  conc <- which(str_detect(coalesce(nm, ""), "^Actual"))[1]
  if (is.na(conc)) stop("No 'Actual Concentration' column in BCA file ", basename(path), call. = FALSE)
  out <- tibble(Sample_ID = str_trim(raw[[which(nm == "Sample ID")[1]]][-(1:hdr)]),
                BCA       = suppressWarnings(as.numeric(raw[[conc]][-(1:hdr)]))) |>
    filter(!is.na(Sample_ID), Sample_ID != "")
  dup <- unique(out$Sample_ID[duplicated(out$Sample_ID)])
  if (length(dup)) {
    warning(basename(path), ": Sample ID(s) listed more than once, values averaged: ",
            paste(dup, collapse = ", "), call. = FALSE)
    out <- out |> group_by(Sample_ID) |> summarise(BCA = mean(BCA, na.rm = TRUE), .groups = "drop")
  }
  out
}

# ---- Sample IDs ------------------------------------------------------------------------
is_bridge <- function(sample_id) if (is.null(BRIDGE_PATTERN)) rep(FALSE, length(sample_id)) else
  str_detect(coalesce(sample_id, ""), BRIDGE_PATTERN)

# "PGEQ-STG" -> Case "PGEQ", Region "STG". Without REGIONS, the Sample ID is the case.
parse_sample_ids <- function(sample_id) {
  if (!has_regions()) return(tibble(Case_ID = sample_id, Region = "All"))
  m <- str_match(sample_id, "^(.*?)[-_ ]+([^-_ ]+)$")
  tibble(Case_ID = coalesce(m[, 2], sample_id), Region = m[, 3])
}

# Case IDs match after removing spaces / hyphens / underscores and ignoring case,
# so "BEMH-0016", "BEMH 0016" and "bemh0016" are the same case.
id_key <- function(x) toupper(str_remove_all(str_trim(as.character(x)), "[^A-Za-z0-9]"))

# ---- Case information file ----------------------------------------------------------------
clean_text <- function(x) {
  x <- str_replace_all(as.character(x), "\u00a0", " ")     # non-breaking spaces from Excel / Airtable
  x <- str_squish(x)
  replace(x, x %in% c("", "NA", "N/A", "n/a", "NaN", "-"), NA)
}

read_case_info <- function() {
  if (is.null(CASE_INFO_FILE)) return(NULL)
  if (!file.exists(CASE_INFO_FILE)) stop("CASE_INFO_FILE not found: ", CASE_INFO_FILE, call. = FALSE)
  meta <- if (str_detect(CASE_INFO_FILE, "(?i)\\.csv$")) {
    # Read as text first: CSVs exported from Excel / Airtable often carry a
    # byte-order mark and non-breaking spaces that make read.csv stop early.
    # Cleaned byte by byte so it behaves the same on every computer and locale.
    b <- readBin(CASE_INFO_FILE, "raw", file.info(CASE_INFO_FILE)$size)
    if (length(b) >= 3 && all(b[1:3] == as.raw(c(0xef, 0xbb, 0xbf)))) b <- b[-(1:3)]   # byte-order mark
    nb <- which(b[-length(b)] == as.raw(0xc2) & b[-1] == as.raw(0xa0))                 # UTF-8 nbsp
    if (length(nb)) { b[nb] <- as.raw(0x20); b <- b[-(nb + 1)] }
    b[b == as.raw(0xa0)] <- as.raw(0x20)                                                # Latin-1 nbsp
    txt <- rawToChar(b); Encoding(txt) <- "UTF-8"
    read.csv(text = txt, check.names = FALSE, colClasses = "character", na.strings = character(0))
  } else as.data.frame(read_excel(CASE_INFO_FILE, sheet = CASE_INFO_SHEET, col_types = "text"))
  names(meta) <- str_trim(str_remove(names(meta), "^\ufeff"))
  need <- c(CASE_ID_COLUMN, unname(CASE_INFO_COLUMNS))
  if (length(setdiff(need, names(meta))))
    stop("Column(s) not found in ", basename(CASE_INFO_FILE), ": ",
         paste(setdiff(need, names(meta)), collapse = ", "),
         "\nColumns in the file are: ", paste(names(meta), collapse = ", "), call. = FALSE)
  meta <- meta[, need] |> setNames(c("meta_id", names(CASE_INFO_COLUMNS))) |>
    mutate(across(everything(), clean_text), .key = id_key(meta_id)) |>
    filter(!is.na(.key), .key != "") |> distinct()
  dup <- meta$meta_id[duplicated(meta$.key)]
  if (length(dup)) stop("Case ID(s) appear more than once in ", basename(CASE_INFO_FILE),
                        " with different values: ", paste(unique(dup), collapse = ", "), call. = FALSE)
  # numbers stay numbers
  for (v in names(CASE_INFO_COLUMNS)) {
    num <- suppressWarnings(as.numeric(meta[[v]]))
    if (all(is.na(meta[[v]]) == is.na(num))) meta[[v]] <- num
  }
  meta
}

# ---- Groups ---------------------------------------------------------------------------------
recode_group <- function(label, case_id) {
  label <- clean_text(label)
  g <- if (is.null(GROUP_RULES)) label else {
    out <- rep(NA_character_, length(label))
    for (k in seq_along(GROUP_RULES)) {
      hit <- is.na(out) & !is.na(label) & str_detect(label, GROUP_RULES[[k]])
      out[hit] <- names(GROUP_RULES)[k]
    }
    out
  }
  if (length(GROUP_OVERRIDES)) {
    ov <- GROUP_OVERRIDES[match(id_key(case_id), id_key(names(GROUP_OVERRIDES)))]
    g <- ifelse(is.na(ov), g, ov)
  }
  g
}

group_levels <- function(present) {
  present <- unique(na.omit(as.character(present)))
  if (is.null(GROUP_LEVELS)) return(sort(present))
  extra <- setdiff(present, GROUP_LEVELS)
  if (length(extra)) stop("Group(s) found in the data but not in GROUP_LEVELS: ", paste(extra, collapse = ", "),
                          "\nAdd them to GROUP_LEVELS (or fix GROUP_RULES) in study_config.R.", call. = FALSE)
  GROUP_LEVELS[GROUP_LEVELS %in% present]
}

plot_order <- function(levels) {
  if (is.null(GROUP_PLOT_ORDER)) return(levels)
  c(GROUP_PLOT_ORDER[GROUP_PLOT_ORDER %in% levels], setdiff(levels, GROUP_PLOT_ORDER))
}

OKABE_ITO <- c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442", "#999999")
group_colors <- function(levels) {
  cols <- setNames(rep(NA_character_, length(levels)), levels)
  if (!is.null(GROUP_COLORS)) cols[intersect(levels, names(GROUP_COLORS))] <- GROUP_COLORS[intersect(levels, names(GROUP_COLORS))]
  free <- setdiff(OKABE_ITO, cols)
  if (sum(is.na(cols)) > length(free)) free <- c(free, grDevices::hcl.colors(sum(is.na(cols)), "Dark 3"))
  cols[is.na(cols)] <- free[seq_len(sum(is.na(cols)))]
  cols
}

# Every pair once, as "later vs earlier" in statistical order.
group_pairs <- function(levels) {
  if (length(levels) < 2) return(tibble(First = character(), Second = character(), Contrast = character()))
  cb <- combn(levels, 2)
  tibble(First = cb[2, ], Second = cb[1, ], Contrast = paste(cb[2, ], "vs", cb[1, ]))
}

# ---- Load the normalized workbook as a long table ------------------------------------------------
analyte_columns <- function(norm) {
  p <- match("Protein_pgmL", names(norm))
  if (is.na(p)) stop("The normalized workbook has no Protein_pgmL column. Was it made by 01_normalize.R?", call. = FALSE)
  names(norm)[seq.int(p + 1, length.out = ncol(norm) - p)]
}

as_covariate <- function(x) {
  x <- clean_text(x)
  num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x) == is.na(num))) num else factor(x)
}

# Reads the normalized sheet and returns one row per sample x analyte with:
#   Value (NA where Bio-Plex gave text), Code (that text), y = log2(Value),
#   Group, Region, Plate and each covariate. Bridges and EXCLUDE_CASES removed.
# observed_only = TRUE keeps measured values only (the univariate analyses).
load_long <- function(min_detect = MIN_DETECT, observed_only = TRUE, complete_covariates = TRUE) {
  need_file(FILES$norm, "01_normalize.R")
  need_file(FILES$qc, "02_qc.R")
  norm <- read_excel(FILES$norm, sheet = "Normalized", col_types = "text")
  analytes <- analyte_columns(norm)
  detect <- read_excel(FILES$qc, sheet = "Detection")
  keep <- detect |> filter(Pct_Detected >= min_detect) |> pull(Analyte)
  analytes_used <- intersect(analytes, keep)
  if (!length(analytes_used)) stop("No analyte reaches ", min_detect, "% detection.", call. = FALSE)

  if (!GROUP_COLUMN %in% names(norm)) stop("Column '", GROUP_COLUMN, "' is not in the normalized workbook.", call. = FALSE)
  long <- norm |>
    select(Plate, Well, Sample_ID, Case_ID, Region, any_of(names(CASE_INFO_COLUMNS)), Protein_pgmL,
           all_of(analytes_used)) |>
    pivot_longer(all_of(analytes_used), names_to = "Analyte", values_to = "Raw") |>
    mutate(Value = suppressWarnings(as.numeric(Raw)),
           Code  = if_else(is.na(Value), str_trim(Raw), NA_character_),
           Value = if_else(Value > 0, Value, NA_real_),
           Plate = factor(as.integer(Plate)),
           Protein_pgmL = as.numeric(Protein_pgmL),
           Group = recode_group(.data[[GROUP_COLUMN]], Case_ID),
           y     = log2(Value)) |>
    filter(!is_bridge(Sample_ID), !id_key(Case_ID) %in% id_key(EXCLUDE_CASES))
  for (v in COVARIATES) long[[v]] <- as_covariate(long[[v]])

  if (has_regions()) {
    other <- setdiff(unique(long$Region), REGIONS)
    if (length(other)) message("Note: samples with region(s) not in REGIONS are left out: ",
                               paste(other, collapse = ", "))
    long <- long |> filter(Region %in% REGIONS) |> mutate(Region = factor(Region, levels = REGIONS))
  } else long <- long |> mutate(Region = factor("All"))

  unassigned <- long |> filter(is.na(Group)) |> distinct(Case_ID, .data[[GROUP_COLUMN]])
  if (nrow(unassigned)) warning(nrow(unassigned), " case(s) have no group and are left out: ",
                                paste(head(paste0(unassigned$Case_ID, " (", unassigned[[GROUP_COLUMN]], ")"), 10),
                                      collapse = ", "), if (nrow(unassigned) > 10) " ...", call. = FALSE)
  long <- long |> filter(!is.na(Group))
  if (complete_covariates && length(COVARIATES)) {
    n0 <- n_distinct(long$Case_ID)
    long <- long |> filter(if_all(all_of(COVARIATES), ~ !is.na(.x)))
    lost <- n0 - n_distinct(long$Case_ID)
    if (lost) message("Note: ", lost, " case(s) dropped for a missing covariate (",
                      paste(COVARIATES, collapse = ", "), ").")
  }
  lev <- group_levels(long$Group)
  long <- long |> mutate(Group = factor(Group, levels = lev))
  if (observed_only) long <- long |> filter(!is.na(Value))

  dups <- long |> count(Region, Analyte, Case_ID) |> filter(n > 1)
  if (nrow(dups)) stop("Some cases have more than one sample in a region, e.g. ",
                       paste(head(unique(dups$Case_ID), 5), collapse = ", "),
                       ". Each case should appear once per region.", call. = FALSE)
  attr(long, "detect") <- detect
  long
}

design_table <- function(long) long |> distinct(Region, Case_ID, Group) |>
  count(Region, Group) |> pivot_wider(names_from = Group, values_from = n, values_fill = 0)

# ---- Model terms that can be fitted in a given subset -----------------------------------------------
# Drops covariates with a single value (e.g. all male) and Plate when only one is present.
usable_terms <- function(d) {
  covs <- COVARIATES[vapply(COVARIATES, function(v) {
    x <- d[[v]]; if (is.factor(x)) nlevels(droplevels(x)) > 1 else length(unique(x)) > 1 }, logical(1))]
  if (ADJUST_FOR_PLATE && nlevels(droplevels(d$Plate)) > 1) covs <- c(covs, "Plate")
  covs
}

# ---- lmPerm ----------------------------------------------------------------------------------------
load_lmperm <- function() {
  ok <- requireNamespace("lmPerm", quietly = TRUE) && packageVersion("lmPerm") == "2.1.0"
  if (!ok && AUTO_INSTALL_LMPERM) {
    message("Installing lmPerm 2.1.0 from the CRAN archive (needs Rtools on Windows) ...")
    if ("lmPerm" %in% loadedNamespaces()) unloadNamespace("lmPerm")
    install.packages("https://cran.r-project.org/src/contrib/Archive/lmPerm/lmPerm_2.1.0.tar.gz",
                     repos = NULL, type = "source")
  }
  if (!requireNamespace("lmPerm", quietly = TRUE) || packageVersion("lmPerm") != "2.1.0")
    stop("lmPerm 2.1.0 is required (2.1.6 gives miscalibrated p-values when a factor has 3+ levels). Install with:\n",
         "  install.packages(\"https://cran.r-project.org/src/contrib/Archive/lmPerm/lmPerm_2.1.0.tar.gz\",\n",
         "                   repos = NULL, type = \"source\")\nthen restart R.", call. = FALSE)
  suppressPackageStartupMessages(library(lmPerm))
}

CA_TINY <- 1e-9   # disables lmPerm's adaptive early stopping

# Reproducibility: every fit gets its own seed, built from SEED plus labels such as
# model, region and analyte. Its permutations then depend only on its own data and
# settings, not on how many other fits ran before it or in what order, so adding a
# case, changing MIN_DETECT or rerunning one region leaves every other result unchanged.
fit_seed <- function(...) {
  h <- 0
  for (b in utf8ToInt(enc2utf8(paste(SEED, ..., sep = "|")))) h <- (h * 31 + b) %% 2147483647
  as.integer(h)
}

run_lmp <- function(form, d, n_perm = N_PERM)
  lmPerm::lmp(form, data = d, perm = "Prob", maxIter = n_perm, Ca = CA_TINY, settings = FALSE)

# Interaction labels in a fixed order: the covariate first ("Age:Group", "Age:g2").
norm_label <- function(x, iv) vapply(strsplit(x, ":"), function(v)
  paste(c(v[v %in% iv], v[!v %in% iv]), collapse = ":"), "")

# Permutation p, iteration count, unique SS and ordinary F p for every term.
lmp_terms <- function(fit, iv = INTERACTION_COVARIATE) {
  tl   <- attr(fit$terms, "term.labels")
  asgn <- fit$assign[fit$qr$pivot][seq_len(fit$rank)]
  lab  <- c("(Intercept)", tl)[1 + unique(asgn)]
  ss   <- tapply(fit$effects[seq_len(fit$rank)]^2, asgn, sum)
  df   <- tapply(asgn, asgn, length)
  rdf  <- df.residual(fit); ms_res <- sum(fit$residuals^2) / rdf
  tibble(Term = norm_label(lab, iv), Df = as.integer(df), SumSq = as.vector(ss),
         Perm_P = as.vector(fit$perm$P[, 1]), Iter = as.integer(fit$perm$Mn[, 1])) |>
    mutate(LM_P = pf((SumSq / Df) / ms_res, Df, rdf, lower.tail = FALSE)) |>
    filter(Term != "(Intercept)")
}

# Fit one outcome (column y) in one subset under Model A or B.
#   1. factor coding     -> term p-values (Group, covariates, Plate, [IV:Group])
#   2. indicator coding  -> each pairwise comparison as its own 1-df term,
#                           one fit per reference group (all groups but the last)
# Groups with fewer than MIN_GROUP cases are left out of that fit.
fit_group_model <- function(d, interaction = FALSE, n_perm = N_PERM) {
  iv <- if (interaction) INTERACTION_COVARIATE else NULL
  d <- d |> filter(!is.na(y)) |> arrange(Case_ID)   # fixed row order: permutations act on row positions
  ng0 <- table(droplevels(d$Group))
  small <- names(ng0)[ng0 < MIN_GROUP]
  d <- d |> filter(!Group %in% small) |> mutate(Group = droplevels(Group))
  lev <- levels(d$Group)
  if (nrow(d) < MIN_N || length(lev) < 2) return(NULL)
  ng <- table(d$Group)
  covs <- usable_terms(d)
  if (interaction && !iv %in% covs) return(NULL)
  for (v in covs) if (is.numeric(d[[v]])) d[[v]] <- d[[v]] - mean(d[[v]])   # centre within the subset
  for (v in covs) if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  codes <- setNames(paste0("g", seq_along(lev)), lev)        # safe names for any group label
  for (g in lev) d[[codes[[g]]]] <- as.numeric(d$Group == g)

  rhs_cov <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""
  f_terms <- as.formula(paste("y ~ Group", rhs_cov, if (interaction) paste0("+ ", iv, ":Group") else ""))
  terms <- lmp_terms(run_lmp(f_terms, d, n_perm), iv) |>
    mutate(N = nrow(d), Groups = paste(lev, collapse = ", "),
           Groups_Dropped = paste(small, collapse = ", "), Covariates = paste(covs, collapse = ", "))

  pair_rows <- list(); slope_rows <- list()
  for (ref in lev[-length(lev)]) {
    later <- lev[seq(match(ref, lev) + 1, length(lev))]
    ind <- codes[setdiff(lev, ref)]
    f <- as.formula(paste("y ~", paste(c(ind, covs, if (interaction) paste0(iv, ":", ind)), collapse = " + ")))
    tp <- lmp_terms(run_lmp(f, d, n_perm), iv)
    cf <- coef(lm(f, data = d)); names(cf) <- norm_label(names(cf), iv)
    for (g in later) {
      code <- codes[[g]]
      pair_rows[[length(pair_rows) + 1]] <- tibble(
        Contrast = paste(g, "vs", ref), Log2FC = unname(cf[code]),
        Perm_P = tp$Perm_P[tp$Term == code], Iter = tp$Iter[tp$Term == code],
        N_First = as.integer(ng[[g]]), N_Second = as.integer(ng[[ref]]))
      if (interaction) {
        it <- paste0(iv, ":", code)
        slope_rows[[length(slope_rows) + 1]] <- tibble(
          Contrast = paste(g, "vs", ref), Slope_Diff_log2_per_unit = unname(cf[it]),
          Perm_P = tp$Perm_P[tp$Term == it], Iter = tp$Iter[tp$Term == it])
      }
    }
  }
  list(terms = terms,
       pairs = bind_rows(pair_rows) |> mutate(Fold_Change = 2^Log2FC, .after = Log2FC),
       slopes = if (interaction) bind_rows(slope_rows) else NULL)
}

# ---- Output helpers -----------------------------------------------------------------------------------
HEADER_STYLE <- function() createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                                       fgFill = "#D9D9D9", border = "Bottom", wrapText = TRUE)

# sheets: named list of data frames. highlight: column holding TRUE for rows to shade.
write_workbook <- function(sheets, file, highlight = "Signif", width = 14) {
  wb <- createWorkbook()
  hs <- HEADER_STYLE(); sig <- createStyle(fgFill = "#C6EFCE")
  for (nm in names(sheets)) {
    d <- sheets[[nm]]
    if (is.null(d)) next
    if (!nrow(d)) d <- tibble(Note = "nothing to report")
    addWorksheet(wb, nm)
    writeData(wb, nm, d, headerStyle = hs)
    freezePane(wb, nm, firstRow = TRUE)
    setColWidths(wb, nm, seq_len(ncol(d)), width)
    if (highlight %in% names(d)) {
      rows <- which(d[[highlight]] %in% TRUE)
      if (length(rows)) addStyle(wb, nm, sig, rows = rows + 1, cols = seq_len(ncol(d)), gridExpand = TRUE, stack = TRUE)
    }
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

chunks <- function(x, n) split(x, ceiling(seq_along(x) / n))
