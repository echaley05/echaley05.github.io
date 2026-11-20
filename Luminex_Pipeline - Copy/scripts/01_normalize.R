# =============================================================================
# STEP 1  Normalize Bio-Plex concentrations to BCA total protein
#
# Input : every Bio-Plex export + BCA file in DATA_DIR (paired by plate number),
#         and the case information file (optional)
# Output: output/Luminex_normalized.xlsx
#           Normalized      one row per sample: IDs, case info, BCA, analytes
#           Raw_pg_per_mL   the same layout before dividing by protein
#           Group_Check     how each diagnosis label was assigned to a group
#           Metadata_check  cases missing from either file
#           Notes           formula, codes and files used
#
#   Normalized = Obs Conc (pg/mL) / BCA protein (pg/mL) x SCALE   (unitless)
#
# Extrapolated values (Bio-Plex "*24.7") are kept and shaded yellow.
# "OOR <", "OOR >", "***" and "---" are kept as text; nothing is imputed.
# =============================================================================

if (!file.exists("R/pipeline_helpers.R"))
  stop("Run this from the Luminex_Pipeline folder: open Luminex_Pipeline.Rproj in RStudio, ",
       "or setwd(\"path/to/Luminex_Pipeline\") first.", call. = FALSE)
source("R/pipeline_helpers.R")

to_pgmL <- switch(BCA_UNITS, "ug/mL" = 1e6, "mg/mL" = 1e9)
bca_col <- paste0("BCA_", str_remove(BCA_UNITS, "/"))

# ---- Read every plate ------------------------------------------------------------
pairs <- find_plate_files()
plates <- pmap(pairs, function(Plate, Luminex, BCA) {
  oc <- read_bioplex_sheet(Luminex, "Obs Conc")
  if (is.null(oc)) stop("No 'Obs Conc' sheet in ", basename(Luminex), call. = FALSE)
  x <- oc$summary |> filter(str_detect(Type, "^X")) |>
    transmute(Plate, Well, Sample_ID, Analyte, Raw,
              Obs = num_of(Raw),
              Extrapolated = str_detect(coalesce(Raw, ""), "^\\*\\s*[0-9.]"),
              Code = if_else(is.na(Obs), coalesce(Raw, "***"), NA_character_))
  bca <- read_bca(BCA)
  ids <- distinct(x, Sample_ID)
  hit <- match(ids$Sample_ID, bca$Sample_ID)
  loose <- is.na(hit)                                    # fall back to a forgiving match
  hit[loose] <- match(id_key(ids$Sample_ID[loose]), id_key(bca$Sample_ID))
  ids$BCA <- bca$BCA[hit]
  miss <- ids$Sample_ID[is.na(ids$BCA)]
  if (length(miss)) warning("Plate ", Plate, ": no BCA value for ", paste(miss, collapse = ", "), call. = FALSE)
  message(sprintf("Plate %d: %d samples, %d analytes  (%s + %s)", Plate, nrow(ids),
                  n_distinct(x$Analyte), basename(Luminex), basename(BCA)))
  left_join(x, ids, by = "Sample_ID")
})
x <- bind_rows(plates) |>
  mutate(Protein_pgmL = BCA * to_pgmL,
         Norm = Obs / Protein_pgmL * SCALE)

analytes <- unique(x$Analyte)
analytes <- analytes[order(suppressWarnings(as.integer(str_match(analytes, "\\((\\d+)\\)\\s*$")[, 2])))]

# ---- Sample IDs -> case and region ------------------------------------------------
id_tbl <- x |> distinct(Plate, Well, Sample_ID, BCA, Protein_pgmL)
id_tbl <- bind_cols(id_tbl, parse_sample_ids(id_tbl$Sample_ID))
if (has_regions()) {
  odd <- id_tbl |> filter(!is_bridge(Sample_ID), !Region %in% REGIONS)
  if (nrow(odd)) warning(nrow(odd), " sample ID(s) do not end in one of REGIONS (",
                         paste(REGIONS, collapse = ", "), "), e.g. ",
                         paste(head(odd$Sample_ID, 5), collapse = ", "),
                         ". They will be left out of the analyses.", call. = FALSE)
}

# ---- Case information -----------------------------------------------------------------
meta <- read_case_info()
meta_check <- NULL
if (!is.null(meta)) {
  id_tbl <- id_tbl |> mutate(.key = id_key(Case_ID)) |>
    left_join(meta |> select(-meta_id), by = ".key") |> select(-.key)
  no_meta <- id_tbl |> filter(!is_bridge(Sample_ID),
                              if_all(all_of(names(CASE_INFO_COLUMNS)), is.na)) |> distinct(Case_ID) |> pull()
  unused <- meta$meta_id[!meta$.key %in% id_key(id_tbl$Case_ID)]
  message("Case info: matched ", n_distinct(id_tbl$Case_ID) - length(no_meta), " of ",
          n_distinct(id_tbl$Case_ID), " cases")
  if (length(no_meta)) warning("No case information for: ", paste(no_meta, collapse = ", "), call. = FALSE)
  meta_check <- bind_rows(
    tibble(Issue = "In Luminex data but not in case info file", Case_ID = no_meta),
    tibble(Issue = "In case info file but not in Luminex data", Case_ID = unused))
} else for (v in names(CASE_INFO_COLUMNS)) id_tbl[[v]] <- NA

id_tbl <- id_tbl |>
  select(Plate, Well, Sample_ID, Case_ID, Region, all_of(names(CASE_INFO_COLUMNS)),
         !!bca_col := BCA, Protein_pgmL) |>
  arrange(Plate, Sample_ID)

group_check <- id_tbl |> filter(!is_bridge(Sample_ID)) |>
  distinct(Case_ID, Label = .data[[GROUP_COLUMN]]) |>
  mutate(Group = recode_group(Label, Case_ID),
         Excluded = id_key(Case_ID) %in% id_key(EXCLUDE_CASES),
         Override = id_key(Case_ID) %in% id_key(names(GROUP_OVERRIDES))) |>
  count(Label, Group, Override, Excluded, name = "N_Cases") |> arrange(Group, Label)
if (any(is.na(group_check$Group) & !is.na(group_check$Label)))
  warning("Some labels match none of GROUP_RULES; see the Group_Check sheet.", call. = FALSE)

# ---- Wide tables ------------------------------------------------------------------------
wide <- function(val) x |> select(Plate, Sample_ID, Analyte, v = all_of(val)) |>
  pivot_wider(names_from = Analyte, values_from = v) |>
  select(Plate, Sample_ID, all_of(analytes))
norm_num <- id_tbl |> left_join(wide("Norm"), by = c("Plate", "Sample_ID"))
raw_num  <- id_tbl |> left_join(wide("Obs"),  by = c("Plate", "Sample_ID"))
codes    <- x |> filter(!is.na(Code)) |> select(Plate, Sample_ID, Analyte, Code)
extrap   <- x |> filter(Extrapolated) |> select(Plate, Sample_ID, Analyte)

# ---- Write -------------------------------------------------------------------------------
wb <- createWorkbook()
hs <- HEADER_STYLE()
yellow <- createStyle(fgFill = "#FFF2CC")
sci <- createStyle(numFmt = "0.000E+00")
write_values <- function(sheet, df, numfmt) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df, headerStyle = hs)
  rows <- match(paste(codes$Plate, codes$Sample_ID), paste(df$Plate, df$Sample_ID)) + 1
  cols <- match(codes$Analyte, names(df))
  for (i in seq_along(rows)) writeData(wb, sheet, codes$Code[i], startRow = rows[i], startCol = cols[i])
  a_cols <- match(analytes, names(df))
  if (!is.null(numfmt)) addStyle(wb, sheet, numfmt, rows = 2:(nrow(df) + 1), cols = a_cols, gridExpand = TRUE)
  er <- match(paste(extrap$Plate, extrap$Sample_ID), paste(df$Plate, df$Sample_ID)) + 1
  ec <- match(extrap$Analyte, names(df))
  if (length(er)) addStyle(wb, sheet, yellow, rows = er, cols = ec, gridExpand = FALSE, stack = TRUE)
  freezePane(wb, sheet, firstRow = TRUE, firstActiveCol = match("Protein_pgmL", names(df)) + 1)
  setColWidths(wb, sheet, seq_len(ncol(df)), 13)
}
write_values("Normalized", norm_num, sci)
write_values("Raw_pg_per_mL", raw_num, NULL)

addWorksheet(wb, "Group_Check")
writeData(wb, "Group_Check", group_check, headerStyle = hs)
setColWidths(wb, "Group_Check", 1:5, c(45, 12, 10, 10, 10))
if (!is.null(meta_check)) {
  addWorksheet(wb, "Metadata_check")
  writeData(wb, "Metadata_check", if (nrow(meta_check)) meta_check else
    tibble(Issue = "All cases matched", Case_ID = ""), headerStyle = hs)
  setColWidths(wb, "Metadata_check", 1:2, c(45, 20))
}
notes <- tibble(Note = c(
  sprintf("Normalized = Obs Conc (pg/mL) / BCA protein (pg/mL) x %g. Unitless (pg analyte per pg protein, rescaled).", SCALE),
  sprintf("%s = BCA 'Actual Concentration' as reported (%s); Protein_pgmL = the same value in pg/mL (x %g).",
          bca_col, BCA_UNITS, to_pgmL),
  "Yellow cells: Bio-Plex extrapolated the value beyond the standard curve (shown as *value in the export). Kept as numbers.",
  "OOR < : below the curve.  OOR > : above the curve.  *** : no value.  --- : outlier removed by Bio-Plex. Not imputed.",
  if (has_regions()) paste0("Case_ID and Region come from the Sample ID (last piece after - _ or space). Regions: ",
                            paste(REGIONS, collapse = ", ")) else "Single-tissue study: Case_ID = Sample ID.",
  if (is.null(BRIDGE_PATTERN)) "No bridge samples defined." else paste0("Bridge samples: Sample IDs matching ", BRIDGE_PATTERN),
  if (is.null(CASE_INFO_FILE)) "No case information file: fill in the case columns by hand."
  else paste0("Case information from ", basename(CASE_INFO_FILE), ", matched on Case_ID."),
  "Do not add columns to the right of the analytes: everything after Protein_pgmL is read as an analyte.",
  paste0("Files: ", paste(basename(c(pairs$Luminex, pairs$BCA)), collapse = " | ")),
  paste0("Made ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by 01_normalize.R")))
addWorksheet(wb, "Notes"); writeData(wb, "Notes", notes, headerStyle = hs); setColWidths(wb, "Notes", 1, 120)
saveWorkbook(wb, FILES$norm, overwrite = TRUE)

message("\nWrote ", FILES$norm, ": ", nrow(norm_num), " samples x ", length(analytes), " analytes on ",
        n_distinct(norm_num$Plate), " plate(s)")
print(as.data.frame(group_check))
