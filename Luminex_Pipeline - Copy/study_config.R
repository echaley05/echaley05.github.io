# =============================================================================
# STUDY SETTINGS: the only file you need to edit for a new experiment
#
# Every script in the pipeline reads this file. Change the values below to
# describe your study, save, then run run_pipeline.R (or the steps one by one).
#
# The values filled in here are for the Cross-Disorder Luminex study
# (CTRL / ASD / SMI, three brain regions, six plates). They double as a worked
# example; see README.md for a 2-group single-tissue example.
# =============================================================================

# ---- 1. Where your files are ---------------------------------------------------
DATA_DIR   <- "data"      # Bio-Plex exports, BCA files and the case info file
OUTPUT_DIR <- "output"    # every workbook and PDF is written here

# File names only need to contain "plex" (Bio-Plex export) or "bca" (BCA file),
# plus "Plate" and a number, e.g. "...48 plex ... Plate 3.xlsx" and
# "9-16-26 Final BCA data Plate 3.xlsx". Files are paired by that plate number.
LUMINEX_PATTERN <- "(?i)plex.*plate[ _-]?\\d+.*\\.xlsx?$"
BCA_PATTERN     <- "(?i)bca.*plate[ _-]?\\d+.*\\.xlsx?$"

# ---- 2. Normalization -------------------------------------------------------------
# Normalized = Obs Conc (pg/mL) / BCA protein (pg/mL) x SCALE
BCA_UNITS <- "ug/mL"      # units of the BCA "Actual Concentration" column: "ug/mL" or "mg/mL"
SCALE     <- 1e10          # multiplier so the unitless ratio is readable

# ---- 3. Sample IDs: case and region ----------------------------------------------
# Bio-Plex Sample IDs are read as  <Case ID><separator><Region>,  e.g. "PGEQ-STG".
# The region is the last piece after "-", "_" or a space.
# If every sample is the same tissue, set REGIONS <- NULL; the whole Sample ID
# is then the Case ID.
REGIONS <- c("STG", "ERC", "ACC")   # order is used in tables and plots

# Bridge / inter-plate control samples. They are used in QC and never analysed.
# A regular expression matched against the Sample ID; NULL if you ran none.
BRIDGE_PATTERN <- "-(FC|OC)$"

# Cases to leave out of every analysis (e.g. bridging-only cases).
EXCLUDE_CASES <- c("BEMH-0016", "BEMH-0046")

# ---- 4. Case information file --------------------------------------------------------
# One row per case. CSV or Excel. Set CASE_INFO_FILE <- NULL to skip; the
# normalized workbook then has blank columns to fill in by hand.
CASE_INFO_FILE  <- "data/Case_Info-Grid_view.csv"
CASE_INFO_SHEET <- 1                 # sheet name or number (Excel only)
CASE_ID_COLUMN  <- "Case ID"         # column in YOUR file that holds the case ID

# Columns to bring in:  name in the output = name in your file.
# Output names must be simple (letters, numbers, underscores; no spaces).
CASE_INFO_COLUMNS <- c(Diagnosis = "Diagnosis",
                       Age       = "Age",
                       Sex       = "Sex",
                       PMI       = "PMI")

# ---- 5. Groups -------------------------------------------------------------------------
# Which output column above holds the group (diagnosis, genotype, treatment...).
GROUP_COLUMN <- "Diagnosis"

# How labels in that column become analysis groups. Rules are checked top to
# bottom and the first match wins (regular expressions; "(?i)" = ignore case).
# "." matches anything, so a last rule of "." catches every remaining label.
# Set GROUP_RULES <- NULL to use the labels exactly as they are written.
GROUP_RULES <- c(ASD  = "(?i)asd|autism",
                 CTRL = "(?i)^control\\*?$",
                 SMI  = ".")

# Cases whose group must be set by hand, whatever their label says.
# BEMH-0041 is labelled "Control (Depression with family history of Bipolar)".
GROUP_OVERRIDES <- c("BEMH-0041" = "SMI")

# Statistical order: the FIRST group is the reference. Pairwise comparisons are
# reported as "later vs earlier", e.g. "ASD vs CTRL". NULL = alphabetical.
GROUP_LEVELS <- c("CTRL", "ASD", "SMI")

# Left-to-right order on plots. NULL = same as GROUP_LEVELS.
GROUP_PLOT_ORDER <- c("ASD", "CTRL", "SMI")

# Plot colours. NULL = automatic colour-blind-safe palette.
GROUP_COLORS <- c(CTRL = "#56B4E9", ASD = "#E69F00", SMI = "#009E73")

# ---- 6. Statistical model ----------------------------------------------------------------
# Fitted separately in each region, one row per case:
#   Model A:  log2(value) ~ Group + COVARIATES [+ Plate]
#   Model B:  Model A + Group:INTERACTION_COVARIATE
# Covariates must be output names from section 4. Numeric columns are treated
# as numbers, text columns (Sex) as categories. A covariate with only one value
# in a region (e.g. all male) is dropped from that region's model automatically.
COVARIATES           <- c("Age", "PMI", "Sex")
ADJUST_FOR_PLATE     <- TRUE        # Plate as a fixed factor (ignored with one plate)
INTERACTION_COVARIATE <- "Age"      # Model B; NULL to skip Model B

MIN_DETECT <- 40          # % of samples in range for an analyte to be analysed
MIN_N      <- 20          # fewest cases for a region model to be fitted
MIN_GROUP  <- 5           # groups with fewer cases in a region are left out there
FDR_ALPHA  <- 0.05        # Benjamini-Hochberg threshold

# ---- 7. Permutation settings (lmPerm) ------------------------------------------------------
N_PERM       <- 100000     # permutation iterations per term
SEED         <- 20261005  # makes results reproducible
RUN_SELFTEST <- TRUE      # checks lmPerm gives ~5% false positives on your design
AUTO_INSTALL_LMPERM <- TRUE   # install lmPerm 2.1.0 from the CRAN archive if missing

# ---- 8. QC settings ----------------------------------------------------------------------------
RECOVERY_MIN <- 70        # % recovery window for a standard to count
RECOVERY_MAX <- 130
BEAD_MIN     <- 35        # beads per region below which a well is unreliable
DETECT_TIERS <- c(0.20, 0.80)

# ---- 9. Plots ------------------------------------------------------------------------------------
P_CUTOFF <- 0.05          # raw permutation p below which a panel is drawn
PER_PAGE <- 2             # panels per page

# ---- 10. Multivariate (optional step 6) -----------------------------------------------------
MV_DETECT    <- 65        # stricter detection floor for the multivariate matrix
IMPUTE_MV    <- TRUE      # below-range -> LLOQ/2, above-range -> ULOQ (multivariate only)
IMPUTE_LOW   <- "half_lloq"   # or "lloq_sqrt2"
IMPUTE_HIGH  <- TRUE
MV_N_PERM    <- 999       # permutations for PERMANOVA and the PLS-DA null
MAX_COMP     <- 4         # PLS-DA components to consider
