# =============================================================================
# Run the whole Luminex pipeline, in order.
#
# 1. Put your Bio-Plex exports, BCA files and case information file in data/
# 2. Edit study_config.R
# 3. Open Luminex_Pipeline.Rproj in RStudio (or setwd() to this folder) and run:
#        source("run_pipeline.R")
#
# Each step can also be run on its own, e.g. source("scripts/02_qc.R").
# Check the outputs of steps 1 and 2 (Group_Check, Metadata_check, Detection)
# before trusting anything from step 4 on.
# =============================================================================

STEPS <- c(
  "scripts/01_normalize.R",            # BCA normalization + case info  -> Luminex_normalized.xlsx
  "scripts/02_qc.R",                   # curves, detection, bridges     -> Luminex_QC.xlsx
  "scripts/03_distribution_check.R",   # transform check, box plots     -> Distribution_Check.xlsx
  "scripts/04_lmperm_models.R",        # main statistics (slow)         -> Luminex_lmPerm_Results.xlsx
  "scripts/05_plot_lmperm.R",          # significance figures           -> Luminex_lmPerm_*.pdf
  "scripts/06_multivariate.R"          # optional: modules, PERMANOVA, PLS-DA
)

for (s in STEPS) {
  message("\n==================== ", s, " ====================")
  t0 <- Sys.time()
  source(s, local = new.env())
  message(sprintf("---- done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
message("\nAll steps finished. Results are in ", OUTPUT_DIR, "/")
