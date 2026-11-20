# Luminex analysis pipeline

A standard pipeline for Bio-Plex (Luminex) multiplex data with BCA protein
normalization. It works for any number of groups, one or more tissue regions,
one or more plates, and whatever covariates you have. Everything specific to
your experiment goes in **one file: `study_config.R`**.

```
Luminex_Pipeline/
├── Luminex_Pipeline.Rproj     open this in RStudio (sets the working folder)
├── study_config.R             ← the only file you edit
├── run_pipeline.R             runs every step in order
├── data/                      your Bio-Plex exports, BCA files, case info file
├── output/                    every workbook and PDF lands here
├── templates/case_info_template.csv
├── scripts/
│   ├── 01_normalize.R         BCA normalization + case info
│   ├── 02_qc.R                standard curves, detection, bridges, bead counts
│   ├── 03_distribution_check.R  is log2 right? box plots
│   ├── 04_lmperm_models.R     main statistics (permutation linear models)
│   ├── 05_plot_lmperm.R       figures with significance brackets
│   └── 06_multivariate.R      optional: modules, PERMANOVA, PLS-DA
└── R/pipeline_helpers.R       shared code (don't edit)
```

---

## 1. One-time setup

1. Install R (4.2 or newer) and RStudio.
2. Install the packages (paste into the R console):
   ```r
   install.packages(c("readxl", "dplyr", "tidyr", "stringr", "purrr", "openxlsx",
                      "ggplot2", "MASS", "vegan", "pls", "cluster"))
   ```
   If `openxlsx` is older than 4.2.6, update it; older versions write files Excel asks to "repair".
3. **lmPerm must be version 2.1.0.** The current CRAN version (2.1.6) gives badly
   wrong p-values whenever a factor has 3 or more levels (3+ groups, or 3+ plates).
   Step 4 installs 2.1.0 automatically. To do it yourself:
   ```r
   install.packages("https://cran.r-project.org/src/contrib/Archive/lmPerm/lmPerm_2.1.0.tar.gz",
                    repos = NULL, type = "source")
   ```
   This compiles from source, so Windows needs **Rtools** (cran.r-project.org/bin/windows/Rtools)
   and a Mac needs the Xcode command line tools (`xcode-select --install` in Terminal).
   Restart R afterwards and check with `packageVersion("lmPerm")`.

## 2. Each new experiment

1. **Copy the whole `Luminex_Pipeline` folder** and rename it for your experiment.
   Each experiment gets its own copy, so settings and results never mix.
2. **Put your files in `data/`:**
   - one Bio-Plex Manager Excel export per plate,
   - one BCA file per plate,
   - one case information file (CSV or Excel).
3. **Edit `study_config.R`** (section 4 below).
4. Double-click `Luminex_Pipeline.Rproj`, then in the console:
   ```r
   source("run_pipeline.R")
   ```
   Or run one step at a time, e.g. `source("scripts/01_normalize.R")`.
5. **Check the checkpoints (section 5) before you look at any p-value.**

## 3. What the input files must look like

**File names.** Each file needs "plex" or "bca" in its name, plus `Plate` and a number:
`... 48 plex brain ext Plate 3.xlsx` and `9-16-26 Final BCA data Plate 3.xlsx`.
Spaces, `_` or `-` before the number all work. Files are paired by plate number,
so the dates in the names don't matter.

**Bio-Plex export.** The standard Bio-Plex Manager Excel report. Step 1 needs the
`Obs Conc` sheet. Step 2 also uses `FI`, `Exp Conc`, `Conc in Range`, `Bead Count` and
`Standard Curve`. Analyte names must keep their bead region, e.g. `IL-6 (25)`.

**BCA file.** Any sheet with a header row that contains `Sample ID` and a column
whose name starts with `Actual` (the Actual Concentration). Sample IDs must match
the Bio-Plex Sample IDs.

**Sample IDs.** Written as `<Case ID>-<Region>`, e.g. `PGEQ-STG`, `P-014_HIP`.
The region is the last piece after `-`, `_` or a space. Single-tissue study: use the
case ID alone and set `REGIONS <- NULL`.

**Case information file.** One row per case: an ID column plus whatever you have
(group, age, sex, PMI, ...). Column names can be anything; you map them in the
config. Copy `templates/case_info_template.csv` to start one. Case IDs match even if
spacing, hyphens or capitals differ (`BEMH-0016` = `bemh 0016`).

## 4. Editing `study_config.R`

The file is commented line by line. The settings that change between experiments:

| Setting | What it means |
|---|---|
| `REGIONS` | Region codes at the end of Sample IDs, or `NULL` for one tissue |
| `BRIDGE_PATTERN` | Pattern that marks bridge / inter-plate samples, or `NULL` |
| `EXCLUDE_CASES` | Cases to leave out of all statistics |
| `CASE_INFO_FILE`, `CASE_ID_COLUMN` | Your case file and its ID column |
| `CASE_INFO_COLUMNS` | Columns to bring in: `OutputName = "Name in your file"` |
| `GROUP_COLUMN` | Which of those columns is the group |
| `GROUP_RULES` | How labels collapse into groups (`NULL` = use labels as written) |
| `GROUP_LEVELS` | Group order; **the first is the reference** |
| `COVARIATES` | Covariates in the model (numbers or categories) |
| `INTERACTION_COVARIATE` | Model B's slope comparison (usually `"Age"`), or `NULL` |

**Example: a 2-group, single-tissue, single-plate study:**
```r
REGIONS           <- NULL
BRIDGE_PATTERN    <- NULL
EXCLUDE_CASES     <- NULL
CASE_INFO_FILE    <- "data/sample_info.csv"
CASE_ID_COLUMN    <- "Subject"
CASE_INFO_COLUMNS <- c(Status = "Status", Age = "Age", PMI = "PMI")
GROUP_COLUMN      <- "Status"
GROUP_RULES       <- NULL                    # labels already say Control / Case
GROUP_OVERRIDES   <- NULL
GROUP_LEVELS      <- c("Control", "Case")    # Control is the reference
GROUP_PLOT_ORDER  <- NULL
GROUP_COLORS      <- NULL                    # automatic colour-blind-safe colours
COVARIATES        <- c("Age", "PMI")
```

**Group rules** are checked top to bottom, and the first match wins. Example: collapsing
many diagnoses into three groups:
```r
GROUP_RULES <- c(ASD  = "(?i)asd|autism",       # anything mentioning ASD or autism
                 CTRL = "(?i)^control\\*?$",    # exactly "Control" or "Control*"
                 SMI  = ".")                     # everything else
```
`(?i)` means "ignore capitals", `^` and `$` mean start and end, and `.` matches anything.
After step 1, the **Group_Check** sheet shows which label went to which group.

## 5. Checkpoints: look at these before trusting results

| After | Open | Check |
|---|---|---|
| Step 1 | `Luminex_normalized.xlsx` → **Metadata_check** | every case found its case info |
| Step 1 | → **Group_Check** | every label went to the right group; no blank group |
| Step 2 | `Luminex_QC.xlsx` → **Std_Curves** | flagged curves (fewer than 4 good standards, poor fit) |
| Step 2 | → **Detection** | which analytes reach `MIN_DETECT` (only those are tested) |
| Step 2 | → **Plate_Effects**, QC plots | bridges: is plate drift small next to biology? |
| Step 3 | `Distribution_Check.xlsx` → **Summary_Residuals** | log2 should be at or near the top |
| Step 4 | console / **Settings** sheet | the self-test passed; the design table has enough cases per group |

## 6. What the statistics do (plain version)

- **Normalization.** Each analyte concentration is divided by that sample's total
  protein (BCA), so a sample that simply had more tissue doesn't look "higher".
  The ratio is multiplied by 10⁹ so the numbers are readable.
- **log2.** Cytokine levels are skewed: a few samples are far higher than the rest.
  log2 makes them roughly symmetric, and differences read as fold changes:
  +1 = double, −1 = half.
- **The model (step 4).** In each region, every analyte is fitted as
  `log2(value) ~ Group + covariates + Plate`. The group comparison is therefore made
  *after* allowing for age, PMI, sex and plate.
- **Permutation p-values.** Instead of assuming a bell curve, lmPerm shuffles the data
  thousands of times to see how often a difference this large turns up by chance.
  `Perm_P = 0.003` means about 3 in 1,000 shuffles did as well.
- **Pairwise comparisons.** Every pair of groups is compared once and labelled
  "later vs earlier" in `GROUP_LEVELS`, e.g. `ASD vs CTRL`. `Log2FC` > 0 means the
  first-named group is higher.
- **FDR.** Testing many analytes produces some false hits by chance. The FDR
  (Benjamini–Hochberg) column corrects for this within each region × comparison.
  Treat **FDR < 0.05** as the finding and raw p < 0.05 as "worth a look".
- **Model B** adds `Age:Group`: does the effect of age differ between groups? Its group
  comparisons are made at the region's average age.
- **Small groups.** A group with fewer than `MIN_GROUP` cases in a region is left out of
  that region's model (listed in the `Groups_Dropped` column). A covariate with only one
  value in a region (e.g. all male) is dropped there (see the `Covariates` column).

## 7. Reproducibility

Permutation p-values come from random shuffles, but the shuffles are seeded, so
**the same data and settings give identical p-values on every run**. Each fit gets
its own seed, built from `SEED` plus the model, region and analyte, and rows are sorted
by case before shuffling. Changing `MIN_DETECT`, adding an analyte or rerunning one
region therefore leaves every other result unchanged.

Results will differ when:
- **the data change**, e.g. a case gains or loses a covariate value and enters or
  leaves the model (step 4 prints a note when cases are dropped);
- **`SEED` or `N_PERM` change.** A different seed shifts p-values slightly
  (about ±0.002 near p = 0.05 at 10,000 iterations); that is normal Monte Carlo noise;
- **the R or lmPerm version differs.** Both are recorded on the Settings sheet.

`SCALE` multiplies every value by the same number, so it changes the printed values
but never a statistic.

## 8. Outputs

| File | What's in it |
|---|---|
| `Luminex_normalized.xlsx` | normalized values (yellow = extrapolated), raw values, group check |
| `Luminex_QC.xlsx` + `_plots.pdf` | curve fits, LLOQ/ULOQ, detection, bridges, bead counts, outliers |
| `Distribution_Check.xlsx`, `Distribution_QQ.pdf`, `Analyte_Boxplots.pdf` | transform check and box plots |
| `Luminex_lmPerm_Results.xlsx` | **main results**: Summary, Terms, Pairwise, Slopes, Data_Used, Settings |
| `Luminex_lmPerm_Pairwise_A/B.pdf` | box plots with brackets for p < 0.05 |
| `Luminex_lmPerm_Slopes_B.pdf`, `Luminex_lmPerm_Effect_A_<covariate>.pdf` | slope and covariate figures |
| `Luminex_Multivariate.xlsx` + `_plots.pdf` | modules, PERMANOVA, PLS-DA, region contrasts |

## 9. Troubleshooting

| Message | Fix |
|---|---|
| "Run this from the Luminex_Pipeline folder" | open `Luminex_Pipeline.Rproj`, or `setwd()` to the folder |
| "No Bio-Plex files found" | check file names contain "plex" and "Plate N"; the message lists what it found |
| "No BCA file for plate(s) N" | missing BCA file, or its name lacks "BCA" / "Plate N" |
| "Column(s) not found in <case file>" | the right-hand names in `CASE_INFO_COLUMNS` must match your headers exactly |
| "Group(s) found in the data but not in GROUP_LEVELS" | add the group, or fix `GROUP_RULES` |
| "lmPerm failed the calibration self-test" | wrong lmPerm version; reinstall 2.1.0 and restart R |
| "has no 'Type \| Well' header row" | not a Bio-Plex Manager report export |
| Step 4 is slow | expected: 10,000 permutations × every analyte × region × comparison. Set `N_PERM <- 1000` for a trial run |
