# =============================================================================
# STEP 5  Figures for the step-4 lmPerm results
#
# Input : output/Luminex_lmPerm_Results.xlsx
# Output: output/Luminex_lmPerm_Pairwise_A.pdf         Model A pairwise comparisons
#         output/Luminex_lmPerm_Pairwise_B.pdf         Model B pairwise comparisons
#         output/Luminex_lmPerm_Slopes_B.pdf           Model B slope comparisons
#         output/Luminex_lmPerm_Effect_A_<covariate>.pdf   numeric covariate terms (Model A)
#
# A panel (one analyte in one region) is drawn when any test in it has a raw
# permutation p below P_CUTOFF. Labels read  p = 0.0032 (FDR = 0.043), using the
# BH FDR from the results workbook. Stars follow the raw p:
#   * p < 0.05   ** p < 0.01   *** p < 0.001
# Perm_P = 0 means no permuted statistic beat the observed one: shown as p < 1/N_PERM.
#
# Boxes and points show the observed log2 values. p-values come from the adjusted
# model; the lines on the covariate plots are simple fits drawn to show direction.
# =============================================================================

if (!file.exists("R/pipeline_helpers.R"))
  stop("Run this from the Luminex_Pipeline folder: open Luminex_Pipeline.Rproj in RStudio, ",
       "or setwd(\"path/to/Luminex_Pipeline\") first.", call. = FALSE)
source("R/pipeline_helpers.R")
need_file(FILES$lmperm, "04_lmperm_models.R")
RES <- FILES$lmperm
sheets   <- readxl::excel_sheets(RES)
settings <- read_excel(RES, sheet = "Settings")
get_set  <- function(k) settings$Value[settings$Setting == k][1]

n_perm <- suppressWarnings(as.numeric(get_set("Iterations per term")))
if (is.na(n_perm)) n_perm <- N_PERM
lmperm_ver <- get_set("lmPerm version")
if (!is.na(lmperm_ver) && lmperm_ver != "2.1.0")
  warning("These results were made with lmPerm ", lmperm_ver, ". Models with a factor of 3+ levels ",
          "need 2.1.0; p-values from other versions are not reliable.", call. = FALSE, immediate. = TRUE)
model_label <- c(A = get_set("Model A"), B = get_set("Model B"))
iv <- get_set("Interaction covariate"); if (is.na(iv) || iv == "none") iv <- NULL

# Group order and colours as recorded with the results, so figures match the tables
lev     <- str_split(get_set("Groups (reference first)"), ",\\s*")[[1]]
order_x <- str_split(get_set("Group plot order"), ",\\s*")[[1]]
cols    <- group_colors(lev)

dat <- read_excel(RES, sheet = "Data_Used") |>
  mutate(Group = factor(Group, levels = order_x), Region = as.character(Region))
regions_in <- unique(dat$Region)
single_region <- length(regions_in) == 1

# ---- Label helpers -------------------------------------------------------------------
floor_txt <- format(1 / n_perm, scientific = FALSE)
fmt_p   <- function(p) ifelse(p == 0, paste("<", floor_txt), sprintf("= %.4f", p))
fmt_fdr <- function(q) ifelse(q == 0, paste("<", floor_txt), ifelse(q < 0.001, sprintf("= %.4f", q), sprintf("= %.3f", q)))
p_label <- function(p, q) paste0("p ", fmt_p(p), " (FDR ", fmt_fdr(q), ")")
stars   <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
panel_name <- function(a, r) if (single_region) a else paste0(a, "  |  ", r)

# Full width for a full page, left part of the page for a short last page
print_pages <- function(plot_fun, panels, file, width = 11, height = 8.5) {
  if (!length(panels)) { message("  nothing below p < ", P_CUTOFF, " for ", basename(file)); return(invisible()) }
  pdf(file, width = width, height = height)
  for (ch in chunks(panels, PER_PAGE)) {
    p <- plot_fun(ch)
    if (length(ch) < PER_PAGE) {
      grid::grid.newpage()
      print(p, vp = grid::viewport(x = 0.5 * length(ch) / PER_PAGE, width = length(ch) / PER_PAGE))
    } else print(p)
  }
  invisible(dev.off())
  message("  wrote ", file, " (", length(panels), " panels)")
}

base_theme <- theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12), plot.title = element_text(size = 13),
        plot.subtitle = element_text(size = 10), panel.grid.minor = element_blank())
y_lab <- expression(log[2]~"(normalized value)")

# ---- 1. Pairwise box plots ---------------------------------------------------------------
plot_pairwise <- function(model) {
  sh <- paste0("Pairwise_", model)
  if (!sh %in% sheets) return(invisible())
  pw <- read_excel(RES, sheet = sh) |> mutate(Region = as.character(Region))
  hits <- pw |> filter(Perm_P < P_CUTOFF)
  if (!nrow(hits)) { message("  nothing below p < ", P_CUTOFF, " for Model ", model); return(invisible()) }
  # bracket span on the plot's x order; short brackets sit lowest
  g <- str_split_fixed(hits$Contrast, " vs ", 2)
  hits <- hits |> mutate(x1 = pmin(match(g[, 1], order_x), match(g[, 2], order_x)),
                         x2 = pmax(match(g[, 1], order_x), match(g[, 2], order_x)))
  order_tbl <- hits |> group_by(Region, Analyte) |> summarise(minP = min(Perm_P), .groups = "drop") |>
    arrange(minP, Analyte, Region) |> mutate(Panel = panel_name(Analyte, Region))

  draw <- function(ch) {
    keys <- order_tbl |> filter(Panel %in% ch)
    d <- dat |> inner_join(keys |> select(Region, Analyte, Panel), by = c("Region", "Analyte"))
    rng <- d |> group_by(Panel) |> summarise(top = max(y), span = diff(range(y)), .groups = "drop")
    br <- hits |> mutate(Panel = panel_name(Analyte, Region)) |> filter(Panel %in% ch) |>
      arrange(Panel, x2 - x1, x1) |> group_by(Panel) |> mutate(level = row_number()) |> ungroup() |>
      left_join(rng, by = "Panel") |>
      mutate(yb = top + span * (0.10 + 0.19 * (level - 1)), tick = span * 0.03,
             star = stars(Perm_P), lab = p_label(Perm_P, FDR))
    room <- br |> group_by(Panel) |> summarise(y = max(yb) + span[1] * 0.14, .groups = "drop")
    fix <- function(x) mutate(x, Panel = factor(Panel, levels = keys$Panel))
    set.seed(1)
    ggplot(fix(d), aes(Group, y)) +
      geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.85, width = 0.6) +
      geom_jitter(width = 0.15, size = 1.3, alpha = 0.6) +
      geom_segment(data = fix(br), aes(x = x1, xend = x2, y = yb, yend = yb), inherit.aes = FALSE) +
      geom_segment(data = fix(br), aes(x = x1, xend = x1, y = yb, yend = yb - tick), inherit.aes = FALSE) +
      geom_segment(data = fix(br), aes(x = x2, xend = x2, y = yb, yend = yb - tick), inherit.aes = FALSE) +
      geom_text(data = fix(br), aes(x = (x1 + x2) / 2, y = yb + span * 0.015, label = star),
                vjust = 0, size = 5, inherit.aes = FALSE) +
      geom_text(data = fix(br), aes(x = (x1 + x2) / 2, y = yb + span * 0.075, label = lab),
                vjust = 0, size = 3.4, inherit.aes = FALSE) +
      geom_blank(data = fix(room), aes(y = y), inherit.aes = FALSE) +
      scale_x_discrete(drop = FALSE) +
      scale_fill_manual(values = cols, guide = "none", drop = FALSE) +
      facet_wrap(~ Panel, ncol = PER_PAGE, scales = "free_y") +
      labs(title = paste0("Model ", model, ":  ", model_label[[model]]),
           subtitle = paste0("Pairwise comparisons with permutation p < ", P_CUTOFF,
                             ". Label: p (BH FDR within region x comparison). * p < 0.05  ** p < 0.01  *** p < 0.001",
                             if (model == "B" && !is.null(iv)) paste0("\nComparisons are at the region's mean ", iv, ".") else ""),
           x = NULL, y = y_lab) +
      base_theme
  }
  print_pages(draw, order_tbl$Panel, out_path(paste0("Luminex_lmPerm_Pairwise_", model, ".pdf")))
}

# ---- 2. Slope comparisons (Model B) ----------------------------------------------------------
plot_slopes <- function() {
  if (!"Slopes_B" %in% sheets || is.null(iv)) return(invisible())
  sl <- read_excel(RES, sheet = "Slopes_B") |> mutate(Region = as.character(Region))
  hits <- sl |> filter(Perm_P < P_CUTOFF)
  if (!nrow(hits)) { message("  nothing below p < ", P_CUTOFF, " for Model B slopes"); return(invisible()) }
  order_tbl <- hits |> group_by(Region, Analyte) |>
    summarise(minP = min(Perm_P), txt = paste0(Contrast, " slope: ", p_label(Perm_P, FDR), " ", stars(Perm_P),
                                               collapse = "\n"), .groups = "drop") |>
    arrange(minP, Analyte, Region) |> mutate(Panel = panel_name(Analyte, Region))
  draw <- function(ch) {
    keys <- order_tbl |> filter(Panel %in% ch)
    d <- dat |> inner_join(keys |> select(Region, Analyte, Panel), by = c("Region", "Analyte")) |>
      mutate(xv = .data[[iv]])
    rng <- d |> group_by(Panel) |> summarise(top = max(y), span = diff(range(y)), .groups = "drop")
    ann <- keys |> left_join(rng, by = "Panel") |>
      mutate(y_room = top + span * (0.10 + 0.09 * (str_count(txt, "\n") + 1)))
    fix <- function(x) mutate(x, Panel = factor(Panel, levels = keys$Panel))
    ggplot(fix(d), aes(xv, y, colour = Group)) +
      geom_point(size = 1.6, alpha = 0.7) +
      geom_smooth(aes(fill = Group), method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15, linewidth = 1) +
      geom_text(data = fix(ann), aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
                hjust = -0.03, vjust = 1.3, size = 3.4, lineheight = 1.05) +
      geom_blank(data = fix(ann), aes(x = NULL, y = y_room), inherit.aes = FALSE) +
      scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
      facet_wrap(~ Panel, ncol = PER_PAGE, scales = "free_y") +
      labs(title = paste0("Model B:  ", model_label[["B"]]),
           subtitle = paste0(iv, ":Group slope comparisons with permutation p < ", P_CUTOFF,
                             ". Label: p (BH FDR within region x comparison).\n",
                             "Lines are simple per-group fits with 95% bands; p-values are from the adjusted model."),
           x = iv, y = y_lab) +
      base_theme + theme(legend.position = "bottom")
  }
  print_pages(draw, order_tbl$Panel, out_path("Luminex_lmPerm_Slopes_B.pdf"))
}

# ---- 3. Numeric covariate terms (Model A) ------------------------------------------------------
plot_covariate <- function(cov) {
  tm <- read_excel(RES, sheet = "Terms_A") |> mutate(Region = as.character(Region))
  order_tbl <- tm |> filter(Term == cov, Perm_P < P_CUTOFF) |> arrange(Perm_P, Analyte, Region) |>
    mutate(Panel = panel_name(Analyte, Region), txt = paste0(cov, ": ", p_label(Perm_P, FDR), " ", stars(Perm_P)))
  draw <- function(ch) {
    keys <- order_tbl |> filter(Panel %in% ch)
    d <- dat |> inner_join(keys |> select(Region, Analyte, Panel), by = c("Region", "Analyte")) |>
      mutate(xv = .data[[cov]])
    rng <- d |> group_by(Panel) |> summarise(top = max(y), span = diff(range(y)), .groups = "drop")
    ann <- keys |> left_join(rng, by = "Panel") |> mutate(y_room = top + span * 0.18)
    fix <- function(x) mutate(x, Panel = factor(Panel, levels = keys$Panel))
    ggplot(fix(d), aes(xv, y)) +
      geom_point(aes(colour = Group), size = 1.6, alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20", fill = "grey80", linewidth = 1) +
      geom_text(data = fix(ann), aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
                hjust = -0.03, vjust = 1.4, size = 3.4) +
      geom_blank(data = fix(ann), aes(y = y_room), inherit.aes = FALSE) +
      scale_colour_manual(values = cols, name = NULL) +
      facet_wrap(~ Panel, ncol = PER_PAGE, scales = "free_y") +
      labs(title = paste0("Model A:  ", model_label[["A"]]),
           subtitle = paste0(cov, " term with permutation p < ", P_CUTOFF, ". Label: p (BH FDR within region x term).\n",
                             "Line is a simple fit across all cases with a 95% band; p-values are from the adjusted model."),
           x = cov, y = y_lab) +
      base_theme + theme(legend.position = "bottom")
  }
  print_pages(draw, order_tbl$Panel, out_path(paste0("Luminex_lmPerm_Effect_A_", cov, ".pdf")))
}

# ---- Run -----------------------------------------------------------------------------------------
message("Results: ", RES, "  (lmPerm ", lmperm_ver, ", ", n_perm, " iterations)")
plot_pairwise("A")
plot_pairwise("B")
plot_slopes()
num_covs <- names(dat)[vapply(dat, is.numeric, logical(1))]
for (cv in intersect(COVARIATES, num_covs)) plot_covariate(cv)
