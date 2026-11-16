############################################
# 14c_volcano.R — publication-style volcano plots from saved DEG tables (no re-testing)
# Run in R from the folder containing 14_DEG/ (after 14a and/or 14b)
# Packages: ggplot2, ggrepel
############################################

library(ggplot2)
library(ggrepel)

############################################
# Settings
############################################

fdr          <- 0.05
lfc_cut      <- 0.25                          # significant = FDR < fdr AND |log2FC| >= lfc_cut
n_labels     <- 20
label_biotypes <- c("protein_coding", "lncRNA") # only these are labelled; all genes are plotted
dpi          <- 600

methods <- list(
  pseudobulk = list(dir = "14_DEG/pseudobulk", prefix = "DEG_pseudobulk_",
                    lfc = "log2FoldChange_shrunk", name = "Pseudobulk DESeq2"),
  MAST       = list(dir = "14_DEG/MAST", prefix = "DEG_MAST_",
                    lfc = "log2FC", name = "MAST (fixed effects)")
)

cols <- c("Up in ASD" = "#E64B35", "Down in ASD" = "#4DBBD5", "NS" = "grey80")

############################################
# Plot function
############################################

# Tick positions for a square-root y axis (real -log10 p values)
sqrt_breaks <- function(lim) {
  top <- lim[2]
  cand <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300, 500)
  b <- cand[cand <= top]
  if (length(b) > 7) b <- b[c(1, seq(length(b) - 5, length(b)))]
  b
}

volcano <- function(tab, lfc_col, title) {
  tab <- tab[!is.na(tab$pvalue) & !is.na(tab$padj) & !is.na(tab[[lfc_col]]), ]
  tab$lfc <- tab[[lfc_col]]

  # p = 0 (numerical underflow): plot at the smallest non-zero p in this table
  tab$zero <- tab$pvalue == 0
  n_zero <- sum(tab$zero)
  if (n_zero > 0) tab$pvalue[tab$pvalue == 0] <- min(tab$pvalue[tab$pvalue > 0])
  tab$y <- -log10(tab$pvalue)

  tab$sig <- "NS"
  tab$sig[tab$padj < fdr & tab$lfc >=  lfc_cut] <- "Up in ASD"
  tab$sig[tab$padj < fdr & tab$lfc <= -lfc_cut] <- "Down in ASD"
  tab$sig <- factor(tab$sig, levels = names(cols))

  # FDR line: p-value of the least significant gene with padj < fdr
  fdr_y <- if (any(tab$padj < fdr)) -log10(max(tab$pvalue[tab$padj < fdr])) else NA

  # Labels: significant protein-coding/lncRNA genes with a real symbol (no Ensembl IDs)
  has_symbol <- !is.na(tab$gene_name) & tab$gene_name != "" & !grepl("^ENSG", tab$gene_name)
  lab <- tab[tab$sig != "NS" & tab$gene_biotype %in% label_biotypes & has_symbol, ]
  lab <- head(lab[order(lab$padj, lab$pvalue, -abs(lab$lfc)), ], n_labels)

  n_up <- sum(tab$sig == "Up in ASD"); n_dn <- sum(tab$sig == "Down in ASD")
  xr <- max(abs(tab$lfc)) * 1.05
  ytop <- max(tab$y) * 1.45                     # headroom so labels can sit above the top genes

  tab <- tab[order(tab$sig != "NS"), ]          # significant points drawn on top

  p <- ggplot(tab, aes(lfc, y)) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    { if (!is.na(fdr_y)) geom_hline(yintercept = fdr_y, linetype = "dashed", linewidth = 0.3, color = "grey40") } +
    geom_point(data = tab[!tab$zero, ], aes(color = sig), size = 0.7, alpha = 0.8, stroke = 0) +
    geom_point(data = tab[tab$zero, ], aes(color = sig), size = 1.4, shape = 17) +
    geom_text_repel(data = lab, aes(label = gene_name), size = 2.3, fontface = "italic",
                    max.overlaps = Inf, segment.size = 0.2, segment.color = "grey50",
                    min.segment.length = 0, box.padding = 0.4, point.padding = 0.15,
                    force = 3, force_pull = 0.5, max.time = 3, max.iter = 1e5, seed = 42) +
    scale_color_manual(values = cols, drop = FALSE) +
    scale_x_continuous(limits = c(-xr, xr)) +
    scale_y_continuous(trans = "sqrt", limits = c(0, ytop), breaks = sqrt_breaks,
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, subtitle = paste0("Down: ", n_dn, "   |   Up: ", n_up),
         x = expression(log[2]~"fold change (ASD vs Control)"),
         y = expression(-log[10]~italic(p)~"(square-root scale)"), color = NULL,
         caption = paste0("Significant: FDR < ", fdr, " and |log2FC| >= ", lfc_cut,
                          "; labels: top ", n_labels, " protein-coding/lncRNA genes with a symbol",
                          if (n_zero > 0) paste0("\nTriangles: ", n_zero, " genes with p below the numerical limit (reported as 0), plotted at the smallest non-zero p") else "")) +
    theme_classic(base_size = 8) +
    theme(legend.position = "top", legend.key.size = unit(3, "mm"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 7),
          plot.caption = element_text(size = 5.5, color = "grey40"),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3)) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1, shape = 16)))

  list(plot = p, n_up = n_up, n_down = n_dn, n_zero = n_zero, max_y = max(tab$y))
}

############################################
# Run over all saved tables
############################################

for (m in names(methods)) {
  cfg <- methods[[m]]
  files <- list.files(cfg$dir, pattern = paste0("^", cfg$prefix, ".*\\.csv$"), full.names = TRUE)
  if (length(files) == 0) { cat("No tables found for", m, "\n"); next }
  vdir <- file.path(cfg$dir, "volcano")
  dir.create(vdir, showWarnings = FALSE)

  rows <- list()
  pdf(file.path(vdir, paste0("volcano_all_", m, ".pdf")), width = 4.5, height = 4.5)
  for (f in files) {
    st <- sub("\\.csv$", "", sub(cfg$prefix, "", basename(f)))
    tab <- read.csv(f, check.names = FALSE)
    v <- volcano(tab, cfg$lfc, paste0(st, " - ", cfg$name))
    print(v$plot)
    safe <- gsub("[^A-Za-z0-9_-]", "", st)          # "PVALB+" -> "PVALB" in file names
    ggsave(file.path(vdir, paste0("volcano_", m, "_", safe, ".pdf")), v$plot, width = 4.5, height = 4.5)
    ggsave(file.path(vdir, paste0("volcano_", m, "_", safe, ".png")), v$plot, width = 4.5, height = 4.5, dpi = dpi)
    rows[[st]] <- data.frame(subtype = st, up = v$n_up, down = v$n_down,
                             p_zero = v$n_zero, max_neglog10p = round(v$max_y, 1))
  }
  dev.off()
  s <- do.call(rbind, rows)
  write.csv(s, file.path(vdir, paste0("volcano_counts_", m, ".csv")), row.names = FALSE)
  cat("\n", cfg$name, " (FDR < ", fdr, ", |log2FC| >= ", lfc_cut, "):\n", sep = "")
  print(s, row.names = FALSE)
}
