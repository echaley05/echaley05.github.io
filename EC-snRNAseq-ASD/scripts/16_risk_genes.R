############################################
# 16_risk_genes.R - enrichment of ASD risk genes among MAST DEGs, per subtype
# Gene sets: SFARI (score 1; score 1-2; all scored 1-3; syndromic), Satterstrom 2020, Fu 2022
# Tests: one-sided Fisher's exact (protein-coding background) + expression-matched permutation
# Run in R from the folder containing 14_DEG/ and risk_genes/
# Packages: ggplot2
#
# Input files (put in risk_genes/):
#   SFARI_Genes.csv          - SFARI Gene Human Gene module download (any file name starting "SFARI")
#   Satterstrom_102_ASD.csv  - Satterstrom et al. 2020 102 ASD genes (column ensembl_gene_id; footnote rows ignored)
#   Fu_158_ASD.csv           - Fu et al. 2022 gene table; set = rows with ASD185 == TRUE (185 genes)
# Genes are matched by Ensembl ID, with gene symbol as a fallback. Any missing file is skipped with a message.
############################################

library(ggplot2)

############################################
# Settings
############################################

out        <- "16_risk_genes"
deg_dir    <- "14_DEG/MAST"
deg_prefix <- "DEG_MAST_"
gs_dir     <- "risk_genes"
fdr        <- 0.05
lfc_cut    <- 0.25
n_perm     <- 10000
n_bins     <- 20          # detection-rate bins for the expression-matched null
seed       <- 42
dir.create(out, showWarnings = FALSE)

############################################
# 1. Gene sets
############################################

sets <- list()       # each: list(ensembl = ..., symbol = ...)

sfari_file <- list.files(gs_dir, pattern = "^SFARI.*\\.csv$", full.names = TRUE)[1]
if (!is.na(sfari_file)) {
  sf <- read.csv(sfari_file, check.names = FALSE)
  cn <- tolower(colnames(sf))
  sym_col <- colnames(sf)[grepl("symbol", cn)][1]
  ens_col <- colnames(sf)[grepl("ensembl", cn)][1]
  sc_col  <- colnames(sf)[grepl("score", cn)][1]
  syn_col <- colnames(sf)[grepl("syndromic", cn)][1]
  cat("SFARI file:", basename(sfari_file), "| columns used:", sym_col, ens_col, sc_col, syn_col, "\n")
  if (any(is.na(c(sym_col, sc_col, syn_col)))) stop("SFARI columns not recognised: ", paste(colnames(sf), collapse = ", "))
  score <- suppressWarnings(as.numeric(sf[[sc_col]]))
  syn   <- suppressWarnings(as.numeric(sf[[syn_col]])) == 1
  ens   <- if (!is.na(ens_col)) sub("\\..*$", "", as.character(sf[[ens_col]])) else rep(NA, nrow(sf))
  mk <- function(i) list(ensembl = ens[i & !is.na(i)], symbol = as.character(sf[[sym_col]])[i & !is.na(i)])
  sets$SFARI_score1       <- mk(!is.na(score) & score == 1)
  sets$SFARI_score1_2     <- mk(!is.na(score) & score <= 2)
  sets$SFARI_scored_1_3   <- mk(!is.na(score) & score <= 3)
  sets$SFARI_syndromic    <- mk(!is.na(syn) & syn)
} else cat("No SFARI CSV in", gs_dir, "- SFARI sets skipped\n")

# Read a CSV that may start with a UTF-8 byte-order mark (Excel exports)
read_bom <- function(f) {
  x <- read.csv(f, check.names = FALSE, encoding = "UTF-8", stringsAsFactors = FALSE)
  colnames(x) <- trimws(sub("^(\ufeff|\u00ef\u00bb\u00bf)", "", colnames(x)))
  x
}

satterstrom_file <- file.path(gs_dir, "Satterstrom_102_ASD.csv")
if (file.exists(satterstrom_file)) {
  sa <- read_bom(satterstrom_file)
  sa <- sa[grepl("^ENSG", sa$ensembl_gene_id), ]                     # drops blank/footnote rows
  sets$Satterstrom2020 <- list(ensembl = sub("\\..*$", "", sa$ensembl_gene_id),
                               symbol  = unique(c(sa$gene, sa$hugoGene)))
} else cat("Missing", satterstrom_file, "- skipped\n")

fu_file <- file.path(gs_dir, "Fu_158_ASD.csv")
if (file.exists(fu_file)) {
  fu <- read_bom(fu_file)
  fu <- fu[!is.na(fu$ASD185) & as.logical(fu$ASD185), ]
  sets$Fu2022 <- list(ensembl = sub("\\..*$", "", fu$gene_id),
                      symbol  = unique(c(fu$gene, fu$gene_gencodeV33)))
} else cat("Missing", fu_file, "- skipped\n")
if (!length(sets)) stop("No gene sets found in ", gs_dir)
n_source <- sapply(sets, function(s) length(unique(s$ensembl[!is.na(s$ensembl) & s$ensembl != ""])))
for (nm in names(sets)) cat(sprintf("%-18s %4d genes in source list\n", nm, n_source[[nm]]))

############################################
# 2. Per subtype: Fisher + expression-matched permutation
############################################

files <- list.files(deg_dir, pattern = paste0("^", deg_prefix, ".*\\.csv$"), full.names = TRUE)
res <- list(); hits <- list(); match_rows <- list()
set.seed(seed)

for (f in files) {
  st <- sub("\\.csv$", "", sub(deg_prefix, "", basename(f)))
  tab <- read.csv(f, check.names = FALSE)
  tab$gene_id <- sub("\\..*$", "", tab$gene_id)
  tab <- tab[!is.na(tab$padj) & !is.na(tab$log2FC) & tab$gene_biotype == "protein_coding", ]
  tab <- tab[!duplicated(tab$gene_id), ]
  N <- nrow(tab)
  det <- (tab$pct_Control + tab$pct_ASD) / 2                     # detection rate (% nuclei)
  bin <- cut(rank(det, ties.method = "first"), breaks = n_bins, labels = FALSE)
  deg <- list(up   = tab$padj < fdr & tab$log2FC >=  lfc_cut,
              down = tab$padj < fdr & tab$log2FC <= -lfc_cut)
  deg$any <- deg$up | deg$down
  cat("\n====", st, "==== background (protein-coding tested):", N,
      "| DEGs up", sum(deg$up), "down", sum(deg$down), "\n")

  for (nm in names(sets)) {
    in_set <- tab$gene_id %in% sets[[nm]]$ensembl | tab$gene_name %in% sets[[nm]]$symbol
    k <- sum(in_set)
    match_rows[[paste(st, nm)]] <- data.frame(subtype = st, gene_set = nm,
                                              source_genes = n_source[[nm]],
                                              in_background = k)
    if (k < 3) next

    # Expression-matched null: sample the same number of genes per detection bin
    counts <- table(factor(bin[in_set], levels = 1:n_bins))
    null <- list(up = numeric(n_perm), down = numeric(n_perm), any = numeric(n_perm))
    for (b in which(counts > 0)) {
      pool <- which(bin == b)
      idx <- matrix(replicate(n_perm, sample(pool, counts[b])), nrow = counts[b])
      for (d in names(null)) null[[d]] <- null[[d]] + colSums(matrix(deg[[d]][idx], nrow = counts[b]))
    }

    for (d in names(deg)) {
      a <- sum(in_set & deg[[d]]); b2 <- sum(in_set & !deg[[d]])
      c2 <- sum(!in_set & deg[[d]]); d2 <- sum(!in_set & !deg[[d]])
      m <- matrix(c(a, c2, b2, d2), 2)
      ft2 <- fisher.test(m)                                      # OR and 95% CI
      ft1 <- fisher.test(m, alternative = "greater")             # one-sided p
      res[[paste(st, nm, d)]] <- data.frame(
        subtype = st, gene_set = nm, direction = d,
        set_in_background = k, DEGs = sum(deg[[d]]), overlap = a,
        expected_random = round(k * sum(deg[[d]]) / N, 2),
        odds_ratio = unname(ft2$estimate), OR_CI_low = ft2$conf.int[1], OR_CI_high = ft2$conf.int[2],
        p_fisher = ft1$p.value,
        expected_matched = round(mean(null[[d]]), 2),
        p_matched = (1 + sum(null[[d]] >= a)) / (n_perm + 1))
      if (a > 0) {
        g <- tab[in_set & deg[[d]], c("gene_id", "gene_name", "log2FC", "padj", "pct_Control", "pct_ASD")]
        if (d != "any") hits[[paste(st, nm, d)]] <- cbind(subtype = st, gene_set = nm, direction = d, g)
      }
    }
  }
}

############################################
# 3. Multiple testing, tables
############################################

R <- do.call(rbind, res)
R$fdr_fisher  <- p.adjust(R$p_fisher, "BH")
R$fdr_matched <- p.adjust(R$p_matched, "BH")
R <- R[order(R$fdr_fisher, R$p_fisher), ]
write.csv(R, file.path(out, "risk_gene_enrichment.csv"), row.names = FALSE)
write.csv(do.call(rbind, match_rows), file.path(out, "gene_set_matching.csv"), row.names = FALSE)

H <- do.call(rbind, hits)
if (!is.null(H)) {
  H <- H[order(H$gene_set, H$subtype, H$padj), ]
  write.csv(H, file.path(out, "risk_genes_that_are_DEGs.csv"), row.names = FALSE)
  # One row per gene: which subtypes/directions
  u <- unique(H[, c("gene_name", "subtype", "direction")])
  where <- tapply(paste0(u$subtype, " (", u$direction, ")"), u$gene_name, paste, collapse = "; ")
  sets_of <- tapply(H$gene_set, H$gene_name, function(x) paste(sort(unique(x)), collapse = "; "))
  agg <- data.frame(gene_name = names(where), DEG_in = as.vector(where),
                    gene_sets = as.vector(sets_of[names(where)]))
  agg$n_subtypes <- as.vector(tapply(u$subtype, u$gene_name, function(x) length(unique(x)))[agg$gene_name])
  agg <- agg[order(-agg$n_subtypes, agg$gene_name), ]
  write.csv(agg, file.path(out, "risk_DEGs_by_gene.csv"), row.names = FALSE)
}

cat("\nSignificant (Fisher FDR <", fdr, "):\n")
print(R[R$fdr_fisher < fdr, c("subtype", "gene_set", "direction", "overlap", "expected_matched",
                               "odds_ratio", "fdr_fisher", "fdr_matched")], row.names = FALSE)

############################################
# 4. Heatmap: log2 odds ratio; * Fisher FDR < 0.05, + expression-matched FDR < 0.05
############################################

subtype_levels <- c("Ex-UL", "Ex-DL", "PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+",
                    "Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Endothelial")
R$subtype   <- factor(R$subtype, levels = intersect(subtype_levels, unique(R$subtype)))
R$gene_set  <- factor(R$gene_set, levels = rev(names(sets)))
R$direction <- factor(R$direction, levels = c("up", "down", "any"),
                      labels = c("Up in ASD", "Down in ASD", "All DEGs"))
R$log2OR <- log2(pmax(R$odds_ratio, 1 / 16))
R$log2OR[is.infinite(R$log2OR)] <- 4
R$log2OR <- pmin(R$log2OR, 4)
R$mark <- paste0(ifelse(R$fdr_fisher < fdr, "*", ""), ifelse(R$fdr_matched < fdr, "+", ""))

p <- ggplot(R, aes(subtype, gene_set, fill = log2OR)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = mark), size = 3) +
  facet_grid(. ~ direction) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0,
                       limits = c(-4, 4), name = "log2 OR") +
  labs(title = "ASD risk-gene enrichment among MAST DEGs",
       caption = paste0("* Fisher FDR < ", fdr, ";  + expression-matched permutation FDR < ", fdr,
                        " (", n_perm, " permutations)\nBackground: protein-coding genes tested in each subtype; ",
                        "DEGs: FDR < ", fdr, ", |log2FC| >= ", lfc_cut, "; log2 OR capped at +/-4"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"), plot.caption = element_text(size = 6, color = "grey40"))
ggsave(file.path(out, "risk_gene_enrichment_heatmap.pdf"), p, width = 10, height = 3.2)
ggsave(file.path(out, "risk_gene_enrichment_heatmap.png"), p, width = 10, height = 3.2, dpi = 600)

writeLines(c(
  paste("DEG source:", deg_dir, "(MAST, fixed effects)"),
  paste0("DEGs: FDR < ", fdr, " and |log2FC| >= ", lfc_cut, "; up, down and all DEGs tested"),
  "Background: protein-coding genes tested in each subtype (non-NA padj)",
  "Gene matching: Ensembl ID or gene symbol; Fu 2022 set = ASD185 == TRUE; Satterstrom 2020 = 102 genes",
  paste("SFARI file:", ifelse(is.na(sfari_file), "none", basename(sfari_file)),
        "; sets: score 1; score <= 2; score <= 3; syndromic = 1"),
  "Fisher: one-sided (greater) p; odds ratio and 95% CI from two-sided test",
  paste0("Expression-matched null: ", n_perm, " random gene sets with the same counts per detection-rate bin (",
         n_bins, " quantile bins of mean % nuclei detected); p = (1 + #null >= observed) / (n + 1); seed ", seed),
  "BH FDR across all subtype x gene set x direction tests, separately for Fisher and permutation p-values"),
  file.path(out, "risk_gene_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_16.txt"))
cat("\nSaved", out, "/\n")
