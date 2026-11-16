############################################
# 16b_risk_overlap.R - descriptive overlap of ASD risk-gene lists with MAST DEGs (no statistics)
# Companion to 16_risk_genes.R (the statistical test). Same gene lists, same DEG definition.
# Run in R from the folder containing 14_DEG/ and risk_genes/
# Packages: ggplot2
############################################

library(ggplot2)

############################################
# Settings
############################################

out        <- "16_risk_genes/overlap"
deg_dir    <- "14_DEG/MAST"
deg_prefix <- "DEG_MAST_"
gs_dir     <- "risk_genes"
fdr        <- 0.05
lfc_cut    <- 0.25
gene_level_sets <- c("SFARI_score1", "Satterstrom2020", "Fu2022")   # sets shown in the gene-level heatmap
dir.create(out, recursive = TRUE, showWarnings = FALSE)

subtype_levels <- c("Ex-UL", "Ex-DL", "PVALB+", "SST+", "VIP+", "LAMP5+", "PAX6+",
                    "Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Endothelial")

############################################
# 1. Gene sets (same definitions as 16_risk_genes.R)
############################################

read_bom <- function(f) {
  x <- read.csv(f, check.names = FALSE, encoding = "UTF-8", stringsAsFactors = FALSE)
  colnames(x) <- trimws(sub("^[^A-Za-z0-9_.]+", "", colnames(x), useBytes = TRUE))   # strip byte-order mark
  x
}
sets <- list()

sfari_file <- list.files(gs_dir, pattern = "^SFARI.*\\.csv$", full.names = TRUE)[1]
if (!is.na(sfari_file)) {
  sf <- read.csv(sfari_file, check.names = FALSE)
  cn <- tolower(colnames(sf))
  sym <- as.character(sf[[colnames(sf)[grepl("symbol", cn)][1]]])
  ens <- sub("\\..*$", "", as.character(sf[[colnames(sf)[grepl("ensembl", cn)][1]]]))
  score <- suppressWarnings(as.numeric(sf[[colnames(sf)[grepl("score", cn)][1]]]))
  syn <- suppressWarnings(as.numeric(sf[[colnames(sf)[grepl("syndromic", cn)][1]]])) == 1
  mk <- function(i) { i <- i & !is.na(i); list(ensembl = ens[i], symbol = sym[i]) }
  sets$SFARI_score1     <- mk(!is.na(score) & score == 1)
  sets$SFARI_score1_2   <- mk(!is.na(score) & score <= 2)
  sets$SFARI_scored_1_3 <- mk(!is.na(score) & score <= 3)
  sets$SFARI_syndromic  <- mk(syn)
  sets_sfari_any <- mk((!is.na(score)) | syn)            # for the list-sharing summary
}
f <- file.path(gs_dir, "Satterstrom_102_ASD.csv")
if (file.exists(f)) {
  sa <- read_bom(f); sa <- sa[grepl("^ENSG", sa$ensembl_gene_id), ]
  sets$Satterstrom2020 <- list(ensembl = sub("\\..*$", "", sa$ensembl_gene_id), symbol = unique(c(sa$gene, sa$hugoGene)))
}
f <- file.path(gs_dir, "Fu_158_ASD.csv")
if (file.exists(f)) {
  fu <- read_bom(f); fu <- fu[!is.na(fu$ASD185) & as.logical(fu$ASD185), ]
  sets$Fu2022 <- list(ensembl = sub("\\..*$", "", fu$gene_id), symbol = unique(c(fu$gene, fu$gene_gencodeV33)))
}
stopifnot(length(sets) > 0)
in_set <- function(ids, syms, s) ids %in% s$ensembl | syms %in% s$symbol

############################################
# 2. DEGs per subtype (protein-coding, tested)
############################################

files <- list.files(deg_dir, pattern = paste0("^", deg_prefix, ".*\\.csv$"), full.names = TRUE)
ov_rows <- list(); deg_all <- list()

for (f in files) {
  st <- sub("\\.csv$", "", sub(deg_prefix, "", basename(f)))
  tab <- read.csv(f, check.names = FALSE)
  tab$gene_id <- sub("\\..*$", "", tab$gene_id)
  tab <- tab[!is.na(tab$padj) & !is.na(tab$log2FC) & tab$gene_biotype == "protein_coding", ]
  tab <- tab[!duplicated(tab$gene_id), ]
  up <- tab$padj < fdr & tab$log2FC >= lfc_cut
  dn <- tab$padj < fdr & tab$log2FC <= -lfc_cut
  deg_all[[st]] <- data.frame(subtype = st, gene_id = tab$gene_id[up | dn], gene_name = tab$gene_name[up | dn],
                              log2FC = tab$log2FC[up | dn], padj = tab$padj[up | dn])
  for (nm in names(sets)) {
    m <- in_set(tab$gene_id, tab$gene_name, sets[[nm]])
    ov_rows[[paste(st, nm)]] <- data.frame(
      subtype = st, gene_set = nm, set_tested = sum(m),
      up = sum(m & up), down = sum(m & dn), total = sum(m & (up | dn)),
      pct_up = round(100 * sum(m & up) / max(sum(m), 1), 1),
      pct_down = round(100 * sum(m & dn) / max(sum(m), 1), 1),
      pct_total = round(100 * sum(m & (up | dn)) / max(sum(m), 1), 1),
      genes_up = paste(sort(tab$gene_name[m & up]), collapse = ", "),
      genes_down = paste(sort(tab$gene_name[m & dn]), collapse = ", "))
  }
}
OV <- do.call(rbind, ov_rows)
OV$subtype <- factor(OV$subtype, levels = intersect(subtype_levels, unique(OV$subtype)))
OV <- OV[order(OV$gene_set, OV$subtype), ]
write.csv(OV, file.path(out, "risk_DEG_overlap_table.csv"), row.names = FALSE)
cat("Overlap (total DEGs, % of tested list genes):\n")
print(reshape(OV[, c("subtype", "gene_set", "total")], idvar = "subtype", timevar = "gene_set",
              direction = "wide"), row.names = FALSE)

############################################
# 3. Overlap heatmap: count in cell, colour = % of tested list genes
############################################

long <- rbind(
  data.frame(OV[, c("subtype", "gene_set")], direction = "Up in ASD",   n = OV$up,    pct = OV$pct_up),
  data.frame(OV[, c("subtype", "gene_set")], direction = "Down in ASD", n = OV$down,  pct = OV$pct_down),
  data.frame(OV[, c("subtype", "gene_set")], direction = "All DEGs",    n = OV$total, pct = OV$pct_total))
long$direction <- factor(long$direction, levels = c("Up in ASD", "Down in ASD", "All DEGs"))
long$gene_set <- factor(long$gene_set, levels = rev(names(sets)))

p <- ggplot(long, aes(subtype, gene_set, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n), size = 2.3) +
  facet_grid(. ~ direction) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "% of list\n(tested genes)") +
  labs(title = "Risk genes that are MAST DEGs (counts)", x = NULL, y = NULL,
       caption = paste0("Number in cell = risk genes that are DEGs (FDR < ", fdr, ", |log2FC| >= ", lfc_cut,
                        "); colour = % of the list's genes tested in that subtype. Descriptive only.")) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"), plot.caption = element_text(size = 6, color = "grey40"))
ggsave(file.path(out, "risk_DEG_overlap_heatmap.pdf"), p, width = 10, height = 3.2)
ggsave(file.path(out, "risk_DEG_overlap_heatmap.png"), p, width = 10, height = 3.2, dpi = 600)

############################################
# 4. Gene-level heatmap (risk genes that are DEGs in >= 1 subtype)
############################################

D <- do.call(rbind, deg_all)
gl_sets <- intersect(gene_level_sets, names(sets))
D$member <- Reduce(`|`, lapply(gl_sets, function(nm) in_set(D$gene_id, D$gene_name, sets[[nm]])))
G <- D[D$member, ]

if (nrow(G)) {
  n_sub <- tapply(G$subtype, G$gene_name, function(x) length(unique(x)))
  mean_fc <- tapply(G$log2FC, G$gene_name, mean)
  ord <- names(n_sub)[order(-n_sub, -mean_fc)]
  G$gene_name <- factor(G$gene_name, levels = rev(ord))
  G$subtype <- factor(G$subtype, levels = intersect(subtype_levels, unique(D$subtype)))
  lim <- min(max(abs(G$log2FC)), 2)

  # Membership strip: which list(s) each gene is in
  ids <- unique(G[, c("gene_name", "gene_id")])
  memb <- do.call(rbind, lapply(gl_sets, function(nm)
    data.frame(gene_name = ids$gene_name, list = nm, member = in_set(ids$gene_id, as.character(ids$gene_name), sets[[nm]]))))
  memb <- memb[memb$member, ]
  memb$gene_name <- factor(memb$gene_name, levels = levels(G$gene_name))
  memb$list <- factor(memb$list, levels = gl_sets)

  nr <- length(ord)
  txt <- if (nr > 120) 3.5 else if (nr > 60) 5 else 6
  sub_lv <- intersect(subtype_levels, unique(D$subtype))
  G$x <- factor(as.character(G$subtype), levels = c(sub_lv, gl_sets)); G$panel <- "log2FC by subtype"
  memb$x <- factor(as.character(memb$list), levels = c(sub_lv, gl_sets)); memb$panel <- "List"
  blank <- data.frame(x = factor(sub_lv, levels = c(sub_lv, gl_sets)), panel = "log2FC by subtype",
                      gene_name = factor(ord[1], levels = levels(G$gene_name)))
  G$panel <- factor(G$panel, levels = c("log2FC by subtype", "List"))
  memb$panel <- factor(memb$panel, levels = levels(G$panel)); blank$panel <- factor(blank$panel, levels = levels(G$panel))

  p1 <- ggplot(G, aes(x, gene_name)) +
    geom_blank(data = blank) +
    geom_tile(aes(fill = pmax(pmin(log2FC, lim), -lim)), color = "grey90", linewidth = 0.2) +
    geom_tile(data = memb, fill = "grey30", color = "white", linewidth = 0.2) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
    scale_y_discrete(drop = FALSE) +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0,
                         limits = c(-lim, lim), name = "log2FC\n(ASD vs Control)") +
    labs(x = NULL, y = NULL,
         title = "Risk genes that are MAST DEGs",
         caption = paste0("Lists: ", paste(gl_sets, collapse = ", "), ". Blank = not a DEG (FDR < ", fdr,
                          ", |log2FC| >= ", lfc_cut, ")\nRows ordered by number of subtypes; colour capped at +/-",
                          round(lim, 2))) +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = txt, face = "italic"),
          panel.grid = element_blank(), plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"), plot.caption = element_text(size = 5.5, color = "grey40"))

  h <- max(4, 0.09 * nr + 1.8)
  ggsave(file.path(out, "risk_DEG_gene_heatmap.pdf"), p1, width = 7, height = h, limitsize = FALSE)
  ggsave(file.path(out, "risk_DEG_gene_heatmap.png"), p1, width = 7, height = h, dpi = 600, limitsize = FALSE)

  wide <- reshape(D[D$member, c("gene_name", "subtype", "log2FC")], idvar = "gene_name",
                  timevar = "subtype", direction = "wide")
  colnames(wide) <- sub("^log2FC\\.", "", colnames(wide))
  wide$n_subtypes <- as.vector(n_sub[as.character(wide$gene_name)])
  wide <- wide[order(-wide$n_subtypes, wide$gene_name), ]
  write.csv(wide, file.path(out, "risk_DEG_gene_by_subtype_log2FC.csv"), row.names = FALSE)
  cat("\nGene-level heatmap:", nr, "risk genes are DEGs in >= 1 subtype\n")
}

############################################
# 5. Sharing across lists (DEG risk genes in any subtype)
############################################

U <- unique(D[, c("gene_id", "gene_name")])
lists3 <- list()
if (exists("sets_sfari_any")) lists3$SFARI <- sets_sfari_any
if (!is.null(sets$Satterstrom2020)) lists3$Satterstrom <- sets$Satterstrom2020
if (!is.null(sets$Fu2022)) lists3$Fu <- sets$Fu2022
M <- sapply(lists3, function(s) in_set(U$gene_id, U$gene_name, s))
keep <- rowSums(M) > 0
combo <- apply(M[keep, , drop = FALSE], 1, function(r) paste(names(lists3)[r], collapse = " + "))
share <- as.data.frame(table(combination = combo), stringsAsFactors = FALSE)
share <- share[order(-share$Freq), ]
share$genes <- sapply(share$combination, function(cmb) paste(sort(U$gene_name[keep][combo == cmb]), collapse = ", "))
write.csv(share, file.path(out, "risk_DEG_list_sharing.csv"), row.names = FALSE)
cat("\nRisk genes that are DEGs in >= 1 subtype, by list membership:\n")
print(share[, c("combination", "Freq")], row.names = FALSE)

share$combination <- factor(share$combination, levels = rev(share$combination))
p3 <- ggplot(share, aes(combination, Freq)) + geom_col(fill = "grey40", width = 0.7) +
  geom_text(aes(label = Freq), hjust = -0.2, size = 2.5) + coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "DEG risk genes by list membership", x = NULL, y = "Genes (DEG in >= 1 subtype)",
       caption = "SFARI = any scored or syndromic gene") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold"),
                                       plot.caption = element_text(size = 6, color = "grey40"))
ggsave(file.path(out, "risk_DEG_list_sharing.pdf"), p3, width = 5, height = 2.5)
ggsave(file.path(out, "risk_DEG_list_sharing.png"), p3, width = 5, height = 2.5, dpi = 600)

cat("\nSaved", out, "/\n")
