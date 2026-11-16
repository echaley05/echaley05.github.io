############################################
# 17_hdWGCNA.R - co-expression modules per cell group (hdWGCNA), donor-level Dx testing, module annotation
# Run in R from the folder containing 12_annotated.rds, gene_table_biotype.csv, 14_DEG/ and risk_genes/
# Packages: Seurat, hdWGCNA, WGCNA, clusterProfiler, org.Hs.eg.db, ggplot2, patchwork
#
# Run ONE network at a time: set `network` below and source the script.
# After all networks are done, set network <- "SUMMARY" to combine results across networks.
############################################

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
enableWGCNAThreads()   # no effect on Windows R GUI, harmless

############################################
# Settings
############################################

network <- "PVALB+"   # one of names(networks), or "SUMMARY"

networks <- list(       # mixed design: inhibitory subtypes, excitatory groups, glia by cell type
  "Ex-UL" = c("subtype", "Ex-UL"),   "Ex-DL" = c("subtype", "Ex-DL"),
  "PVALB+" = c("subtype", "PVALB+"), "SST+" = c("subtype", "SST+"),   "VIP+" = c("subtype", "VIP+"),
  "LAMP5+" = c("subtype", "LAMP5+"), "PAX6+" = c("subtype", "PAX6+"),
  "Astrocyte" = c("cell_type", "Astrocyte"), "Oligodendrocyte" = c("cell_type", "Oligodendrocyte"),
  "OPC" = c("cell_type", "OPC"), "Microglia" = c("cell_type", "Microglia"))

out          <- "17_hdWGCNA"
min_frac     <- 0.05      # gene expressed in >= 5% of nuclei in the group
biotypes     <- c("protein_coding", "lncRNA")
k_meta       <- 25        # nuclei per metacell
max_shared   <- 10
min_cells    <- 50        # donor x group combinations with fewer nuclei get no metacells
min_metacells_flag <- 100
min_nuclei_donor   <- 10  # donor included in the Dx test if it has >= 10 nuclei in the group
n_hubs       <- 10
fdr          <- 0.05
lfc_cut      <- 0.25      # MAST DEG definition (same as 16)
seed         <- 42
gs_dir       <- "risk_genes"

dir.create(out, showWarnings = FALSE)
gt <- read.csv("gene_table_biotype.csv")
sym <- function(ids) { s <- gt$gene_name[match(ids, gt$gene_id)]; ifelse(is.na(s), ids, s) }

############################################
# Risk gene sets (same definitions as 16_risk_genes.R)
############################################

read_bom <- function(f) {
  x <- read.csv(f, check.names = FALSE, encoding = "UTF-8", stringsAsFactors = FALSE)
  colnames(x) <- trimws(sub("^[^A-Za-z0-9_.]+", "", colnames(x), useBytes = TRUE))
  x
}
risk <- list()
sf_file <- list.files(gs_dir, pattern = "^SFARI.*\\.csv$", full.names = TRUE)[1]
if (!is.na(sf_file)) {
  sf <- read.csv(sf_file, check.names = FALSE); cn <- tolower(colnames(sf))
  ens <- sub("\\..*$", "", sf[[colnames(sf)[grepl("ensembl", cn)][1]]])
  sc  <- suppressWarnings(as.numeric(sf[[colnames(sf)[grepl("score", cn)][1]]]))
  syn <- suppressWarnings(as.numeric(sf[[colnames(sf)[grepl("syndromic", cn)][1]]])) == 1
  risk$SFARI_score1 <- ens[!is.na(sc) & sc == 1]
  risk$SFARI_all    <- ens[(!is.na(sc) & sc <= 3) | (!is.na(syn) & syn)]
}
if (file.exists(file.path(gs_dir, "Satterstrom_102_ASD.csv"))) {
  sa <- read_bom(file.path(gs_dir, "Satterstrom_102_ASD.csv"))
  risk$Satterstrom2020 <- sub("\\..*$", "", sa$ensembl_gene_id[grepl("^ENSG", sa$ensembl_gene_id)])
}
if (file.exists(file.path(gs_dir, "Fu_158_ASD.csv"))) {
  fu <- read_bom(file.path(gs_dir, "Fu_158_ASD.csv"))
  risk$Fu2022 <- sub("\\..*$", "", fu$gene_id[!is.na(fu$ASD185) & as.logical(fu$ASD185)])
}

fisher_row <- function(set, hits, universe) {
  set <- intersect(set, universe); a <- sum(hits %in% set)
  m <- matrix(c(a, length(hits) - a, length(set) - a, length(universe) - length(set) - length(hits) + a), 2)
  ft <- fisher.test(m, alternative = "greater")
  c(overlap = a, set_size = length(set), odds_ratio = unname(fisher.test(m)$estimate), p = ft$p.value)
}

############################################
# A. One network
############################################

if (network != "SUMMARY") {
  stopifnot(network %in% names(networks))
  col <- networks[[network]][1]; val <- networks[[network]][2]
  tag <- gsub("[^A-Za-z0-9_-]", "", network)
  nd <- file.path(out, tag); dir.create(nd, showWarnings = FALSE)
  set.seed(seed)

  seu <- readRDS("12_annotated.rds")
  seu <- subset(seu, cells = colnames(seu)[seu[[col, drop = TRUE]] == val])
  seu$wgcna_group <- network
  cat("Network:", network, "|", ncol(seu), "nuclei,", length(unique(seu$Sample_ID)), "donors\n")

  # ---- Genes: protein-coding + lncRNA, >= 5% of nuclei, no MT or RPS/RPL ----
  counts <- LayerData(seu, assay = "RNA", layer = "counts")
  frac <- Matrix::rowMeans(counts > 0)
  g <- gt[match(rownames(counts), gt$gene_id), ]
  keep <- frac >= min_frac & g$gene_biotype %in% biotypes &
          !grepl("^MT-", g$gene_name) & !grepl("^RPS|^RPL", g$gene_name)
  genes <- rownames(counts)[keep]
  cat("Network genes:", length(genes), "| by biotype:\n"); print(table(g$gene_biotype[keep]))
  rm(counts); gc()

  # ---- Metacells (within donor) ----
  seu <- SetupForWGCNA(seu, gene_select = "custom", gene_list = genes, wgcna_name = tag)
  seu <- MetacellsByGroups(seu, group.by = c("wgcna_group", "Sample_ID"), reduction = "scvi",
                           k = k_meta, max_shared = max_shared, min_cells = min_cells,
                           ident.group = "wgcna_group")
  seu <- NormalizeMetacells(seu)
  n_meta <- ncol(GetMetacellObject(seu))
  cat("Metacells:", n_meta, if (n_meta < min_metacells_flag) "  <-- FEWER THAN 100: treat network with caution" else "", "\n")

  seu <- SetDatExpr(seu, group_name = network, group.by = "wgcna_group")

  # ---- Soft power and network ----
  seu <- TestSoftPowers(seu, networkType = "signed")
  pt <- GetPowerTable(seu)
  write.csv(pt, file.path(nd, "soft_power_table.csv"), row.names = FALSE)
  seu <- ConstructNetwork(seu, tom_name = tag, tom_outdir = file.path(nd, "TOM"), overwrite_tom = TRUE)
  power_used <- tryCatch(GetNetworkData(seu)$params$power, error = function(e) NA)
  cat("Soft power used:", power_used, "\n")

  # ---- Module eigengenes and connectivity ----
  seu <- ScaleData(seu, features = GetWGCNAGenes(seu), verbose = FALSE)
  seu <- ModuleEigengenes(seu)
  seu <- ModuleConnectivity(seu, group.by = "wgcna_group", group_name = network)

  mods <- GetModules(seu)
  mods$gene_id <- mods$gene_name
  mods$gene_symbol <- sym(mods$gene_id)
  mods$gene_biotype <- gt$gene_biotype[match(mods$gene_id, gt$gene_id)]
  mods <- mods[, c("gene_id", "gene_symbol", "gene_biotype", setdiff(colnames(mods), c("gene_id", "gene_symbol", "gene_biotype", "gene_name")))]
  write.csv(mods, file.path(nd, "modules_genes.csv"), row.names = FALSE)
  module_names <- setdiff(levels(factor(mods$module)), "grey")
  cat("Modules (excluding grey):", length(module_names), "\n")
  print(table(mods$module))

  hubs <- GetHubGenes(seu, n_hubs = n_hubs)
  hubs$gene_symbol <- sym(hubs$gene_name)
  write.csv(hubs, file.path(nd, "hub_genes.csv"), row.names = FALSE)

  # ---- Donor-level Dx test: mean eigengene per donor ~ Dx + Age + PMI ----
  MEs <- GetMEs(seu, harmonized = FALSE)
  MEs <- MEs[, module_names, drop = FALSE]
  md <- seu@meta.data[rownames(MEs), c("Sample_ID", "Dx", "Age", "PMI")]
  n_donor <- table(md$Sample_ID)
  ok <- md$Sample_ID %in% names(n_donor)[n_donor >= min_nuclei_donor]
  sums <- rowsum(as.matrix(MEs[ok, , drop = FALSE]), as.character(md$Sample_ID[ok]))
  dme <- data.frame(Sample_ID = rownames(sums), sums / as.vector(table(as.character(md$Sample_ID[ok]))[rownames(sums)]),
                    check.names = FALSE)
  di <- unique(md[, c("Sample_ID", "Dx", "Age", "PMI")])
  dme <- cbind(di[match(dme$Sample_ID, di$Sample_ID), ], dme[, -1, drop = FALSE])
  dme$Dx <- factor(dme$Dx, levels = c("Control", "ASD"))
  write.csv(dme, file.path(nd, "module_eigengenes_donor_means.csv"), row.names = FALSE)

  dx <- do.call(rbind, lapply(module_names, function(m) {
    f <- lm(dme[[m]] ~ Dx + scale(Age) + scale(PMI), data = dme)
    s <- summary(f)$coefficients["DxASD", ]
    data.frame(network = network, module = m, n_genes = sum(mods$module == m),
               donors_Control = sum(dme$Dx == "Control"), donors_ASD = sum(dme$Dx == "ASD"),
               estimate = s[1], se = s[2], t = s[3], p = s[4])
  }))
  dx$fdr_within_network <- p.adjust(dx$p, "BH")
  dx$hub_genes <- sapply(dx$module, function(m) paste(head(hubs$gene_symbol[hubs$module == m], n_hubs), collapse = ", "))

  # ---- Module annotation: risk genes, MAST DEGs, GO BP ----
  universe <- mods$gene_id
  universe_pc <- mods$gene_id[mods$gene_biotype == "protein_coding"]
  deg_file <- file.path("14_DEG/MAST", paste0("DEG_MAST_", val, ".csv"))
  deg <- NULL
  if (file.exists(deg_file)) {
    dt <- read.csv(deg_file); dt$gene_id <- sub("\\..*$", "", dt$gene_id)
    deg <- list(up   = dt$gene_id[!is.na(dt$padj) & dt$padj < fdr & dt$log2FC >=  lfc_cut],
                down = dt$gene_id[!is.na(dt$padj) & dt$padj < fdr & dt$log2FC <= -lfc_cut])
  } else cat("No MAST table for", val, "- DEG overlap skipped\n")

  ann <- list(); go_all <- list()
  for (m in module_names) {
    mg <- mods$gene_id[mods$module == m]
    for (nm in names(risk)) {
      r <- fisher_row(risk[[nm]], intersect(mg, universe_pc), universe_pc)
      ann[[paste(m, nm)]] <- data.frame(module = m, test = nm, t(r))
    }
    if (!is.null(deg)) for (d in names(deg)) {
      r <- fisher_row(intersect(deg[[d]], universe), mg, universe)
      ann[[paste(m, d)]] <- data.frame(module = m, test = paste0("MAST_DEG_", d), t(r))
    }
    ego <- tryCatch(enrichGO(mg, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                             pvalueCutoff = fdr, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = TRUE),
                    error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego))) go_all[[m]] <- cbind(module = m, as.data.frame(ego))
  }
  A <- do.call(rbind, ann)
  A$fdr <- p.adjust(A$p, "BH")
  write.csv(A, file.path(nd, "module_enrichment_risk_and_DEGs.csv"), row.names = FALSE)
  GO <- do.call(rbind, go_all)
  if (!is.null(GO)) write.csv(GO, file.path(nd, "module_GO_BP.csv"), row.names = FALSE)
  dx$top_GO_BP <- sapply(dx$module, function(m) if (!is.null(GO) && any(GO$module == m))
    paste(head(GO$Description[GO$module == m], 3), collapse = "; ") else "")
  write.csv(dx, file.path(nd, "module_Dx_association.csv"), row.names = FALSE)
  cat("\nModule ~ Dx (donor level):\n")
  print(dx[, c("module", "n_genes", "estimate", "p", "fdr_within_network", "hub_genes")], row.names = FALSE)

  # ---- Figures ----
  pdf(file.path(nd, paste0("hdWGCNA_", tag, ".pdf")), width = 9, height = 6)
  print(wrap_plots(PlotSoftPowers(seu), ncol = 2) + plot_annotation(title = paste(network, "- soft power")))
  PlotDendrogram(seu, main = paste(network, "hdWGCNA dendrogram"))

  long <- do.call(rbind, lapply(module_names, function(m)
    data.frame(module = m, Sample_ID = dme$Sample_ID, Dx = dme$Dx, ME = dme[[m]])))
  long$label <- paste0(long$module, "\nFDR = ", signif(dx$fdr_within_network[match(long$module, dx$module)], 2))
  print(ggplot(long, aes(Dx, ME, fill = Dx)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
          geom_jitter(width = 0.15, size = 1) +
          facet_wrap(~ label, scales = "free_y") +
          scale_fill_manual(values = c(Control = "#4DBBD5", ASD = "#E64B35")) +
          labs(title = paste(network, "- module eigengene per donor (mean over nuclei)"),
               subtitle = "Test: donor mean ME ~ Dx + Age + PMI", x = NULL, y = "Module eigengene") +
          theme_classic(base_size = 8) + theme(legend.position = "none"))

  hp <- hubs; hp$gene_symbol <- factor(hp$gene_symbol, levels = rev(unique(hp$gene_symbol)))
  print(ggplot(hp, aes(kME, gene_symbol, fill = module)) + geom_col(show.legend = FALSE) +
          facet_wrap(~ module, scales = "free_y") +
          scale_fill_manual(values = with(mods[!duplicated(mods$module), ], setNames(color, module))) +
          labs(title = paste(network, "- top", n_hubs, "hub genes per module (kME)"), y = NULL) +
          theme_classic(base_size = 7) + theme(axis.text.y = element_text(face = "italic")))
  dev.off()

  saveRDS(seu, file.path(nd, paste0("hdWGCNA_", tag, ".rds")))
  writeLines(c(
    paste("Network:", network, "(", col, "==", val, ")"),
    paste("Nuclei:", ncol(seu), "| donors:", length(unique(seu$Sample_ID)), "| metacells:", n_meta),
    paste0("Genes: ", paste(biotypes, collapse = " + "), ", expressed in >= ", 100 * min_frac,
           "% of nuclei, MT and RPS/RPL removed; n = ", length(genes)),
    paste0("Metacells: MetacellsByGroups(group.by = group x Sample_ID, reduction = scvi, k = ", k_meta,
           ", max_shared = ", max_shared, ", min_cells = ", min_cells, "); NormalizeMetacells"),
    paste("Network: signed; TestSoftPowers; ConstructNetwork defaults; soft power used:", power_used),
    "Module eigengenes: ModuleEigengenes (no harmonization); kME via ModuleConnectivity",
    paste0("Dx test: donor mean ME (donors with >= ", min_nuclei_donor, " nuclei) ~ Dx + scale(Age) + scale(PMI); BH within network"),
    "Annotation: Fisher (one-sided) for risk genes (protein-coding network genes as background) and MAST DEGs (network genes as background); GO BP enrichGO with network genes as background",
    paste("hdWGCNA version:", packageVersion("hdWGCNA"), "| WGCNA:", packageVersion("WGCNA"), "| seed:", seed)),
    file.path(nd, "hdWGCNA_settings.txt"))
  writeLines(capture.output(sessionInfo()), file.path(nd, "sessionInfo.txt"))
  cat("\nSaved", nd, "/\n")
}

############################################
# B. SUMMARY across networks
############################################

if (network == "SUMMARY") {
  fl <- list.files(out, pattern = "^module_Dx_association\\.csv$", recursive = TRUE, full.names = TRUE)
  stopifnot(length(fl) > 0)
  S <- do.call(rbind, lapply(fl, read.csv))
  S$fdr_global <- p.adjust(S$p, "BH")
  S <- S[order(S$fdr_global, S$p), ]
  write.csv(S, file.path(out, "all_modules_Dx_association.csv"), row.names = FALSE)
  cat("Networks:", length(unique(S$network)), "| modules:", nrow(S), "\n")
  print(S[S$fdr_global < 0.1, c("network", "module", "n_genes", "estimate", "p", "fdr_global", "hub_genes")], row.names = FALSE)

  S$network <- factor(S$network, levels = intersect(names(networks), unique(S$network)))
  S$mark <- ifelse(S$fdr_global < 0.05, "*", ifelse(S$fdr_global < 0.1, ".", ""))
  p <- ggplot(S, aes(module, network, fill = t)) +
    geom_tile(color = "white") + geom_text(aes(label = mark), size = 4) +
    facet_grid(network ~ ., scales = "free", space = "free") +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0, name = "t (ASD vs\nControl)") +
    labs(title = "Module-Dx association (donor-level, ~ Dx + Age + PMI)",
         caption = "* global FDR < 0.05; . global FDR < 0.1", x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_blank(), panel.grid = element_blank())
  ggsave(file.path(out, "all_modules_Dx_heatmap.pdf"), p, width = 10, height = 0.35 * length(levels(S$network)) + 1.5)
  cat("Saved", out, "/all_modules_Dx_association.csv and heatmap\n")
}
