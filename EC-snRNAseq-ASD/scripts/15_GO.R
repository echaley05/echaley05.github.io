############################################
# 15_GO.R - GO, Reactome and SynGO enrichment of ASD vs Control DEGs, per subtype
# Over-representation (up and down separately) + GSEA (all tested genes, ranked), for pseudobulk and MAST results
# Run in R from the folder containing 14_DEG/ (after 14a and 14b)
# Packages: clusterProfiler, ReactomePA, org.Hs.eg.db, enrichplot, readxl, ggplot2
# SynGO: download the bulk release ZIP from https://syngoportal.org (Download page), unzip, and put
#        syngo_ontologies.xlsx in the path set below; SynGO is skipped if the file is missing
# Runtime: GSEA is the slow part (roughly 1-2 min per subtype x ontology); use run_sources / run_subtypes to split runs
############################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(readxl)
  library(ggplot2)
})

############################################
# Settings
############################################

out          <- "15_GO"
fdr          <- 0.05      # DEG threshold and term threshold (BH)
lfc_cut      <- 0.25      # DEG |log2FC| threshold for over-representation (same as volcano plots)
ontologies   <- c("BP", "MF", "CC")
min_gs       <- 10        # GO term size limits
max_gs       <- 500
min_gs_syngo <- 5         # SynGO terms are smaller and more specific
min_degs     <- 10        # over-representation run only if a direction has >= 10 DEGs
simplify_cut <- 0.7       # semantic-similarity cutoff for collapsing redundant GO terms
top_n_plot   <- 15        # terms per dot plot
seed         <- 42
syngo_file   <- "SynGo/ontologies.xlsx"

sources <- list(
  pseudobulk = list(dir = "14_DEG/pseudobulk", prefix = "DEG_pseudobulk_", lfc = "log2FoldChange_shrunk"),
  MAST       = list(dir = "14_DEG/MAST",       prefix = "DEG_MAST_",       lfc = "log2FC")
)
run_sources  <- names(sources)   # e.g. "pseudobulk" to run one source
run_subtypes <- NULL             # NULL = all; or e.g. c("PVALB+", "SST+")

############################################
# SynGO gene sets (optional)
############################################

syngo <- NULL
if (file.exists(syngo_file)) {
  sg <- as.data.frame(read_excel(syngo_file))
  cn <- tolower(colnames(sg))
  cat("SynGO columns:", paste(colnames(sg), collapse = ", "), "\n")
  id_col   <- colnames(sg)[cn == "id"][1]
  name_col <- colnames(sg)[cn == "name"][1]
  gene_col <- colnames(sg)[grepl("ensembl", cn)][1]
  if (is.na(gene_col)) gene_col <- colnames(sg)[grepl("hgnc_symbol|symbol", cn)][1]
  if (any(is.na(c(id_col, name_col, gene_col))))
    stop("Could not identify SynGO id/name/gene columns; check the column names printed above")
  t2g <- do.call(rbind, lapply(seq_len(nrow(sg)), function(i) {
    g <- trimws(unlist(strsplit(as.character(sg[[gene_col]][i]), "[,;|]")))
    g <- g[g != "" & !is.na(g)]
    if (length(g)) data.frame(term = sg[[id_col]][i], gene = g) else NULL
  }))
  syngo <- list(TERM2GENE = t2g,
                TERM2NAME = unique(data.frame(term = sg[[id_col]], name = sg[[name_col]])),
                id_type = if (all(grepl("^ENSG", head(t2g$gene, 50)))) "ensembl" else "symbol")
  cat("SynGO:", length(unique(t2g$term)), "terms,", length(unique(t2g$gene)), "genes; gene IDs:",
      syngo$id_type, "(column", gene_col, ")\n")
} else {
  cat("SynGO file not found at", syngo_file, "- SynGO analyses skipped\n")
}

############################################
# Helpers
############################################

as_df <- function(x) if (is.null(x)) NULL else as.data.frame(x)
n_sig <- function(x) if (is.null(x)) 0L else sum(as.data.frame(x)$p.adjust < fdr)

safe <- function(expr, what) tryCatch(expr, error = function(e) {
  cat("  ", what, "failed:", conditionMessage(e), "\n"); NULL })

save_tab <- function(x, file) if (!is.null(x) && nrow(as.data.frame(x)) > 0)
  write.csv(as.data.frame(x), file, row.names = FALSE)

plot_ora <- function(x, title) {
  if (is.null(x) || n_sig(x) == 0) return(NULL)
  dotplot(x, showCategory = top_n_plot) + ggtitle(title) +
    theme(plot.title = element_text(size = 9, face = "bold"), axis.text.y = element_text(size = 7))
}
plot_gsea <- function(x, title) {
  if (is.null(x) || n_sig(x) == 0) return(NULL)
  dotplot(x, showCategory = top_n_plot, split = ".sign") + facet_grid(. ~ .sign) + ggtitle(title) +
    theme(plot.title = element_text(size = 9, face = "bold"), axis.text.y = element_text(size = 7))
}

############################################
# Main loop
############################################

set.seed(seed)
summary_rows <- list()
coverage_rows <- list()
go_ids <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
go_annot <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_ids, keytype = "ENSEMBL", columns = "GO"))
has_go <- unique(go_annot$ENSEMBL[!is.na(go_annot$GO)])

for (src in run_sources) {
  cfg <- sources[[src]]
  tdir <- file.path(out, src, "tables"); fdir <- file.path(out, src, "figures")
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(cfg$dir, pattern = paste0("^", cfg$prefix, ".*\\.csv$"), full.names = TRUE)
  subtypes <- sub("\\.csv$", "", sub(cfg$prefix, "", basename(files)))
  if (!is.null(run_subtypes)) { keep <- subtypes %in% run_subtypes; files <- files[keep]; subtypes <- subtypes[keep] }

  for (k in seq_along(files)) {
    st <- subtypes[k]; fs <- gsub("[^A-Za-z0-9_-]", "", st)
    cat("\n====", src, "-", st, "====\n")
    tab <- read.csv(files[k], check.names = FALSE)
    tab$gene_id <- sub("\\..*$", "", tab$gene_id)                 # drop version suffix if present
    tab <- tab[!is.na(tab$padj) & !is.na(tab$pvalue) & !is.na(tab[[cfg$lfc]]), ]
    tab <- tab[!duplicated(tab$gene_id), ]
    tab$lfc <- tab[[cfg$lfc]]

    universe <- tab$gene_id                                        # genes tested in this subtype
    degs <- list(up   = tab$gene_id[tab$padj < fdr & tab$lfc >=  lfc_cut],
                 down = tab$gene_id[tab$padj < fdr & tab$lfc <= -lfc_cut])

    # Ranking for GSEA: DESeq2 Wald statistic; MAST: sign(log2FC) * -log10(p)
    p <- tab$pvalue; p[p == 0] <- min(p[p > 0])
    rank_stat <- if ("stat" %in% colnames(tab) && !all(is.na(tab$stat))) tab$stat else sign(tab$lfc) * -log10(p)
    gl <- sort(setNames(rank_stat, tab$gene_id)[!is.na(rank_stat)], decreasing = TRUE)

    for (d in names(degs)) coverage_rows[[paste(src, st, d)]] <- data.frame(
      source = src, subtype = st, direction = d, DEGs = length(degs[[d]]),
      DEGs_with_GO = sum(degs[[d]] %in% has_go),
      protein_coding = sum(tab$gene_biotype[tab$gene_id %in% degs[[d]]] == "protein_coding"),
      lncRNA = sum(tab$gene_biotype[tab$gene_id %in% degs[[d]]] == "lncRNA"))
    cat("DEGs up:", length(degs$up), "| down:", length(degs$down), "| universe:", length(universe), "\n")

    plots <- list()
    rec <- function(analysis, direction, ont, full, simp) {
      summary_rows[[paste(src, st, analysis, direction, ont)]] <<- data.frame(
        source = src, subtype = st, analysis = analysis, direction = direction, ontology = ont,
        sig_terms = n_sig(full), sig_terms_simplified = n_sig(simp))
    }

    # ---- GO over-representation (up / down separately) ----
    for (d in names(degs)) for (ont in ontologies) {
      if (length(degs[[d]]) < min_degs) { rec("Over-representation", d, ont, NULL, NULL); next }
      r <- safe(enrichGO(gene = degs[[d]], universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                         ont = ont, pAdjustMethod = "BH", pvalueCutoff = fdr, qvalueCutoff = 1,
                         minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE), paste("Over-representation", d, ont))
      s <- if (n_sig(r) > 1) safe(simplify(r, cutoff = simplify_cut, by = "p.adjust", select_fun = min),
                                  paste("simplify over-representation", d, ont)) else r
      save_tab(r, file.path(tdir, paste0("OverRep_", fs, "_", d, "_", ont, ".csv")))
      save_tab(s, file.path(tdir, paste0("OverRep_", fs, "_", d, "_", ont, "_simplified.csv")))
      rec("Over-representation", d, ont, r, s)
      plots[[length(plots) + 1]] <- plot_ora(s, paste0(st, " - GO ", ont, " over-representation, ", d, " in ASD (", src, ")"))
    }

    # ---- GO GSEA (all tested genes) ----
    for (ont in ontologies) {
      set.seed(seed)
      r <- safe(gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = ont,
                      minGSSize = min_gs, maxGSSize = max_gs, pvalueCutoff = fdr, pAdjustMethod = "BH",
                      eps = 0, seed = TRUE, verbose = FALSE), paste("GSEA", ont))
      if (!is.null(r) && n_sig(r) > 0) r <- setReadable(r, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
      s <- if (n_sig(r) > 1) safe(simplify(r, cutoff = simplify_cut, by = "p.adjust", select_fun = min),
                                  paste("simplify GSEA", ont)) else r
      save_tab(r, file.path(tdir, paste0("GSEA_", fs, "_", ont, ".csv")))
      save_tab(s, file.path(tdir, paste0("GSEA_", fs, "_", ont, "_simplified.csv")))
      rec("GSEA", "both", ont, r, s)
      plots[[length(plots) + 1]] <- plot_gsea(s, paste0(st, " - GO ", ont, " GSEA (", src, ")"))
    }

    # ---- SynGO (over-representation up/down + GSEA) ----
    if (!is.null(syngo)) {
      key <- if (syngo$id_type == "ensembl") tab$gene_id else tab$gene_name
      map <- setNames(key, tab$gene_id)
      for (d in names(degs)) {
        if (length(degs[[d]]) < min_degs) { rec("Over-representation", d, "SynGO", NULL, NULL); next }
        r <- safe(enricher(gene = unique(map[degs[[d]]]), universe = unique(map[universe]),
                           TERM2GENE = syngo$TERM2GENE, TERM2NAME = syngo$TERM2NAME,
                           pAdjustMethod = "BH", pvalueCutoff = fdr, qvalueCutoff = 1,
                           minGSSize = min_gs_syngo, maxGSSize = max_gs), paste("SynGO over-representation", d))
        save_tab(r, file.path(tdir, paste0("SynGO_OverRep_", fs, "_", d, ".csv")))
        rec("Over-representation", d, "SynGO", r, r)
        plots[[length(plots) + 1]] <- plot_ora(r, paste0(st, " - SynGO over-representation, ", d, " in ASD (", src, ")"))
      }
      gl_s <- gl; names(gl_s) <- map[names(gl)]
      gl_s <- gl_s[!is.na(names(gl_s)) & !duplicated(names(gl_s))]
      set.seed(seed)
      r <- safe(GSEA(geneList = gl_s, TERM2GENE = syngo$TERM2GENE, TERM2NAME = syngo$TERM2NAME,
                     minGSSize = min_gs_syngo, maxGSSize = max_gs, pvalueCutoff = fdr,
                     pAdjustMethod = "BH", eps = 0, seed = TRUE, verbose = FALSE), "SynGO GSEA")
      save_tab(r, file.path(tdir, paste0("SynGO_GSEA_", fs, ".csv")))
      rec("GSEA", "both", "SynGO", r, r)
      plots[[length(plots) + 1]] <- plot_gsea(r, paste0(st, " - SynGO GSEA (", src, ")"))
    }

    # ---- Reactome (over-representation up/down + GSEA); Ensembl -> Entrez ----
    ez <- suppressMessages(suppressWarnings(
      bitr(universe, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)))
    ez <- ez[!duplicated(ez$ENSEMBL), ]
    to_ez <- setNames(ez$ENTREZID, ez$ENSEMBL)
    universe_ez <- unique(unname(to_ez[universe[universe %in% names(to_ez)]]))
    for (d in names(degs)) {
      g_ez <- unique(unname(to_ez[degs[[d]][degs[[d]] %in% names(to_ez)]]))
      if (length(g_ez) < min_degs) { rec("Over-representation", d, "Reactome", NULL, NULL); next }
      r <- safe(enrichPathway(gene = g_ez, universe = universe_ez, organism = "human",
                              pAdjustMethod = "BH", pvalueCutoff = fdr, qvalueCutoff = 1,
                              minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE),
                paste("Reactome over-representation", d))
      save_tab(r, file.path(tdir, paste0("OverRep_", fs, "_", d, "_Reactome.csv")))
      rec("Over-representation", d, "Reactome", r, r)
      plots[[length(plots) + 1]] <- plot_ora(r, paste0(st, " - Reactome over-representation, ", d, " in ASD (", src, ")"))
    }
    gl_e <- gl[names(gl) %in% names(to_ez)]
    names(gl_e) <- to_ez[names(gl_e)]
    gl_e <- sort(gl_e[!duplicated(names(gl_e))], decreasing = TRUE)
    set.seed(seed)
    r <- safe(gsePathway(geneList = gl_e, organism = "human", minGSSize = min_gs, maxGSSize = max_gs,
                         pvalueCutoff = fdr, pAdjustMethod = "BH", eps = 0, seed = TRUE, verbose = FALSE),
              "Reactome GSEA")
    if (!is.null(r) && n_sig(r) > 0) r <- setReadable(r, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    save_tab(r, file.path(tdir, paste0("GSEA_", fs, "_Reactome.csv")))
    rec("GSEA", "both", "Reactome", r, r)
    plots[[length(plots) + 1]] <- plot_gsea(r, paste0(st, " - Reactome GSEA (", src, ")"))

    plots <- Filter(Negate(is.null), plots)
    if (length(plots)) {
      pdf(file.path(fdir, paste0("GO_", fs, ".pdf")), width = 8, height = 6.5)
      for (pl in plots) print(pl)
      dev.off()
    }
    # Save progress after each subtype
    write.csv(do.call(rbind, summary_rows), file.path(out, "GO_summary.csv"), row.names = FALSE)
    write.csv(do.call(rbind, coverage_rows), file.path(out, "GO_DEG_annotation_coverage.csv"), row.names = FALSE)
  }

  ############################################
  # Cross-subtype comparison: top 5 over-represented terms per subtype x direction
  # (GO BP simplified, and Reactome)
  ############################################
  for (db in c("BP_simplified", "Reactome")) {
    pat <- paste0("^OverRep_(.*)_(up|down)_", db, "\\.csv$")
    fl <- list.files(tdir, pattern = pat, full.names = TRUE)
    if (!length(fl)) next
    all <- do.call(rbind, lapply(fl, function(f) {
      x <- read.csv(f); if (!nrow(x)) return(NULL)
      m <- regmatches(basename(f), regexec(pat, basename(f)))[[1]]
      x$subtype <- m[2]; x$direction <- m[3]
      x$ratio <- sapply(strsplit(x$GeneRatio, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
      x <- x[x$p.adjust < fdr, ]
      head(x[order(x$p.adjust), c("ID", "Description", "GeneRatio", "ratio", "p.adjust", "Count",
                                  "subtype", "direction")], 5)
    }))
    if (is.null(all) || !nrow(all)) next
    lab <- if (db == "Reactome") "Reactome" else "GO BP"
    all$Description <- factor(all$Description,
                              levels = rev(unique(all$Description[order(all$direction, all$p.adjust)])))
    all$direction <- factor(paste(all$direction, "in ASD"), levels = c("up in ASD", "down in ASD"))
    p <- ggplot(all, aes(subtype, Description, size = ratio, color = -log10(p.adjust))) +
      geom_point() + facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
      scale_color_gradient(low = "#9ECAE1", high = "#B2182B", name = "-log10 FDR") +
      scale_size(range = c(1, 4), name = "Gene ratio") +
      labs(title = paste0(lab, " over-representation across subtypes (", src, "; top 5 terms each)"),
           x = NULL, y = NULL) +
      theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                       plot.title = element_text(face = "bold"))
    tag <- gsub(" ", "_", lab)
    ggsave(file.path(fdir, paste0(tag, "_OverRep_across_subtypes.pdf")), p, width = 9,
           height = max(4, 0.16 * nrow(all) + 2))
    write.csv(all, file.path(tdir, paste0(tag, "_OverRep_across_subtypes_top5.csv")), row.names = FALSE)
  }
}

############################################
# Summary and settings
############################################

summary_tab <- do.call(rbind, summary_rows)
coverage_tab <- do.call(rbind, coverage_rows)
write.csv(summary_tab, file.path(out, "GO_summary.csv"), row.names = FALSE)
write.csv(coverage_tab, file.path(out, "GO_DEG_annotation_coverage.csv"), row.names = FALSE)
cat("\nSignificant terms (FDR <", fdr, "):\n"); print(summary_tab, row.names = FALSE)
cat("\nDEG GO annotation coverage:\n"); print(coverage_tab, row.names = FALSE)

writeLines(c(
  paste("Sources:", paste(run_sources, collapse = ", ")),
  paste0("DEGs for over-representation: padj < ", fdr, " and |log2FC| >= ", lfc_cut,
         " (pseudobulk: shrunk log2FC; MAST: MAST log2FC); up and down tested separately; run if >= ", min_degs, " DEGs"),
  "Universe: all genes tested in that subtype (non-NA padj), all biotypes",
  "GSEA ranking: pseudobulk = DESeq2 Wald statistic; MAST = sign(log2FC) * -log10(p) (p = 0 set to smallest non-zero p)",
  paste("GO: clusterProfiler enrichGO / gseGO, org.Hs.eg.db, keyType ENSEMBL, ontologies", paste(ontologies, collapse = ", ")),
  paste0("Term size ", min_gs, "-", max_gs, " (SynGO ", min_gs_syngo, "-", max_gs, "); BH FDR < ", fdr,
         "; GSEA eps = 0, seed ", seed),
  paste("Redundancy: clusterProfiler::simplify, cutoff", simplify_cut, "(GO only; full and simplified tables saved)"),
  "Reactome: ReactomePA enrichPathway / gsePathway (human), Ensembl -> Entrez via bitr (first match); not simplified",
  paste("SynGO:", if (is.null(syngo)) "not run (file missing)" else paste("from", syngo_file, "; gene IDs:", syngo$id_type)),
  paste("org.Hs.eg.db version:", packageVersion("org.Hs.eg.db"), "| clusterProfiler:", packageVersion("clusterProfiler"),
        "| ReactomePA:", packageVersion("ReactomePA"), "| reactome.db:", packageVersion("reactome.db"))),
  file.path(out, "GO_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_15.txt"))
cat("\nSaved", out, "/\n")
