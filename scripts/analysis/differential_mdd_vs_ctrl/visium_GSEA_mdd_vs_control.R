###################################
# Visium v1 and v2 Gene Set Enrichment Analysis
# Requires v1 and v2 unfiltered DEGs .csv files produced from visium_mdd_vs_control_degs.R
# By Victor Anosike
####################################

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
  library(Polychrome)
  library(scCustomize)
  library(reticulate)
  library(dplyr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(openxlsx)
})


########################### Visium v1 ######################################################
############################################################################################

## Ensure that v1 cluster names to match final decided names
v1_cluster_fullnames <- c(
  "0" = "ecm",
  "1" = "luc-rad",
  "2" = "axon",
  "3" = "sgz-ml",
  "4" = "ca.1-4.1",
  "5" = "ca.1-4.2",
  "6" = "dendr",
  "7" = "sgz-pl",
  "8" = "gcl",
  "9" = "ca.1-2",
  "10" = "inn",
  "11" = "vasc",
  "12" = "cp"
)

# Define the function for cluster naming
cluster_namer <- function(i) {
  i <- as.character(i)  # ensure numeric input works
  if (i %in% names(v1_cluster_fullnames)) {
    return(v1_cluster_fullnames[i])
  } else {
    warning(paste("Cluster", i, "not found in v1 cluster list."))
    return(NA)
  }
}

### GSEA pipeline (GO BP/CC/MF, KEGG, Reactome)
### Uses ranking_score = sign(logFC) * -log10(adj.P.Val) (your formula)
### One PDF per cell_type containing sections: CC, MF, BP, KEGG, Reactome

# --- required packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("clusterProfiler","org.Hs.eg.db","fgsea","ReactomePA","dplyr","enrichplot","ggplot2","data.table")
for (p in pkgs) {
  if (!suppressWarnings(require(p, character.only = TRUE))) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
    library(p, character.only = TRUE)
  }
}

# --- user inputs: replace these with your objects if needed ---
# degs: data.frame / tibble with at least columns: gene (SYMBOL), logFC, adj.P.Val, cell_type
# cell_types: vector/list of cell types you want to iterate through (e.g., unique(degs$cell_type))
# Example placeholders (REMOVE or REPLACE if you already have objects):
# Insert v1 unfiltered DEGs .csv generated from visium_mdd_vs_control_degs.R
degs <- read.csv('/path/to/visium_v1_degs_unfiltered_date.csv', stringsAsFactors = FALSE)
cell_types <- unique(degs$cell_type)

# For safety, require the user to have degs and cell_types already defined.
if (!exists("degs")) stop("Please provide 'degs' data.frame with columns: gene, logFC, adj.P.Val, cell_type")
if (!exists("cell_types")) cell_types <- unique(degs$cell_type)

# --- parameters ---
min_genes <- 15   # minimum gene set size
max_genes <- 500  # maximum gene set size
pvalue_cutoff <- 0.05
out_prefix <- "GSEA_results"   # prefix for saved files
#output_dir <- "./gsea_output"
output_dir <- '/path/to/folder/'
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# output accumulators
cc_all <- data.frame()
mf_all <- data.frame()
bp_all <- data.frame()
kegg_all <- data.frame()
reactome_all <- data.frame()

# reproducibility
set.seed(42)

# --- main loop ---
for (ct in cell_types) {
  
  #convert ct to name of the actual cluster
  ct_new <- cluster_namer(ct)
  
  message("\n========== Doing GSEA for cell_type: ", ct_new, " ==========")
  
  de <- degs %>% dplyr::filter(cell_type == ct)
  if (nrow(de) == 0) {
    message("No rows for ", ct_new, " — skipping.")
    next
  }
  
  # compute ranking score exactly as in your code
  de <- de %>% dplyr::mutate(ranking_score = sign(logFC) * (-log10(adj.P.Val)))
  
  # remove NA genes / NA ranking_score
  de <- de %>% dplyr::filter(!is.na(gene) & !is.na(ranking_score))
  if (nrow(de) < 3) {
    message("Too few genes after filtering for ", ct_new, " — skipping.")
    next
  }
  
  # map SYMBOL -> ENTREZID
  entrez_map <- tryCatch({
    bitr(unique(de$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("Error mapping to ENTREZ for ", ct_new, ": ", e$message)
    return(NULL)
  })
  if (is.null(entrez_map) || nrow(entrez_map) == 0) {
    message("No ENTREZ mapping for ", ct_new, " — skipping.")
    next
  }
  
  # keep only genes present in mapping, and preserve the ranking order
  de2 <- de %>% dplyr::filter(gene %in% entrez_map$SYMBOL)
  # ensure unique mapping per SYMBOL (drop duplicates)
  entrez_map <- entrez_map %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  # align entrez_map to de2 order
  de2 <- de2 %>% dplyr::left_join(entrez_map, by = c("gene" = "SYMBOL"))
  # remove any NA ENTREZID just in case
  de2 <- de2 %>% dplyr::filter(!is.na(ENTREZID))
  
  # build named numeric vector (names = ENTREZID) ordered descending by ranking_score
  geneList <- de2$ranking_score
  names(geneList) <- de2$ENTREZID
  # remove duplicates by keeping first occurrence (highest rank because de2 is in original order)
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  
  if (length(geneList) < min_genes) {
    message("Less than min_genes (", min_genes, ") after mapping for ", ct_new, " — skipping.")
    next
  }
  
  # --- run GSEA for GO: CC, MF, BP ---
  gsea_cc <- tryCatch({
    gseGO(geneList = geneList,
          OrgDb = org.Hs.eg.db,
          ont = "CC",
          minGSSize = min_genes,
          maxGSSize = max_genes,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE)
  }, error = function(e) { message("gseGO CC error: ", e$message); return(NULL) })
  
  gsea_mf <- tryCatch({
    gseGO(geneList = geneList,
          OrgDb = org.Hs.eg.db,
          ont = "MF",
          minGSSize = min_genes,
          maxGSSize = max_genes,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE)
  }, error = function(e) { message("gseGO MF error: ", e$message); return(NULL) })
  
  gsea_bp <- tryCatch({
    gseGO(geneList = geneList,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          minGSSize = min_genes,
          maxGSSize = max_genes,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE)
  }, error = function(e) { message("gseGO BP error: ", e$message); return(NULL) })
  
  # make readable (ENTREZ -> SYMBOL) if results exist
  if (!is.null(gsea_cc) && nrow(gsea_cc@result) > 0) gsea_cc <- setReadable(gsea_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  if (!is.null(gsea_mf) && nrow(gsea_mf@result) > 0) gsea_mf <- setReadable(gsea_mf, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  if (!is.null(gsea_bp) && nrow(gsea_bp@result) > 0) gsea_bp <- setReadable(gsea_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # attach cell_type meta and accumulate
  if (!is.null(gsea_cc) && nrow(gsea_cc@result) > 0) {
    tmp <- gsea_cc@result; tmp$cell_type <- ct; cc_all <- rbind(cc_all, tmp)
  }
  if (!is.null(gsea_mf) && nrow(gsea_mf@result) > 0) {
    tmp <- gsea_mf@result; tmp$cell_type <- ct; mf_all <- rbind(mf_all, tmp)
  }
  if (!is.null(gsea_bp) && nrow(gsea_bp@result) > 0) {
    tmp <- gsea_bp@result; tmp$cell_type <- ct; bp_all <- rbind(bp_all, tmp)
  }
  
  # --- run GSEA for KEGG (requires ENTREZ named vector) ---
  gsea_kegg <- tryCatch({
    gseKEGG(geneList = geneList,
            organism = "hsa",
            minGSSize = min_genes,
            maxGSSize = max_genes,
            pvalueCutoff = pvalue_cutoff,
            verbose = FALSE)
  }, error = function(e) { message("gseKEGG error: ", e$message); return(NULL) })
  
  if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
    gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    tmp <- gsea_kegg@result; tmp$cell_type <- ct; kegg_all <- rbind(kegg_all, tmp)
  }
  
  # --- run GSEA for Reactome ---
  gsea_reactome <- tryCatch({
    gsePathway(geneList = geneList,
               organism = "human",
               minGSSize = min_genes,
               maxGSSize = max_genes,
               pvalueCutoff = pvalue_cutoff,
               verbose = FALSE)
  }, error = function(e) { message("gsePathway error: ", e$message); return(NULL) })
  
  if (!is.null(gsea_reactome) && nrow(gsea_reactome@result) > 0) {
    gsea_reactome <- setReadable(gsea_reactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    tmp <- gsea_reactome@result; tmp$cell_type <- ct; reactome_all <- rbind(reactome_all, tmp)
  }
  
  # --- create one PDF for this cell type with sections ---
  
  #12-23-2025 --> Changing pdf name
  pdf_file <- file.path(output_dir, paste0(ct_new, "_v1_GSEA_BP_CC_MF_KEGG_Reactome_date.pdf"))
  
  pdf(pdf_file, width = 10, height = 8)
  message("Writing PDF: ", pdf_file)
  
  # helper to safe-plot dotplot and example enrichment curve
  safe_plot_section <- function(gsea_obj, title_label, show_n = 25) {
    if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
      plot.new(); title(main = paste(title_label, "- no significant terms"))
      return(NULL)
    }
    res_df <- gsea_obj@result %>% dplyr::arrange(p.adjust) 
    
    
    # dotplot of top terms by NES
    dp <- tryCatch({
      dotplot(gsea_obj, x = "NES", showCategory = show_n, label_format = 60) + ggtitle(paste(title_label)) +
        expand_limits(x = 0) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40")
    }, error = function(e) {
      message("dotplot error for ", title_label, ": ", e$message)
      NULL
    })
    if (!is.null(dp)) print(dp)
    # example enrichment plot (top pathway)
    top_id <- res_df$ID[1]
    if (!is.na(top_id) && nzchar(top_id)) {
      gp <- tryCatch({
        gseaplot2(gsea_obj, geneSetID = top_id, title = paste(title_label, "-", res_df$Description[1]))
      }, error = function(e) {
        message("gseaplot2 error for ", title_label, ": ", e$message)
        NULL
      })
      if (!is.null(gp)) print(gp)
    }
  }
  
  # Section: CC
  safe_plot_section(gsea_cc, paste(ct_new, " - GO CC"))
  # Section: MF
  safe_plot_section(gsea_mf, paste(ct_new, " - GO MF"))
  # Section: BP
  safe_plot_section(gsea_bp, paste(ct_new, " - GO BP"))
  # Section: KEGG
  safe_plot_section(gsea_kegg, paste(ct_new, " - KEGG"))
  # Section: Reactome
  safe_plot_section(gsea_reactome, paste(ct_new, " - Reactome"))
  
  dev.off()
  message("Finished cell_type: ", ct_new, " — PDF saved.")
}

# --- save combined CSVs for all cell types ---
#12-16-2025 --> Updated location of .csv files (although since the only the graph dimensions are changing, these should be exactly the same 
# as the 12-2-2025 .csv files)
write.csv(cc_all, file.path(output_dir, paste0(out_prefix, "_v1_GO_CC_all_celltypes_date.csv")), row.names = FALSE)
write.csv(mf_all, file.path(output_dir, paste0(out_prefix, "_v1_GO_MF_all_celltypes_date.csv")), row.names = FALSE)
write.csv(bp_all, file.path(output_dir, paste0(out_prefix, "_v1_GO_BP_all_celltypes_date.csv")), row.names = FALSE)
write.csv(kegg_all, file.path(output_dir, paste0(out_prefix, "_v1_KEGG_all_celltypes_date.csv")), row.names = FALSE)
write.csv(reactome_all, file.path(output_dir, paste0(out_prefix, "_v1_Reactome_all_celltypes_date.csv")), row.names = FALSE)

message("\nAll done. Results and PDFs are in: ", normalizePath(output_dir))

########################### Visium v2 ######################################################
############################################################################################

### GSEA pipeline (GO BP/CC/MF, KEGG, Reactome)
### Uses ranking_score = sign(logFC) * -log10(adj.P.Val) (your formula)
### One PDF per cell_type containing sections: CC, MF, BP, KEGG, Reactome

# --- required packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("clusterProfiler","org.Hs.eg.db","fgsea","ReactomePA","dplyr","enrichplot","ggplot2","data.table")
for (p in pkgs) {
  if (!suppressWarnings(require(p, character.only = TRUE))) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
    library(p, character.only = TRUE)
  }
}

# --- user inputs: replace these with your objects if needed ---
# degs: data.frame / tibble with at least columns: gene (SYMBOL), logFC, adj.P.Val, cell_type
# cell_types: vector/list of cell types you want to iterate through (e.g., unique(degs$cell_type))
# Example placeholders (REMOVE or REPLACE if you already have objects):
# Insert v2 unfiltered DEGs csv generated from visium_mdd_vs_control_degs.R
degs <- read.csv('/path/to/visium_v2_degs_unfiltered_date.csv', stringsAsFactors = FALSE)
cell_types <- unique(degs$cell_type)

# For safety, require the user to have degs and cell_types already defined.
if (!exists("degs")) stop("Please provide 'degs' data.frame with columns: gene, logFC, adj.P.Val, cell_type")
if (!exists("cell_types")) cell_types <- unique(degs$cell_type)

# --- parameters ---
min_genes <- 15   # minimum gene set size
max_genes <- 500  # maximum gene set size
pvalue_cutoff <- 0.05
out_prefix <- "GSEA_results"   # prefix for saved files
#output_dir <- "./gsea_output"
output_dir <- '/path/to/folder/'
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# output accumulators
cc_all <- data.frame()
mf_all <- data.frame()
bp_all <- data.frame()
kegg_all <- data.frame()
reactome_all <- data.frame()

# reproducibility
set.seed(42)

# --- main loop ---
for (ct in cell_types) {
  message("\n========== Doing GSEA for cell_type: ", ct, " ==========")
  
  de <- degs %>% dplyr::filter(cell_type == ct)
  if (nrow(de) == 0) {
    message("No rows for ", ct, " — skipping.")
    next
  }
  
  # compute ranking score exactly as in your code
  de <- de %>% dplyr::mutate(ranking_score = sign(logFC) * (-log10(adj.P.Val)))
  
  # remove NA genes / NA ranking_score
  de <- de %>% dplyr::filter(!is.na(gene) & !is.na(ranking_score))
  if (nrow(de) < 3) {
    message("Too few genes after filtering for ", ct, " — skipping.")
    next
  }
  
  # map SYMBOL -> ENTREZID
  entrez_map <- tryCatch({
    bitr(unique(de$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("Error mapping to ENTREZ for ", ct, ": ", e$message)
    return(NULL)
  })
  if (is.null(entrez_map) || nrow(entrez_map) == 0) {
    message("No ENTREZ mapping for ", ct, " — skipping.")
    next
  }
  
  # keep only genes present in mapping, and preserve the ranking order
  de2 <- de %>% dplyr::filter(gene %in% entrez_map$SYMBOL)
  # ensure unique mapping per SYMBOL (drop duplicates)
  entrez_map <- entrez_map %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  # align entrez_map to de2 order
  de2 <- de2 %>% dplyr::left_join(entrez_map, by = c("gene" = "SYMBOL"))
  # remove any NA ENTREZID just in case
  de2 <- de2 %>% dplyr::filter(!is.na(ENTREZID))
  
  # build named numeric vector (names = ENTREZID) ordered descending by ranking_score
  geneList <- de2$ranking_score
  names(geneList) <- de2$ENTREZID
  # remove duplicates by keeping first occurrence (highest rank because de2 is in original order)
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  
  if (length(geneList) < min_genes) {
    message("Less than min_genes (", min_genes, ") after mapping for ", ct, " — skipping.")
    next
  }
  
  # --- run GSEA for GO: CC, MF, BP ---
  gsea_cc <- tryCatch({
    gseGO(geneList = geneList,
          OrgDb = org.Hs.eg.db,
          ont = "CC",
          minGSSize = min_genes,
          maxGSSize = max_genes,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE)
  }, error = function(e) { message("gseGO CC error: ", e$message); return(NULL) })
  
  gsea_mf <- tryCatch({
    gseGO(geneList = geneList,
          OrgDb = org.Hs.eg.db,
          ont = "MF",
          minGSSize = min_genes,
          maxGSSize = max_genes,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE)
  }, error = function(e) { message("gseGO MF error: ", e$message); return(NULL) })
  
  gsea_bp <- tryCatch({
    gseGO(geneList = geneList,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          minGSSize = min_genes,
          maxGSSize = max_genes,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE)
  }, error = function(e) { message("gseGO BP error: ", e$message); return(NULL) })
  
  # make readable (ENTREZ -> SYMBOL) if results exist
  if (!is.null(gsea_cc) && nrow(gsea_cc@result) > 0) gsea_cc <- setReadable(gsea_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  if (!is.null(gsea_mf) && nrow(gsea_mf@result) > 0) gsea_mf <- setReadable(gsea_mf, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  if (!is.null(gsea_bp) && nrow(gsea_bp@result) > 0) gsea_bp <- setReadable(gsea_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # attach cell_type meta and accumulate
  if (!is.null(gsea_cc) && nrow(gsea_cc@result) > 0) {
    tmp <- gsea_cc@result; tmp$cell_type <- ct; cc_all <- rbind(cc_all, tmp)
  }
  if (!is.null(gsea_mf) && nrow(gsea_mf@result) > 0) {
    tmp <- gsea_mf@result; tmp$cell_type <- ct; mf_all <- rbind(mf_all, tmp)
  }
  if (!is.null(gsea_bp) && nrow(gsea_bp@result) > 0) {
    tmp <- gsea_bp@result; tmp$cell_type <- ct; bp_all <- rbind(bp_all, tmp)
  }
  
  # --- run GSEA for KEGG (requires ENTREZ named vector) ---
  gsea_kegg <- tryCatch({
    gseKEGG(geneList = geneList,
            organism = "hsa",
            minGSSize = min_genes,
            maxGSSize = max_genes,
            pvalueCutoff = pvalue_cutoff,
            verbose = FALSE)
  }, error = function(e) { message("gseKEGG error: ", e$message); return(NULL) })
  
  if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
    gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    tmp <- gsea_kegg@result; tmp$cell_type <- ct; kegg_all <- rbind(kegg_all, tmp)
  }
  
  # --- run GSEA for Reactome ---
  gsea_reactome <- tryCatch({
    gsePathway(geneList = geneList,
               organism = "human",
               minGSSize = min_genes,
               maxGSSize = max_genes,
               pvalueCutoff = pvalue_cutoff,
               verbose = FALSE)
  }, error = function(e) { message("gsePathway error: ", e$message); return(NULL) })
  
  if (!is.null(gsea_reactome) && nrow(gsea_reactome@result) > 0) {
    gsea_reactome <- setReadable(gsea_reactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    tmp <- gsea_reactome@result; tmp$cell_type <- ct; reactome_all <- rbind(reactome_all, tmp)
  }
  
  # --- create one PDF for this cell type with sections ---
  pdf_file <- file.path(output_dir, paste0(ct, "_v2_GSEA_BP_CC_MF_KEGG_Reactome_date.pdf"))
  pdf(pdf_file, width = 10, height = 8)
  message("Writing PDF: ", pdf_file)
  
  # helper to safe-plot dotplot and example enrichment curve
  safe_plot_section <- function(gsea_obj, title_label, show_n = 25) {
    if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
      plot.new(); title(main = paste(title_label, "- no significant terms"))
      return(NULL)
    }
    res_df <- gsea_obj@result %>% dplyr::arrange(p.adjust) 
    # dotplot of top terms by NES
    dp <- tryCatch({
      dotplot(gsea_obj, x = "NES", showCategory = show_n, label_format = 60) + ggtitle(paste(title_label))
    }, error = function(e) {
      message("dotplot error for ", title_label, ": ", e$message)
      NULL
    })
    if (!is.null(dp)) print(dp)
    # example enrichment plot (top pathway)
    top_id <- res_df$ID[1]
    if (!is.na(top_id) && nzchar(top_id)) {
      gp <- tryCatch({
        gseaplot2(gsea_obj, geneSetID = top_id, title = paste(title_label, "-", res_df$Description[1]))
      }, error = function(e) {
        message("gseaplot2 error for ", title_label, ": ", e$message)
        NULL
      })
      if (!is.null(gp)) print(gp)
    }
  }
  
  # Section: CC
  safe_plot_section(gsea_cc, paste(ct, " - GO CC"))
  # Section: MF
  safe_plot_section(gsea_mf, paste(ct, " - GO MF"))
  # Section: BP
  safe_plot_section(gsea_bp, paste(ct, " - GO BP"))
  # Section: KEGG
  safe_plot_section(gsea_kegg, paste(ct, " - KEGG"))
  # Section: Reactome
  safe_plot_section(gsea_reactome, paste(ct, " - Reactome"))
  
  dev.off()
  message("Finished cell_type: ", ct, " — PDF saved.")
}

# --- save combined CSVs for all cell types ---
write.csv(cc_all, file.path(output_dir, paste0(out_prefix, "_v2_GO_CC_all_celltypes_date.csv")), row.names = FALSE)
write.csv(mf_all, file.path(output_dir, paste0(out_prefix, "_v2_GO_MF_all_celltypes_date.csv")), row.names = FALSE)
write.csv(bp_all, file.path(output_dir, paste0(out_prefix, "_v2_GO_BP_all_celltypes_date.csv")), row.names = FALSE)
write.csv(kegg_all, file.path(output_dir, paste0(out_prefix, "_v2_KEGG_all_celltypes_date.csv")), row.names = FALSE)
write.csv(reactome_all, file.path(output_dir, paste0(out_prefix, "_v2_Reactome_all_celltypes_date.csv")), row.names = FALSE)

message("\nAll done. Results and PDFs are in: ", normalizePath(output_dir))
