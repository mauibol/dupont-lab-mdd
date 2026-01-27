#Gene Set Enrichment Analysis on DEGs
# Maddy Peng

library(org.Hs.eg.db)
library(clusterProfiler)

degs <- read.csv('/data/Share/pengmad/multiome_official/mdd_vs_ctrl_degs/unsupervised_degs_batch_as_fixed_062525.csv', row.names = 1)
cell_types <- unique(degs$cell_type)

cc <- data.frame()
mf <- data.frame()
bp <- data.frame()

for (c in cell_types) {
  
  print(paste("Doing GSEA for", c))
  
  
  de <- degs %>% dplyr::filter(cell_type == c) 
  de$ranking_score <- sign(de$logFC)*(-log10(de$adj.P.Val))
  
  de <- de %>% dplyr::arrange(desc(ranking_score))
  
  if (nrow(de)<3) {
    next
  }
  
  entrez <- bitr(de$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  print(nrow(entrez))
  entrez <- entrez %>% distinct(SYMBOL, .keep_all = TRUE)
  
  de_genes <- de[de$gene %in% entrez$SYMBOL,]
  stopifnot(de_genes$gene == entrez$SYMBOL)
  de_genes$ENTREZID <- entrez$ENTREZID
  
  entrez <- as.vector(de_genes$ranking_score)
  names(entrez) <- de_genes$ENTREZID
  
  genes_cc <- gseGO(geneList = entrez,
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    minGSSize = 15,
                    maxGSSize = 400,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
  
  genes_cc <- setReadable(genes_cc, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  genes_mf <- gseGO(geneList = entrez,
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    minGSSize = 15,
                    maxGSSize = 400,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
  
  genes_mf <- setReadable(genes_mf, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  genes_bp <- gseGO(geneList = entrez,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 15,
                    maxGSSize = 400,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
  
  genes_bp <- setReadable(genes_bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  print('saving...')
  
  if (length(genes_cc@result$ID > 0)){
    genes_cc@result$cell_type <- c
    cc <- rbind(cc, genes_cc@result)
  }
  
  
  if (length(genes_mf@result$ID > 0)){
    genes_mf@result$cell_type <- c
    mf <- rbind(mf, genes_mf@result)
  }
  
  
  if (length(genes_bp@result$ID > 0)){
    genes_bp@result$cell_type <- c
    bp <- rbind(bp, genes_bp@result)
  }
  
  print('plotting')
  pdf(paste0('/data/Share/pengmad/multiome_official/GSEA_on_degs/', c, '_gsea_10.24.25.pdf'), width = 10, height = 10, onefile = T)
  
  direction_colors <- c("Upregulated" = "#E64B35FF", "Downregulated" = "#4DBBD5FF")
  
  if (nrow(genes_cc@result) > 0) {
    genes_cc@result$Direction <- ifelse(genes_cc@result$NES > 0, "Upregulated", "Downregulated")
    p1 <- dotplot(genes_cc, x = "NES", showCategory = 50, label_format = 60) +
      ggtitle(paste(c, "GSEA for CC"))
    print(p1)
  }
  
  if (nrow(genes_mf@result) > 0) {
    genes_mf@result$Direction <- ifelse(genes_mf@result$NES > 0, "Upregulated", "Downregulated")
    p2 <- dotplot(genes_mf, x = "NES", showCategory = 50, label_format = 60) +
      ggtitle(paste(c, "GSEA for MF"))
    print(p2)
  }
  
  if (nrow(genes_bp@result) > 0) {
    genes_bp@result$Direction <- ifelse(genes_bp@result$NES > 0, "Upregulated", "Downregulated")
    p3 <- dotplot(genes_bp, x = "NES", showCategory = 50, label_format = 60) +
      ggtitle(paste(c, "GSEA for BP"))
    print(p3)
  }
  
  dev.off()
  
}
dev.off()

# Save the excels
save_excel <- function(df, path){
  df_sort <- df %>% dplyr::arrange(p.adjust, .by_group = T)
  df_list <- split(df_sort, as.factor(df_sort$cell_type))
  export(df_list, path)
}

save_excel(cc, '/data/Share/pengmad/multiome_official/GSEA_on_degs/cc_10.24.25.xlsx')
save_excel(mf, '/data/Share/pengmad/multiome_official/GSEA_on_degs/mf_10.24.25.xlsx')
save_excel(bp, '/data/Share/pengmad/multiome_official/GSEA_on_degs/bp_10.24.25.xlsx')

