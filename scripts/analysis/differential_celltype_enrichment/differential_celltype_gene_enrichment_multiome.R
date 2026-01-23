## Differential Enrichment on unsupervised leiden clusters
## Using LIBD packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(spatialLIBD)
  library(here)
  library(edgeR)
  library(scuttle)
  library(scater)
  library(scran)
  library(dplyr)
  library(gridExtra)
  library(ggforce)
  library(pheatmap)
  library(scater)
  library(scran)
  library(readxl)
})

################# CUSTOMIZED THE LIBD FUNCTION ###############################
registration_stats_enrichment_custom <-
  function(
    sce_pseudo,
    block_cor,
    covars = NULL,
    var_registration = "registration_variable",
    var_sample_id = "registration_sample_id",
    gene_ensembl = NULL,
    gene_name = NULL) {
    ## For each cluster, test it against the rest
    cluster_idx <- split(seq(along = sce_pseudo[[var_registration]]), sce_pseudo[[var_registration]])
    print('I am in my function')
    message(Sys.time(), " computing enrichment statistics")
    eb0_list_cluster <- lapply(cluster_idx, function(x) {
      res <- rep(0, ncol(sce_pseudo))
      res[x] <- 1
      if (!is.null(covars)) {
        res_formula <-
          eval(str2expression(paste(
            "~", "res", "+", paste(covars, collapse = " + ")
          )))
      } else {
        res_formula <- eval(str2expression(paste("~", "res")))
      }
      m <- model.matrix(res_formula, data = colData(sce_pseudo))
      
      if (is.finite(block_cor)) {
        res <- limma::eBayes(limma::lmFit(
          logcounts(sce_pseudo),
          design = m,
          block = sce_pseudo[[var_sample_id]],
          correlation = block_cor
        ))
      } else {
        res <- limma::eBayes(limma::lmFit(logcounts(sce_pseudo),
                                          design = m
        ))
      }
      return(res)
    })
    
    return(eb0_list_cluster)
  }
#############################################################################################

# load object
multi <- LoadSeuratRds('/data/Share/pengmad/multiome_official/current_working_multimodal_bpcells.rds')
meta <- multi@meta.data

# Create single cell experiment with counts
counts <- multi[['RNA']]@layers$counts
counts <- as(counts, "dgCMatrix")
rownames(counts) <- rownames(multi)
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = meta)
dim(counts(sce)) #36601 x 495037

# Remove genes expressed in less than 0.1% of cells
sce <- sce[rowSums(counts(sce) > 0) >= 495, ]
dim(counts(sce)) #27510 x 49503

#Pseudobulk
sce_pseudo <- aggregateAcrossCells(sce, 
                                   DataFrame(
                                     cluster = sce$cell_annotation_leiden,
                                     sample = sce$sample
                                   ))

# Format the pseudobulk metadata!!
sce_pseudo$RIN <- as.numeric(sce_pseudo$RIN)
sce_pseudo$PMI <- as.numeric(sce_pseudo$PMI)
sce_pseudo$ncells <- as.numeric(sce_pseudo$ncells)


# Winsorize RIN
colData(sce_pseudo)$RIN_new <- dplyr::case_when(
  colData(sce_pseudo)$RIN < 5.95 ~ 5.95,
  .default =colData(sce_pseudo)$RIN
)


# Winsorize PMI
colData(sce_pseudo)$PMI_new <- dplyr::case_when(
  colData(sce_pseudo)$PMI > 32.5 ~ 32.5,
  .default = colData(sce_pseudo)$PMI
)

colData(sce_pseudo) <- colData(sce_pseudo)[, c("sample","Diagnosis","Age",'ncells','ID4','batch','Suicide','RIN_new','PMI_new', 'Sex', 'cluster')]
colData(sce_pseudo)$batch <- as.factor(colData(sce_pseudo)$batch)
colData(sce_pseudo)$ID4 <- as.factor(colData(sce_pseudo)$ID4)

# Scale age and define superfine.cell.class
sce_pseudo$age_scaled<-scales::rescale(sce_pseudo$Age,to=c(0,1))
sce_pseudo$superfine.cell.class<-factor(make.names(sce_pseudo$cluster))

# log normalize
sce_pseudo <- logNormCounts(sce_pseudo)

rowData(sce_pseudo)$gene_name <- rownames(sce_pseudo)
rowData(sce_pseudo)$ensemble_id <- rep('none', times=length(rownames(sce_pseudo)))

# Call LIBD functions on our pseudobulk
mod_rna<-registration_model(
  sce_pseudo,
  covars = c('batch','Sex','age_scaled', 'ncells', 'RIN_new', 'PMI_new'),
  var_registration = "superfine.cell.class"
)

cors_rna<-registration_block_cor(
  sce_pseudo,
  mod_rna,
  var_sample_id = "ID4"
)

reg<-registration_stats_enrichment_custom(
  sce_pseudo,
  block_cor=cors_rna,
  covars = c('batch','Sex','age_scaled','ncells', 'RIN_new', 'PMI_new'),
  var_registration = "superfine.cell.class",
  var_sample_id = "ID4",
  gene_ensembl = 'ensemble_id',
  gene_name = 'gene_name'
)

results_pairwise_rna <- registration_stats_pairwise(
  sce_pseudo,
  mod_rna,
  block_cor = cors_rna,
  var_registration = "superfine.cell.class",
  var_sample_id = "ID4",
  gene_ensembl = 'ensemble_id',
  gene_name = 'gene_name'
)


# Create top table for each cell type and merge results
dif_enrich_leiden <- data.frame()
for (type in names(reg)) {
  fit <- reg[[type]]
  topTable_clusters <- topTable(fit, number = length(fit$F.p.value), coef = 'res')
  topTable_clusters$cell_type <- type
  topTable_clusters$gene <- rownames(topTable_clusters)
  
  table <- topTable_clusters %>% dplyr::filter(adj.P.Val< 1)
  
  dif_enrich_leiden <- rbind(dif_enrich_leiden, topTable_clusters)
}
dif_enrich_leiden_sig <- dif_enrich_leiden %>% dplyr::filter(adj.P.Val<0.05)


# Save as a tab-separated excel
df_sort <- dif_enrich_leiden_sig %>% arrange(adj.P.Val, desc(logFC),.by_group = T)
df_list <- split(df_sort, as.factor(df_sort$cell_type))
export(df_list, 'leiden_celltype_enrichment.xlsx')

# Filter out non coding genes
protein_only <- dif_enrich_leiden_sig %>% 
  dplyr::filter(!grepl("^AC\\d{3,}", gene) & !grepl("^LINC\\d{3,}", gene) & !grepl("^AL\\d{3,}", gene)& !grepl("^BX\\d{3,}", gene)) %>% 
  dplyr::filter(!grepl("^AP\\d{3,}", gene)) %>%
  dplyr::filter(!grepl("LNC", gene)) %>% 
  dplyr::filter(!grepl("^AJ\\d{3,}", gene)) %>% dplyr::filter(!grepl("-IT\\d+$", gene) & !grepl("-AS\\d+$", gene) 
                                                              & !grepl("-DT$", gene) & !grepl("^Z\\d{3,}", gene))
