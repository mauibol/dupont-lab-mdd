#############################################
# Visium v1 and v2 Differential Gene Enrichment
# By Victor Anosike
# Adapted from code By Maddy Peng
#############################################

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
  library(Matrix)
  library(rio)
  library(openxlsx)
  
})

############################# Visium v1 ##################################################################
##########################################################################################################

##CUSTOMIZED LIBD FUNCTION ##
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
#############################################################################################

# Differential Enrichment on the subset annotations official
# PERFORM DIFFERENTIAL ENRICHMENT ON THESE SUBCLUSTERS

####Retrieve combined, integrated and processed v1 clusters
v1 <- readRDS('insert/path/to/v1.rds')

#Add diagnosis and sex columns to object
colnames(v1@meta.data)[colnames(v1@meta.data) == "Psych1A"] <- "diagnosis"
v1@meta.data$sex <- rep("m", 65133) 

# - All CM samples for V1 belong to one batch, and all MM samples to another batch. That means that we can use
# - batch as a covariable, but since there are only two experimenters, the experimenter column in the V1 metadata effectively acts as 
# - a batch covariable for V1. However, this is not the case for V2, since V2 has six different batches. Therefore, V1 will simply 
# - have the experimenter column's NAME replaced to "batch", but V2 will have the experimenter column removed as a covariable ENTIRELY 
# - and replaced with the variable batch. 
v1$batch <- as.factor(as.character(substr(v1$sample, 1, 2)))
# The above line is there for V1 ONLY because V1 does not have batch info added to the metadata naturally


counts <- v1@assays[["Spatial"]]@data
counts <- as(counts, "dgCMatrix")
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = v1@meta.data)
dim(counts(sce)) #2055 x 190282

#Pseudobulk
sce_pseudo <- aggregateAcrossCells(sce, 
                                   DataFrame(
                                     cluster = v1@meta.data[["seurat_clusters"]],
                                     sample = sce$sample
                                   ))

# Format the pseudobulk metadata!!
sce_pseudo$RIN <- as.numeric(sce_pseudo$RIN)
sce_pseudo$PMI <- as.numeric(sce_pseudo$PMI)
sce_pseudo$ncells <- as.numeric(sce_pseudo$ncells)

colData(sce_pseudo)$RIN_new <- dplyr::case_when(
  colData(sce_pseudo)$RIN < 5.95 ~ 5.95,
  .default =colData(sce_pseudo)$RIN
)

colData(sce_pseudo)$PMI_new <- dplyr::case_when(
  colData(sce_pseudo)$PMI > 32.5 ~ 32.5,
  .default = colData(sce_pseudo)$PMI
)


colData(sce_pseudo) <- colData(sce_pseudo)[, c("sample","diagnosis","age",'ncells','donor','batch','suicide','RIN_new','PMI_new', 'cluster')]  
#For the above line, 'sex' was removed as a coavariate because all of the v1 samples are from male subjects


colData(sce_pseudo)$batch <- as.factor(colData(sce_pseudo)$batch)
colData(sce_pseudo)$donor <- as.factor(colData(sce_pseudo)$donor)

# Scale age and define superfine.cell.class
sce_pseudo$age_scaled<-scales::rescale(sce_pseudo$age,to=c(0,1))
sce_pseudo$superfine.cell.class<-factor(make.names(sce_pseudo$cluster))

# log normalize
sce_pseudo <- logNormCounts(sce_pseudo)

rowData(sce_pseudo)$gene_name <- rownames(sce_pseudo)
rowData(sce_pseudo)$ensemble_id <- rep('none', times=length(rownames(sce_pseudo)))

# Call LIBD functions on our pseudobulk
mod<-registration_model(
  sce_pseudo,
  covars = c('batch','age_scaled', 'ncells', 'RIN_new', 'PMI_new'),
  var_registration = "superfine.cell.class"
)

cors<-registration_block_cor(
  sce_pseudo,
  mod,
  var_sample_id = "donor"
)

reg<-registration_stats_enrichment_custom(
  sce_pseudo,
  block_cor=cors,
  covars = c('batch','age_scaled','ncells', 'RIN_new', 'PMI_new'),
  var_registration = "superfine.cell.class",
  var_sample_id = "donor",
  gene_ensembl = 'ensemble_id',
  gene_name = 'gene_name'
)

# Create top table for each cell type and merge results
de_imgn <- data.frame()
for (type in names(reg)) {
  fit <- reg[[type]]
  topTable_clusters <- topTable(fit, number = length(fit$F.p.value), coef = 'res')
  topTable_clusters$cell_type <- type
  topTable_clusters$gene <- rownames(topTable_clusters)
  
  table <- topTable_clusters %>% dplyr::filter(adj.P.Val< 1)
  
  de_imgn <- rbind(de_imgn, topTable_clusters)
}

#Export de_imgn (which is not limited by P-Value) to a csv file
write.csv(de_imgn, '/path/to/v1_celltype_enrichment_sig_and_nonsig_date_batch_effect.csv', row.names = T)

#Export de_imgn_sig (which only contains pvalues less than 0.05) to a csv file
de_imgn_sig <- de_imgn %>% dplyr::filter(adj.P.Val<0.05, logFC>0)
write.csv(de_imgn_sig,'/path/to/v1_celltype_enrichment_date_batch_effect.csv', row.names = T)

# Create excel sheet to be used in visium_diff_enrichment_ncbi_and_uniprot_naming.R and visium_diff_enrichment_kegg_reactome_naming.py
df_sort <- de_imgn_sig %>% dplyr::arrange(adj.P.Val, desc(logFC),.by_group = T)
df_list <- split(df_sort, as.factor(df_sort$cell_type))
write.xlsx(df_list,
           file = "/path/to/v1_differential_enrichment_date_batch_effect.xlsx",
           asTable = TRUE)


#Save Gene Enrichment Data for Later
saveRDS(df_list, "/path/to/v1_gene_enrichment_data_date_batch_effect.rds")



########################## Visium v2 ############################################################
#################################################################################################

###CUSTOMIZED LIBD FUNCTION ####
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
#############################################################################################

# Differential Enrichment on the subset annotations official
# PERFORM DIFFERENTIAL ENRICHMENT ON THESE SUBCLUSTERS

####Retrieve combined, integrated and processed v2 clusters
v2 <- readRDS('insert/path/to/v2.rds')

# - All CM samples for V1 belong to one batch, and all MM samples to another batch. That means that we can use
# - batch as a covariable, but since there are only two experimenters, the experimenter column in the V1 metadata effectively acts as 
# - a batch covariable for V1. However, this is not the case for V2, since V2 has six different batches. Therefore, V1 will simply 
# - have the experimenter column's NAME replaced to "batch", but V2 will have the experimenter column removed as a covariable ENTIRELY 
# - and replaced with the variable batch, that already exists in V2. As a result, the below line is no longer needed since v2.rds has a batch factor included already
#combined$experimenter <- as.factor(as.character(substr(combined$sample, 1, 2)))


#Make sure gene names are preserved
rownames(v2@assays[["Spatial"]]@layers[["counts.1"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.2"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.3"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.4"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.5"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.6"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.7"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.8"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.9"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.10"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.11"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.12"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.13"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.14"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.15"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.16"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.17"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.18"]]) <- rownames(v2@assays[["Spatial"]])

#Create counts
counts <- cbind(v2@assays[["Spatial"]]@layers[["counts.1"]],
                v2@assays[["Spatial"]]@layers[["counts.2"]],
                v2@assays[["Spatial"]]@layers[["counts.3"]],
                v2@assays[["Spatial"]]@layers[["counts.4"]],
                v2@assays[["Spatial"]]@layers[["counts.5"]],
                v2@assays[["Spatial"]]@layers[["counts.6"]],
                v2@assays[["Spatial"]]@layers[["counts.7"]],
                v2@assays[["Spatial"]]@layers[["counts.8"]],
                v2@assays[["Spatial"]]@layers[["counts.9"]],
                v2@assays[["Spatial"]]@layers[["counts.10"]],
                v2@assays[["Spatial"]]@layers[["counts.11"]],
                v2@assays[["Spatial"]]@layers[["counts.12"]],
                v2@assays[["Spatial"]]@layers[["counts.13"]],
                v2@assays[["Spatial"]]@layers[["counts.14"]],
                v2@assays[["Spatial"]]@layers[["counts.15"]],
                v2@assays[["Spatial"]]@layers[["counts.16"]],
                v2@assays[["Spatial"]]@layers[["counts.17"]],
                v2@assays[["Spatial"]]@layers[["counts.18"]])


counts <- as(counts, "dgCMatrix")
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = v2@meta.data)
dim(counts(sce)) #18085 x 177566


#Pseudobulk
sce_pseudo <- aggregateAcrossCells(sce, 
                                   DataFrame(
                                     cluster = v2@meta.data[["seurat_clusters"]],
                                     sample = sce$sample
                                   ))

# Format the pseudobulk metadata!!
sce_pseudo$RIN <- as.numeric(sce_pseudo$RIN)
sce_pseudo$PMI <- as.numeric(sce_pseudo$PMI)
sce_pseudo$ncells <- as.numeric(sce_pseudo$ncells)

colData(sce_pseudo)$RIN_new <- dplyr::case_when(
  colData(sce_pseudo)$RIN < 5.95 ~ 5.95,
  .default =colData(sce_pseudo)$RIN
)

colData(sce_pseudo)$PMI_new <- dplyr::case_when(
  colData(sce_pseudo)$PMI > 32.5 ~ 32.5,
  .default = colData(sce_pseudo)$PMI
)


colData(sce_pseudo) <- colData(sce_pseudo)[, c("sample","diagnosis","age",'ncells','donor','batch','suicide','RIN_new','PMI_new', 'sex', 'cluster')]


colData(sce_pseudo)$batch <- as.factor(colData(sce_pseudo)$batch)
colData(sce_pseudo)$donor <- as.factor(colData(sce_pseudo)$donor)

# Scale age and define superfine.cell.class
sce_pseudo$age_scaled<-scales::rescale(sce_pseudo$age,to=c(0,1))
sce_pseudo$superfine.cell.class<-factor(make.names(sce_pseudo$cluster))

# log normalize
sce_pseudo <- logNormCounts(sce_pseudo)

rowData(sce_pseudo)$gene_name <- rownames(sce_pseudo)
rowData(sce_pseudo)$ensemble_id <- rep('none', times=length(rownames(sce_pseudo)))

# Call LIBD functions on our pseudobulk
mod<-registration_model(
  sce_pseudo,
  covars = c('batch','sex','age_scaled', 'ncells', 'RIN_new', 'PMI_new'),
  var_registration = "superfine.cell.class"
)

cors<-registration_block_cor(
  sce_pseudo,
  mod,
  var_sample_id = "donor"
)

reg<-registration_stats_enrichment_custom(
  sce_pseudo,
  block_cor=cors,
  covars = c('batch','sex','age_scaled','ncells', 'RIN_new', 'PMI_new'),
  var_registration = "superfine.cell.class",
  var_sample_id = "donor",
  gene_ensembl = 'ensemble_id',
  gene_name = 'gene_name'
)

# Create top table for each cell type and merge results
de_imgn <- data.frame()
for (type in names(reg)) {
  fit <- reg[[type]]
  topTable_clusters <- topTable(fit, number = length(fit$F.p.value), coef = 'res')
  topTable_clusters$cell_type <- type
  topTable_clusters$gene <- rownames(topTable_clusters)
  
  table <- topTable_clusters %>% dplyr::filter(adj.P.Val< 1)
  
  de_imgn <- rbind(de_imgn, topTable_clusters)
}

#Export de_imgn (which is not limited by P-Value) to a csv file for later use
write.csv(de_imgn, '/path/to/v2_celltype_enrichment_sig_and_nonsig_date_batch_effect.csv', row.names = T)

#Export de_imgn_sig (which only contains pvalues less than 0.05) to a csv file
de_imgn_sig <- de_imgn %>% dplyr::filter(adj.P.Val<0.05, logFC>0)
write.csv(de_imgn_sig,'/path/to/v2_celltype_enrichment_date_batch_effect.csv', row.names = T)

# Create excel sheet to be used in visium_diff_enrichment_ncbi_and_uniprot_naming.R and visium_diff_enrichment_kegg_reactome_naming.py
df_sort <- de_imgn_sig %>% dplyr::arrange(adj.P.Val, desc(logFC),.by_group = T)
df_list <- split(df_sort, as.factor(df_sort$cell_type))
write.xlsx(df_list,
           file = '/path/to/v2_differential_enrichment_date_batch_effect.xlsx',
           asTable = TRUE)

#Save Gene Enrichment Data for Later
saveRDS(df_list, "/path/to/v2_gene_enrichment_data_11_6_2025_batch_effect.rds")





























