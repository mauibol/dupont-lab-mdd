## Differential Expression MDD vs CTRL on leiden clusters
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

################# I CUSTOMIZED THE LIBD FUNCTION ###############################
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
multi <- readRDS(paste0(MERGED_RDS_PATH, "multimodal_bpcells.rds"))
meta <- multi@meta.data
meta$pseudo_group <- paste(meta$cell_annotation_leiden, meta$sample, sep='_')

# Create single cell experiment with counts
counts <- multi[['RNA']]@layers$counts
counts <- as(counts, "dgCMatrix")
rownames(counts) <- rownames(imgns)
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = meta)
dim(counts(sce)) #36601 x 495037

# Remove genes expressed in less than 0.1% of cells
sce <- sce[rowSums(counts(sce) > 0) >= 495, ]
dim(counts(sce)) #27510 x 49503

#Pseudobulk the RNA
sce_pseudo <- aggregateAcrossCells(sce, 
                                   DataFrame(
                                     cluster = sce$cell_annotation_leiden,
                                     sample = sce$sample
                                   ))

# create the rna pseudobulk object
sce_pseudo_5 <- sce_pseudo[,sce_pseudo$ncells >= 5]

# bpcells pseudobulk method for ATAC data
atac_counts <- open_matrix_dir(dir=paste0(MERGED_RDS_PATH, 'merged/bpcells_atac'))
atac_counts <- atac_counts[, colnames(multi)]
atac_pseudo <- pseudobulk_matrix(atac_counts, cell_groups = meta$pseudo_group, method='sum')

# create sce obj from ATAC
atac_pseudo <- SingleCellExperiment(assays = list(counts = atac_pseudo))
dim(counts(atac_pseudo))

# Make sure the RNA pseudobulked object has rownames -- we will be using this to filter our atac pseudobulk
colData(sce_pseudo_5)$pseudo_group <- paste(colData(sce_pseudo_5)$cluster, colData(sce_pseudo_5)$sample, sep='_')
rownames(colData(sce_pseudo_5)) <- colData(sce_pseudo_5)$pseudo_group

# Filter based on the RNA pseudobulked object metadata
atac_pseudo_5 <- atac_pseudo[,colnames(sce_pseudo_5)]
dim(counts(atac_pseudo_5))

colData(atac_pseudo_5) <- colData(sce_pseudo_5)

# Format the pseudobulk metadata!!
atac_pseudo_5$RIN <- as.numeric(atac_pseudo_5$RIN)
atac_pseudo_5$PMI <- as.numeric(atac_pseudo_5$PMI)
atac_pseudo_5$ncells <- as.numeric(atac_pseudo_5$ncells)

colData(atac_pseudo_5)$RIN_new <- dplyr::case_when(
  colData(atac_pseudo_5)$RIN < 5.95 ~ 5.95,
  .default =colData(atac_pseudo_5)$RIN
)

colData(atac_pseudo_5)$PMI_new <- dplyr::case_when(
  colData(atac_pseudo_5)$PMI > 32.5 ~ 32.5,
  .default = colData(atac_pseudo_5)$PMI
)

colData(atac_pseudo_5) <- colData(atac_pseudo_5)[, c("sample","Diagnosis","Age",'ncells','ID4','batch','Suicide','RIN_new','PMI_new', 'Sex', 'cluster')]
colData(atac_pseudo_5)$batch <- as.factor(colData(atac_pseudo_5)$batch)
colData(atac_pseudo_5)$ID4 <- as.factor(colData(atac_pseudo_5)$ID4)

colData(atac_pseudo_5)$Diagnosis <- factor(colData(atac_pseudo_5)$Diagnosis, levels=c('NONE', 'MDD'))
saveRDS(atac_pseudo_5, paste0(METADATA_PATH, 'atac_pseudo_5_leiden.rds'))

cls <- unique(atac_pseudo_5$cluster)
cls


# Loop through all the clusters and run differential expression limma voom pipeline
dif_chromatin <- data.frame()
for (type in cls) {
  print(type)
  
  # if (type %in% unique(dif_chromatin$cell_type)) {
  #   next
  # }
  
  pse <- atac_pseudo_5[, atac_pseudo_5$cluster == type]
  meta <- colData(pse)[, c("sample","Diagnosis","Age",'ncells','ID4','batch','Suicide','RIN_w','PMI_w', 'Sex', 'cluster')]
  
  
  if (length(unique(meta$Diagnosis)) ==1) {
    print(paste('Skipping', type, 'because only 1 diagnosis in this observation'))
    next
  }
  
  # meta <- meta %>% column_to_rownames(var = 'sample_id')
  dge <- DGEList(counts = assays(pse)$counts, group = meta$Diagnosis)
  print(dim(dge))
  keep.exprs <- filterByExpr(dge, group='group')
  dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
  genes <- rownames(dge)
  dge <- calcNormFactors(dge)
  print(dim(dge))
  
  # Model design with conditional batch handling
  if (length(unique(meta$Sex)) == 1) {
    design <- model.matrix(~ Diagnosis + batch + scale(Age) + RIN_w + PMI_w + ncells, meta)
  } else {
    design <- model.matrix(~ Diagnosis + batch + Sex + scale(Age) + RIN_w + PMI_w + ncells, meta)
  }
  
  nest <- nonEstimable(design)
  if (!is.null(nest)) {
    print("Not Estimable:")
    print(nest)
    if (any(grepl("^batch", nest))) {
      if (length(unique(meta$Sex)) == 1) {
        design <- model.matrix(~ Diagnosis + scale(Age) + RIN_w + PMI_w + ncells, meta)
      } else {
        design <- model.matrix(~ Diagnosis + Sex + scale(Age) + RIN_w + PMI_w + ncells, meta)
      }
    }
  }
  # weights estimated without the correlation
  vobj <- voomWithQualityWeights(dge, design, plot = FALSE)
  dupcor <- duplicateCorrelation(vobj, design, block = meta$ID4)
  print(dupcor$consensus.correlation)
  
  # this step uses the genome-wide average for random effect
  # Fit model
  print('fitting linear model...')
  fitDupCor <- tryCatch({
    fit <- lmFit(vobj, design, block = meta$ID4, correlation = dupcor$consensus)
    fit <- eBayes(fit)
    topTable_dupcor <- topTable(fit, number = length(genes), coef = 'DiagnosisMDD')
    topTable_dupcor$cell_type <- type
    topTable_dupcor$gene <- rownames(topTable_dupcor)
    dif_chromatin <- rbind(dif_chromatin, topTable_dupcor)
  }, error = function(e) {
    message("Model failed for ", type, ": ", e$message)
    return(NULL)
  })
}
write.csv(dif_chromatin,'/data/Share/pengmad/multiome_official/mdd_vs_ctrl_degs/dars_unsupervised_batch_as_fixed_09.04.25.csv', row.names = T)


# CODE FOR DOT PLOTS
chromatin_sig <- dif_chromatin %>% dplyr::filter(adj.P.Val < 0.05 & abs(logFC) >= 0.1)
chromatin_sig$cell_type <- factor(chromatin_sig$cell_type, levels = rev(unique(chromatin_sig$cell_type)))

chromatin_sig %>%
  ggplot(aes(x = logFC,y = cell_type, color  = adj.P.Val), alpha = 0.5)+
  geom_point()+
  scale_color_distiller(palette = "Blues")+
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  ggtitle(label = 'DARs MDD vs. Ctrl')

write.csv(dif_chromatin, paste0(METADATA_PATH, 'diff_chromatin_mdd_vs_ctrl_leiden_clusters.csv'), row.names = T)


