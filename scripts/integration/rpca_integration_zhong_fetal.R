# Rpca Integration with zhong et al fetal data and neurogenic clusters

library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(SeuratDisk)
library(Signac)
library(tidyverse)
library(scCustomize)
library(Matrix)
library(BPCells)
library(dplyr)
library(harmony)
library(viridis)
library(future)
library(openxlsx)
library(rio)


# Read in neurogenic clusters subset
rna_sub <- readRDS(paste0(MERGED_RDS_PATH, 'neuro_sub_no_oligo_10.10.25.rds'))
rna_sub <- subset(rna_sub, subset = annotations %in% c('NSC.a','NSC.b','NB','INP','ImGC.1','ImGC.2'))

rna_sub <- NormalizeData(rna_sub)
rna_sub <- FindVariableFeatures(rna_sub)
rna_sub <- ScaleData(rna_sub)
rna_sub <- RunPCA(rna_sub, npcs = 40)

# Must run SCTransform for rpca integration
rna_sub <- SCTransform(rna_sub, vst.flavor = "v2", vars.to.regress = NULL, verbose = FALSE)

# Read in zhong fetal data
zhong <- readRDS('/data/Share/marimad/preprocessing_home/lasso/zhong_for_merge.rds')

# Normalize, dim reduction, cluster, etc.
zhong <- NormalizeData(zhong)
zhong <- FindVariableFeatures(zhong, nfeatures = 2000)
zhong <- ScaleData(zhong)
zhong <- RunPCA(zhong, npcs = 40)
zhong <- RunUMAP(zhong, reduction = "pca", dims = 1:20, reduction.name = 'umap', reduction.key = 'UMAP_')

ElbowPlot(zhong)

zhong <- FindNeighbors(zhong, reduction = "pca", dims = 1:20)
zhong <- FindClusters(zhong, resolution = 1, cluster.name = 'zhong1_20pcs', algorithm=4)

# Plot cannonical markers for cell type annotation
DimPlot_scCustom(zhong, group.by = "zhong1_20pcs", label = TRUE, repel = T, raster=F)
zhong_markers <- c('ASCL1', 'NEUROD2', 'GAD1', 'GAD2','OLIG2','MBP','AQP4', 'SPARC', 'PTPRC', 'RELN', 'CALB2', 'CALB1','DCX', 'EOMES', 'MEIS2', 'PROX1', 
                   'PDGFRA', 'DLX1', 'DLX2', 'SOX2', 'PAX6', 'HOPX','NES', 'GFAP', 'HES5',
                   'VIM','PCNA', 'MCM2', 'TOP2A', 'MKI67', 'NR2E1', 'NEUROD1','NEUROG2','FST')

DotPlot(zhong, assay='RNA', features = rev(zhong_markers), cols = c("yellow","blue"), group.by ='zhong1_20pcs')+ coord_flip() +
  theme(axis.text.x = element_text(angle = 90))

# Subset to only astrocytes, GCs, and progenitors
zhong_sub <- subset(zhong, subset = zhong1_20pcs %in% c(1,4,5,11,13,15,8,14))

zhong_sub <- NormalizeData(zhong_sub)
zhong_sub <- FindVariableFeatures(zhong_sub)
zhong_sub <- ScaleData(zhong_sub)
zhong_sub <- RunPCA(zhong_sub, npcs = 40)
zhong_sub <- SCTransform(zhong_sub, vst.flavor = "v2", vars.to.regress = NULL, verbose = FALSE)

zhong_sub <- FindNeighbors(zhong_sub, reduction = "pca", dims = 1:20)
zhong_sub <- FindClusters(zhong_sub, resolution = 0.8, cluster.name = 'zhong_sub_clusters', algorithm=4)

zhong_sub$fetal_cls <- dplyr::case_when(
  zhong_sub$zhong_sub_clusters == 5 ~ 'Zhong_Astro',
  zhong_sub$zhong_sub_clusters %in% c(2,4,6) ~ 'Zhong_GC',
  zhong_sub$zhong_sub_clusters == 3 ~ 'Zhong_P2',
  zhong_sub$zhong_sub_clusters == 9 ~ 'Zhong_P3',
  zhong_sub$zhong_sub_clusters == 1 ~ 'Zhong_OPC1',
  zhong_sub$zhong_sub_clusters == 8 ~ 'Zhong_OPC2',
  zhong_sub$zhong_sub_clusters == 7 ~ 'Zhong_P1',
  .default = as.character(zhong_sub$zhong_sub_clusters)
)

p_order <- c('Zhong_P1','Zhong_P2', 'Zhong_P3','Zhong_OPC1',
             'Zhong_OPC2', 'Zhong_GC', 'Zhong_Astro')
zhong_sub$fetal_cls <- factor(zhong_sub$fetal_cls, levels = rev(p_order))


immature_genes <- c('SOX2','PAX6','NES','ASCL1','HES6','SOX11','SOX4','HEY2',
                    'LMX1A','OTX2','ESM1','RSPO2','HSPB8','FABP7','FABP5','TUBA1A',
                    'NEUROD6','STMN1','STMN2','NNAT','TUBB2A','CCK','CALB2','DCX',"RELN",'NELL1','ELMO1','EPHA3','ST8SIA2',
                    'NEUROD1','CAMK2N2','GPR85','INSM1', 'YPEL4',
                    'CALB1','NEUROD2',
                    'BHLHE22','FST','STC1','TRPC6','PCOLCE2','PROX1','RBFOX3','SEMA5A')

DotPlot(zhong_sub, assay='RNA', features = zhong_markers, cols = c("yellow","blue"), group.by ='fetal_cls')+
  theme(axis.text.x = element_text(angle = 90))

DimPlot(zhong_sub, group.by = "fetal_cls", label = TRUE, repel = T, raster=F, reduction='umap')


# Select Integration features and prep integration
features <- SelectIntegrationFeatures(object.list = c(zhong_sub, rna_sub), nfeatures = 3000)
features <- unique(c(features, imm_genes, enrichment_markers))

objs <- PrepSCTIntegration(object.list = c(zhong_sub, rna_sub), anchor.features = features)

# Find integration anchors and integrate data
rpca_anchors <- FindIntegrationAnchors(object.list = objs, anchor.features = features, 
                                       reduction = "rpca", k.anchor = 20 , dims=1:30, normalization.method = 'SCT')
rpca_combined <- IntegrateData(anchorset = rpca_anchors)

# Scale and dimension reduction combined object
rpca_combined <- ScaleData(rpca_combined, verbose = FALSE)
rpca_combined <- RunPCA(rpca_combined, npcs = 30, verbose = FALSE)

ElbowPlot(rpca_combined)

rpca_combined <- RunUMAP(rpca_combined, reduction = "pca", dims = 1:10)
rpca_combined <- FindNeighbors(rpca_combined, reduction = "pca", dims = 1:10)
rpca_combined <- FindClusters(rpca_combined, resolution = 1.0, algorithm = 4)


rpca_combined$investigator <- dplyr::case_when(
  is.na(rpca_combined$investigator) ~ 'Boldrini',
  .default = as.character(rpca_combined$investigator)
)

rpca_combined$merged_annotations <- dplyr::case_when(
  is.na(rpca_combined$fetal_cls) ~ as.character(rpca_combined$annotations),
  .default = as.character(rpca_combined$fetal_cls)
)

saveRDS(rpca_combined, '/data/Share/pengmad/quest_for_immature/zhong_bold_rpca2_10.03.25.rds')
saveRDS(zhong_sub, '/data/Share/pengmad/quest_for_immature/zhong_apgs_10.03.25.rds')


# Some plots
DimPlot(rpca_combined, group.by = "merged_annotations", label = TRUE, repel = T, raster=F, reduction='umap')
DimPlot(rpca_combined, group.by = "fetal_cls", label = TRUE, repel = T, raster=F, reduction='umap')
DimPlot(rpca_combined, group.by = "annotations", label = TRUE, repel = T, raster=F, reduction='umap') + ggtitle('boldrini_subset')
DimPlot(rpca_combined, group.by = "investigator", label = TRUE, repel = T, raster=F, reduction='umap')
DimPlot(rpca_combined, group.by = "fetal_cls", label = TRUE, repel = T, raster=F, reduction='umap')

VlnPlot(zhong_sub, features=c('ASCL1','PCNA','MCM2','MKI67','TOP2A','GFAP','AQP4','HOPX','EOMES','NEUROG2','NEUROD1','NEUROD2','TBR1','SOX10','PDGFRA','MBP','PROX1','CALB2'),stack=T, pt.size=F,group.by = 'fetal_cls')


## dendrogram to see how clusters are related
Idents(rpca_combined) <- 'merged_annotations'
embeddings <- Embeddings(object = rpca_combined, reduction = 'pca')
data.dims <- lapply(
  X = levels(rpca_combined),
  FUN = function(x) {
    cells <- WhichCells(object = rpca_combined, expression=merged_annotations ==x)
    if (length(x = cells) == 1) {
      cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
  }
)
data.dims <- do.call(what = 'cbind', args = data.dims)
colnames(x = data.dims) <- levels(x = rpca_combined)
data.dist <- dist(x = t(x = data.dims))
data.tree <- ape::as.phylo(x = hclust(d = data.dist, method='complete'))
ape::plot.phylo(x = data.tree, direction = 'downwards')
ape::nodelabels()
