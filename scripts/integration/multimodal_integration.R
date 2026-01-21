# Multimodal integration of ATAC and RNA with BPCells
# Maddy Peng

library(Seurat)
library(SeuratObject)
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

# Read in both RNA and ATAC objects
rna <- readRDS(paste0(MERGED_RDS_PATH, "merged_rna_bpcells.rds"))
atac <- readRDS(paste0(MERGED_RDS_PATH, "merged_atac_bpcells.rds"))

# Find common cells between the atac and rna objects and subset both
common_cells <- intersect(Cells(rna), Cells(atac))
rna <- subset(rna, cells=common_cells)
atac <- subset(atac_merge, cells=common_cells)

# Add the two assays to the same object -- we will use the RNA metadata
atac@meta.data <- rna@meta.data
multi <- rna
multi[['ATAC']] <- atac[['ATAC']]

# Run Seurat Pipeline on each assay separately
#### RNA ####
DefaultAssay(multi) <- 'RNA'
multi <- NormalizeData(multi)
multi <- FindVariableFeatures(multi)
multi <- ScaleData(multi)
multi <- RunPCA(multi, npcs = 40)
multi <- RunUMAP(multi, reduction = 'pca', dims = 1:40, reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

multi <- RunHarmony(multi, group.by.vars = c('sample', 'batch', 'ID4'),
                  plot_convergence = T, reduction.save = "harmony.rna")

multi <- FindNeighbors(multi, reduction = "harmony.rna", dims = 1:40)
multi <- FindClusters(multi, resolution = 0.8, cluster.name = 'harmony_rna_0.8')
multi <- RunUMAP(multi, reduction = "harmony.rna", dims = 1:40, reduction.name = 'harmony.rna.umap', reduction.key = 'HarmonyRnaUMAP_')



# Because Signac does not have compatibility with BPCells, must run TFIDF on the ATAC matrix itself
# Normalize and LSI ATAC BPCells matrix 
mat_lsi <- mtx %>%
  multiply_cols(1 / Matrix::colSums(mtx)) %>%
  multiply_rows(1 / Matrix::rowMeans(mtx))

mat_lsi <- log1p(10000 * mat_lsi)

cell_peak_stats <- matrix_stats(mat_lsi, col_stats="variance")$col_stats
cell_means <- cell_peak_stats["mean",]
cell_vars <- cell_peak_stats["variance",]
mat_lsi_norm <- mat_lsi %>%
  add_cols(-cell_means) %>%
  multiply_cols(1 / cell_vars)

svd_atac <- BPCells::svds(mat_lsi_norm, k=40)
pca_atac <- multiply_cols(svd_atac$v, svd_atac$d)


# Add lsi embedding into Seurat Object
rownames(pca_atac) <- colnames(atac)
atac[['lsi']] <- Seurat::CreateDimReducObject(
  embeddings = pca_atac,
  key = "lsi_",
  assay = 'ATAC'
)


#UMAP and clustering PRE-integration
atac <- RunUMAP(atac, reduction = "lsi", dims = 1:30, 
                reduction.name = 'umap', 
                reduction.key = 'umap_')

atac_merge <- FindNeighbors(atac_merge, reduction='lsi', dims=1:30)
atac_merge <- FindClusters(atac_merge, reduction='lsi', resolution=0.8, cluster.name = 'atac_clusters_0.8')


#Harmony Integration
atac <- RunHarmony(atac, group.by.vars = c('sample', 'batch', 'ID4'), reduction = 'lsi', dims.use = 2:40,
                   plot_convergence = T, reduction.save = "harmony.atac", assay.use = "ATAC", verbose = TRUE, project.dim = FALSE)



atac <- RunUMAP(atac, reduction = "harmony.atac", dims = 1:29, 
                reduction.name = 'harmony.atac.umap', 
                reduction.key = 'harmony.umap_')

atac <- FindNeighbors(atac, reduction='harmony.atac', dims=1:29)
atac <- FindClusters(atac, reduction='harmony.atac', resolution=0.8, cluster.name = 'harmony_atac_0.8')

# PLOTS

#lsi
p1 <- DimPlot(atac_merge, reduction = "lsi", dims = c(1, 2), group.by=c('batch'), raster = T)
p2 <- DimPlot(atac_merge, reduction = "lsi", dims = c(1, 2), group.by=c('sample'), raster = T)
p3 <- DimPlot(atac_merge, reduction = "lsi", dims = c(1, 2), group.by=c('seq_study_ID'), raster = T)
p4 <- DimPlot(atac_merge, reduction = "harmony.atac", dims = c(1, 2), group.by=c('batch'), raster = T)
p5 <- DimPlot(atac_merge, reduction = "harmony.atac", dims = c(1, 2), group.by=c('sample'), raster = T)
p6 <- DimPlot(atac_merge, reduction = "harmony.atac", dims = c(1, 2), group.by=c('seq_study_ID'), raster = T)

#umap
p7 <- DimPlot(atac_merge, reduction = "umap", dims = c(1, 2), group.by=c('batch'), raster = T)
p8 <- DimPlot(atac_merge, reduction = "umap", dims = c(1, 2), group.by=c('sample'), raster = T)
p9 <- DimPlot(atac_merge, reduction = "umap", dims = c(1, 2), group.by=c('ID4'), raster = T)
p10 <- DimPlot(atac_merge, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('batch'), raster = T)
p11 <- DimPlot(atac_merge, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('sample'), raster = T)
p12 <- DimPlot(atac_merge, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('ID4'), raster = T)


#Clusters
cluster_rna <- DimPlot(atac_merge, reduction = "umap", dims = c(1, 2), group.by=c('harmony_atac_clusters_0.8'), raster = T, label=T)
cluster_h <- DimPlot(atac, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('predictions'), raster = T, label=T, label.size = 2, repel = T)

cluster_h

grid.arrange(cluster_rna, cluster_h, ncol=2)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)
grid.arrange(p7, p8, p9, p10, p11, p12, ncol=3, nrow=2)