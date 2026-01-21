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
############### RNA #############
DefaultAssay(multi) <- 'RNA'
multi <- NormalizeData(multi)
multi <- FindVariableFeatures(multi)
multi <- ScaleData(multi)
multi <- RunPCA(multi, npcs = 40)

multi <- RunHarmony(multi, group.by.vars = c('sample', 'batch', 'ID4'),
                  plot_convergence = T, reduction.save = "harmony.rna")

multi <- FindNeighbors(multi, reduction = "harmony.rna", dims = 1:40)
multi <- FindClusters(multi, resolution = 0.8, cluster.name = 'harmony_rna_0.8')
multi <- RunUMAP(multi, reduction = "harmony.rna", dims = 1:40, reduction.name = 'harmony.rna.umap', reduction.key = 'HarmonyRnaUMAP_')


# Run Seurat Pipeline on each assay separately
########### ATAC ##############
DefaultAssay('ATAC')

# Because Signac does not have compatibility with BPCells, must run TFIDF on the ATAC matrix itself
mtx <- multi[['ATAC']]$counts

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

# Add lsi embedding into Seurat Object atac assay
rownames(pca_atac) <- colnames(multi)
multi[['lsi']] <- Seurat::CreateDimReducObject(
  embeddings = pca_atac,
  key = "lsi_",
  assay = 'ATAC'
)

#UMAP, harmony, clustering
multi <- RunHarmony(multi, group.by.vars = c('sample', 'batch', 'ID4'), reduction = 'lsi', dims.use = 1:40,
                   plot_convergence = T, reduction.save = "harmony.atac", assay.use = "ATAC", verbose = TRUE, project.dim = FALSE)

multi <- FindNeighbors(multi, reduction='harmony.atac', dims=1:40)
multi <- FindClusters(multi, reduction='harmony.atac', resolution=0.8, cluster.name = 'harmony_atac_0.8')

multi <- RunUMAP(multi, reduction = "harmony.atac", dims = 1:40, 
                 reduction.name = 'harmony.atac.umap', 
                 reduction.key = 'HarmonyAtacUMAP_')


############### Multimodal WNN Analysis #######################
multi <- FindMultiModalNeighbors(
  multi,
  reduction.list = list('harmony.rna', 'harmony.atac'),
  dims.list = list(1:40, 1:39),
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  modality.weight.name = "RNA.weight",
  verbose = TRUE,
  prune.SNN = 1/20,
  k.nn = 20
)

multi <- RunUMAP(
  multi,
  nn.name='weighted.nn',
  reduction.name='wnn.umap'
)

# Algorithm=4 is leiden 
multi <- FindClusters(multi, graph.name = 'wsnn', algorithm = 4, 
                      verbose = T, resolution = 0.8, 
                      cluster.name='wnn_0.8', random.seed = 1)

saveRDS(multi, paste0(MERGED_RDS_PATH, 'multimodal_bpcells.rds'))

################################## PLOTS ######################################

# rna
p1 <- DimPlot(multi, reduction = "pca", dims = c(1, 2), group.by=c('batch'), raster = T)
p2 <- DimPlot(multi, reduction = "pca", dims = c(1, 2), group.by=c('sample'), raster = T)
p3 <- DimPlot(multi, reduction = "pca", dims = c(1, 2), group.by=c('ID4'), raster = T)
p4 <- DimPlot(multi, reduction = "harmony.rna", dims = c(1, 2), group.by=c('batch'), raster = T)
p5 <- DimPlot(multi, reduction = "harmony.rna", dims = c(1, 2), group.by=c('sample'), raster = T)
p6 <- DimPlot(multi, reduction = "harmony.rna", dims = c(1, 2), group.by=c('ID4'), raster = T)
p10 <- DimPlot(multi, reduction = "harmony.rna.umap", dims = c(1, 2), group.by=c('batch'), raster = T)
p11 <- DimPlot(multi, reduction = "harmony.rna.umap", dims = c(1, 2), group.by=c('sample'), raster = T)
p12 <- DimPlot(multi, reduction = "harmony.rna.umap", dims = c(1, 2), group.by=c('ID4'), raster = T)
p13 <- DimPlot(multi, reduction = "harmony.rna.umap", dims = c(1, 2), group.by=c('harmony_rna_0.8')
               , raster = T, label=T, label.size = 2, repel = T)

# atac
p1 <- DimPlot(multi, reduction = "lsi", dims = c(1, 2), group.by=c('batch'), raster = T)
p2 <- DimPlot(multi, reduction = "lsi", dims = c(1, 2), group.by=c('sample'), raster = T)
p3 <- DimPlot(multi, reduction = "lsi", dims = c(1, 2), group.by=c('ID4'), raster = T)
p4 <- DimPlot(multi, reduction = "harmony.atac", dims = c(1, 2), group.by=c('batch'), raster = T)
p5 <- DimPlot(multi, reduction = "harmony.atac", dims = c(1, 2), group.by=c('sample'), raster = T)
p6 <- DimPlot(multi, reduction = "harmony.atac", dims = c(1, 2), group.by=c('ID4'), raster = T)
p10 <- DimPlot(multi, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('batch'), raster = T)
p11 <- DimPlot(multi, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('sample'), raster = T)
p12 <- DimPlot(multi, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('ID4'), raster = T)
p13 <- DimPlot(multi, reduction = "harmony.atac.umap", dims = c(1, 2), group.by=c('harmony_atac_0.8')
               , raster = T, label=T, label.size = 2, repel = T)

# multimodal
plot_genes <- c('PDGFRA','OLIG1','MBP','MOG','MAL','PTPRC','P2RY12','APBB1IP','CD163','CD247',
           'PECAM1','FLT1','NES','CEMIP',
           'TTR','HTR2C','CFAP54','AQP4','SLC1A2','ALDH1L1',
           'GLUL','GFAP','PAX6','SOX2','MYT1L','RBFOX3',
           'SOX1','DCX','GAD1','GAD2','SLC32A1',
           'RELN','CALB2','NR2E1','PROX1',
           'CCK','SLC17A8','VIP','PVALB','LHX6','SST','PENK',
           'BHLHE22','POSTN','CALB1','SEMA5A','SLC17A7',
           'CAMK2A','SV2B','SATB2')

p1 <- DimPlot(multi, reduction = "wnn.umap", dims = c(1, 2), group.by=c('batch'), raster = T)
p2 <- DimPlot(multi, reduction = "wnn.umap", dims = c(1, 2), group.by=c('sample'), raster = T)
p3 <- DimPlot(multi, reduction = "wnn.umap", dims = c(1, 2), group.by=c('ID4'), raster = T)
p4 <- DimPlot(multi, reduction = "wnn.umap", dims = c(1, 2), group.by=c('wnn_0.8'),
              raster = T, label=T, label.size = 2, repel = T)

p5 <- DotPlot(multi, assay='RNA', features = rev(plot_genes), cols = c("yellow","blue"), group.by ='wnn_0.8' ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90))


# dendrogram to see similarity between clusters
embeddings <- Embeddings(object = multi, reduction = 'harmony.rna')
data.dims <- lapply(
  X = levels(multi),
  FUN = function(x) {
    cells <- WhichCells(object = multi, expression=wnn_0.8 ==x)
    if (length(x = cells) == 1) {
      cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
  }
)
data.dims <- do.call(what = 'cbind', args = data.dims)
colnames(x = data.dims) <- levels(x = multi)
data.dist <- dist(x = t(x = data.dims))
data.tree <- ape::as.phylo(x = hclust(d = data.dist))
ape::plot.phylo(x = data.tree, direction = 'downwards')
ape::nodelabels()


# Annotate the clusters based off gene expression dot plot
multi$cell_annotation_leiden <- dplyr::case_when(
  multi$wsnn_0.8 ==1 ~ "Oligo.1",
  multi$wsnn_0.8 ==3 ~ "Oligo.2",
  multi$wsnn_0.8 ==30 ~ "Oligo.3",
  multi$wsnn_0.8 ==10 ~ "Astro.2",
  multi$wsnn_0.8 == 6 ~ "Astro.1",
  multi$wsnn_0.8 ==2 ~ "GC.1",
  multi$wsnn_0.8 ==13 ~ "GC.2",
  multi$wsnn_0.8 ==22 ~ "GC.3",
  multi$wsnn_0.8 ==24 ~ "GC.4",
  multi$wsnn_0.8 ==29 ~ "GC.5",
  multi$wsnn_0.8 ==5 ~ "Micro",
  multi$wsnn_0.8 ==31 ~ "TC",
  multi$wsnn_0.8 ==21 ~ "Epe",
  multi$wsnn_0.8 ==9 ~ "Endo",
  multi$wsnn_0.8 ==27 ~ "VLMC",
  multi$wsnn_0.8 ==7 ~ "OPC",
  multi$wsnn_0.8 ==19 ~ "CP",
  multi$wsnn_0.8 ==14 ~ "InN.2",
  multi$wsnn_0.8 ==12 ~ "InN.1",
  multi$wsnn_0.8 ==20 ~ "InN.4",
  multi$wsnn_0.8 ==18 ~ "InN.3",
  multi$wsnn_0.8 ==26 ~ "InN.6",
  multi$wsnn_0.8 ==28 ~ "InN.7",
  multi$wsnn_0.8 ==25 ~ "InN.5",
  multi$wsnn_0.8 ==4 ~ "ExN.1",
  multi$wsnn_0.8 ==8 ~ "ExN.2",
  multi$wsnn_0.8 ==11 ~ "ExN.3",
  multi$wsnn_0.8 ==15 ~ "ExN.4",
  multi$wsnn_0.8 ==16 ~ "ExN.5",
  multi$wsnn_0.8 ==17 ~ "ExN.6",
  multi$wsnn_0.8 ==23 ~ "ExN.7",
)

# Plot proportions of cluster in samples
props <- prop.table(table(multi$sample, multi$cell_annotation_leiden)) %>% as.data.frame()
random_colors <- sample(colors(), length(unique(props$Var2)))
ggplot(props, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = 'fill') +
  labs(title = "WNN Clusters",
       x = "Cluster", y = "Proportion", fill = "Sample") +
  theme_minimal() +
  scale_fill_manual(values = random_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
