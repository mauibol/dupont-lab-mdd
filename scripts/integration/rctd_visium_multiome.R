# RCTD -- developed from Anthony's Script
# Maddy Peng

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(Matrix)
  library(doParallel)
  library(scatterpie)
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

######## V1 RCTD ###############


# load single-nucleus multiome that overlap in donors
multi <- readRDS(paste0(MERGED_RDS_PATH, 'multimodal_bpcells.rds'))
Idents(multi) <- "cell_annotation_leiden"


# Create single cell experiment
counts <- multi[['RNA']]@layers$counts
counts <- as(counts, "dgCMatrix")
rownames(counts) <- rownames(multi)
meta <- multi@meta.data
meta$nUMI <- meta$nCount_RNA


reference_sce <- SummarizedExperiment(assays = list(counts = counts), 
                                      colData = meta[,c('cell_annotation_leiden', 'nUMI', 'broad_type_leiden')])

counts_2 <- imgns_sub[['RNA']]@layers$counts
rownames(counts_2) <- rownames(imgns_sub)
meta_2 <- imgns_sub@meta.data
meta_2$nUMI <- meta_2$nCount_RNA


#load v1 obj
v1 <- readRDS('/data/Share/Victor_Anosike/V1_Seurat_Spatial/integrated_visium_V1_PRECAST_and_Seurat.rds')
v1$nUMI <- v1$nCount_Spatial

v1_sub <- subset(v1, subset = sample == 'CM014')
v1_sub <- AddMetaData(v1_sub, metadata = v1_sub@images$CM014@coordinates[,2], col.name = 'row')
v1_sub <- AddMetaData(v1_sub, metadata = v1_sub@images$CM014@coordinates[,3], col.name = 'col')

v_meta <- v1_sub@meta.data


#Create Spatial Experiment for v1 obj
spatial_spe_v1 <- SpatialExperiment(
  assay= v1_sub@assays$Spatial@counts,
  spatialCoordsNames = c('col','row'),
  colData = v_meta[,c('nUMI', 'seurat_clusters','row','col')]
)

ggplot(spatialCoords(spatial_spe_v1), aes(col, row)) + 
  geom_point(size = 5) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# Create RCTD
rctd_data <- createRctd(spatial_spe, reference_sce, cell_type_col = 'cell_annotation_leiden', 
                        require_int = F)

# Run RCTD

results_spe <- runRctd(rctd_data, rctd_mode = "multi", max_cores = 4, max_multi_types = 4)
saveRDS(results_spe, '/data/Share/pengmad/visium/rctd_results_v1.rds')

colData(results_spe)

classification_df <- data.frame(
  pixel = colnames(assay(results_spe)),
  spot_class = colData(results_spe)$spot_class,
  first_type = colData(results_spe)$first_type,
  second_type = colData(results_spe)$second_type
)

plotAllWeights(results_spe, 'weights_full', r = 1, lwd = 0, title = "Rctd Results")

cell_types <- rownames(results_spe_traj1)
palette_named <- setNames(palette1[seq_along(cell_types)], cell_types)
palette_named["reject"] <- "#BEBEBE"  # or "#A0A0A0" for a slightly darker gray


# check
length(palette_named)
tail(palette_named)


p1 <- plotAllWeights(
  results_spe,
  r = .7,
  lwd = 0,
  title = "v2 Cell Type Decomposition with Multiome Clusters"
) +
  scale_fill_manual(values = palette_named) +
  scale_color_manual(values = palette_named)

print(p1)

pdf("/data/Share/pengmad/visium/rctd_results.pdf", width = 16, height = 6)
plotAllWeights(results_spe, r = .7, lwd = 0, title = "V1 Cell Type Decomposition with Multiome Clusters")
dev.off()

SpatialDimPlot_scCustom(v1_sub, group.by = 'seurat_clusters',label.size = 5, label = T)




######################### v2 RCTD ########################
#load v2 obj
v2 <- readRDS('/data/Share/Victor_Anosike/V2_Seurat_Spatial/combined_object_integrated_processed_updated_metadata_V2.rds')
v2$nUMI <- v2$nCount_Spatial

v2_sub <- subset(v2, subset = sample == 'MM131')

coords <- v2_sub@images$MM131@boundaries$centroids@coords
rownames(coords) <- colnames(v2_sub)

coords

y_coords <- v2_sub@images$MM131@boundaries$centroids@coords[,'y']
rownames(y_coords) <- colnames(v2_sub)

v2_sub <- AddMetaData(v2_sub, metadata = coords[,'x'], col.name = 'x')
v2_sub <- AddMetaData(v2_sub, metadata = coords[,'y'], col.name = 'y')

v_meta <- v2_sub@meta.data

#Create Spatial Experiment
spatial_spe <- SpatialExperiment(
  assay= v2_sub@assays$Spatial@layers$counts.14,
  spatialCoordsNames = c('x','y'),
  colData = v_meta[,c('nUMI', 'seurat_clusters','x','y')]
)

rownames(spatial_spe) <- rownames(v2_sub)

ggplot(spatialCoords(spatial_spe), aes(x, y)) + 
  geom_point(size = 5) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# Create RCTD
rctd_data <- createRctd(spatial_spe, reference_sce, cell_type_col = 'cell_annotation_leiden', 
                        require_int = F)

# Run RCTD
results_spe_v2 <- runRctd(rctd_data, rctd_mode = "multi", max_cores = 8, max_multi_types = 4)
saveRDS(results_spe_v2, '/data/Share/pengmad/visium/rctd_results_v2.rds')

colData(results_spe_v2)

classification_df <- data.frame(
  pixel = colnames(assay(results_spe)),
  spot_class = colData(results_spe)$spot_class,
  first_type = colData(results_spe)$first_type,
  second_type = colData(results_spe)$second_type
)

test <- plotAllWeights(results_spe_v2, 'weights_full', r = 3, lwd = 0, title = "Rctd Results")
test

# Assign colors by cell type order
cell_types <- rownames(results_spe_v2)
palette_named <- setNames(palette3[seq_along(cell_types)], cell_types)
palette_named["reject"] <- "#BEBEBE"  # or "#A0A0A0" for a slightly darker gray
palette_named['ExN.CA1-2.FIBCD1.FNDC1'] <- NULL
palette_named['ExN.CA1.2.FIBCD1.FNDC1'] <- "#CD73FF"
palette_named['ExN.CA3-4.GUCA1C.MYPN'] <- NULL
palette_named['ExN.CA3.4.GUCA1C.MYPN'] <- "#B53DCC"


# check
length(palette_named)
tail(palette_named)


p <- plotAllWeights(
  results_spe_v2,
  r = 100,
  lwd = 0,
  title = "v2 Cell Type Decomposition with Multiome Clusters"
) +
  scale_fill_manual(values = palette_named) +
  scale_color_manual(values = palette_named)

pdf("/data/Share/pengmad/visium/rctd_results_v2_MM131.pdf", width = 16, height = 8)
print(p)
dev.off()

SpatialDimPlot_scCustom(v2, group.by = 'seurat_clusters',label.size = 5, label = T,combine = F)

v2$cluster <- dplyr::case_when(
  v2$seurat_clusters == 0 ~ 'wm',
  v2$seurat_clusters == 1 ~ 'ca1',
  v2$seurat_clusters == 2 ~ 'sm',
  v2$seurat_clusters == 3 ~ 'ca4.ca3',
  v2$seurat_clusters == 4 ~ 'ca3.sm',
  v2$seurat_clusters == 5 ~ 'sgz.ml',
  v2$seurat_clusters == 6 ~ 'sm.int',
  v2$seurat_clusters == 7 ~ 'ca1.out',
  v2$seurat_clusters == 8 ~ 'so',
  v2$seurat_clusters == 9 ~ 'ca1.sr',
  v2$seurat_clusters == 10 ~ 'vasc',
  v2$seurat_clusters == 11 ~ 'gcl',
  v2$seurat_clusters == 12 ~ 'cp',
  v2$seurat_clusters == 13 ~ 'sgz.int',
  v2$seurat_clusters == 14 ~ 'ca1.rostral',
  v2$seurat_clusters == 15 ~ 'ca2',
)

SpatialDimPlot_scCustom(v2, group.by = 'cluster',label.size = 5, label = T,combine = F)

v2_clusters <- read.csv('visium/v2_cluster_names.csv',row.names = 1)
v2_clusters$sample_barrcode <- rownames(v2_clusters)
v2_clusters <- v2_clusters[grepl("^MM131", rownames(v2_clusters)), ]

all(rownames(v2_clusters) == colnames(results_spe_v2))   # should be TRUE
colData(results_spe_v2)$cluster <- v2_clusters$x


# Can filter based on cluster number v2
spe_subset_v2 <- results_spe_v2[, results_spe_v2$cluster == 12]
p <- plotAllWeights(
  spe_subset_v2,
  r = 100,
  lwd = 0,
  title = "v2 Cell Type Decomposition with Multiome Clusters"
) +
  scale_fill_manual(values = palette_named) +
  scale_color_manual(values = palette_named)
p


# Same for v1
all(colnames(v1_sub) == colnames(results_spe))
colData(results_spe)$cluster <- v1_sub$seurat_clusters
spe_subset_v1 <- results_spe[, results_spe$cluster == 0]


p1 <- plotAllWeights(
  results_spe,
  r = .7,
  lwd = 0,
  title = "v1 Cell Type Decomposition with Multiome Clusters"
) +
  scale_fill_manual(values = palette_named) +
  scale_color_manual(values = palette_named)

print(p1)

pdf("/data/Share/pengmad/visium/rctd_results_v1.pdf", width = 16, height = 8)
print(p1)
dev.off()
