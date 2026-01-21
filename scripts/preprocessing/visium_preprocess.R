# Preprocess visium samples

library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(Signac)
library(ggplot2)
library(cowplot)
library(patchwork)
library(scDblFinder)
library(SoupX)
library(ggridges)
library(tidyverse)
library(hdf5r)
library(future)
library(harmony)

source('config.R')
source('get_raw_data_paths.R')

############################# HELPER FUNCTIONS #################################

#Input: sample id (string)
read_visium <- function(sample_id) {
  
  path <- get_raw_path_visium(sample_id)
  if (!file.exists(paste0(path, sample_id,'_spaceranger_count_outs'))) {
    print('Decompressing tar file')
    untar(paste0(path, sample_id,'_spaceranger_count_outs.tar'), 
          exdir = path)
  }
  
  dir <- glue('{path}{sample_id}_spaceranger_count_outs')
  
  # Hard coded
  if (substr(sample_id, 1,4)=='MM07'){
    dir <- path
  }
  
  print(dir)
  
  visium_srat <- Load10X_Spatial(data.dir=dir, filename="filtered_feature_bc_matrix.h5", slice=sample_id)
  visium_srat$sample <- sample_id
  qc_umi_size <- visium_srat$nCount_Spatial < 500
  visium_srat$qc_umi_size <- qc_umi_size
  
  saveRDS(visium_srat, paste0(VISIUM_RDS_PATH, sample_id,'_new_visium.rds'))
  return(visium_srat)
}


#Gets clusters
#input: seurat obj, assay (string)
#output: seurat obj
get_clusters <- function(srat, res = 0.8, dims = 30) {
  srat <- NormalizeData(srat, assay = 'Spatial')
  srat <- FindVariableFeatures(srat)
  srat <- ScaleData(srat)
  srat <- RunPCA(srat, assay = 'Spatial')
  srat <- FindNeighbors(srat, dims = 1:dims)
  srat <- FindClusters(srat, resolution = res)
  return(srat)
  
}
################################################################################

# List the sample names 
sample_ids <- c('MM128', 'MM131', 'RM002')

# Loop over each sample id, run processing pipeline
for (id in sample_ids) {
  print(paste("Working on sample", id, "!"))
  if (file.exists(paste0(VISIUM_RDS_PATH, id, "_visium_processed.rds"))) {
    print("Finished!!")
    next
  }
  srat <- read_visium(id)
  
  #Run normalization, dimension reduction, clustering algorithm, UMAP
  srat <- get_clusters(srat)
  srat <- RunUMAP(srat, dims = 1:30)
  
  #Plots for QC
  violin_qc <- VlnPlot(srat, features = c("nCount_Spatial", "nFeature_Spatial"), 
                       stack=TRUE, log = TRUE, pt.size = 0) + NoLegend()
  spatial_qc <- SpatialDimPlot(srat, group.by = "qc_umi_size")
  umap_plot <- DimPlot(srat,reduction = 'umap', label = TRUE, pt.size = 0.3, group.by = 'Spatial_snn_res.0.8')
  spatial_plot <- SpatialDimPlot(srat, group.by = 'Spatial_snn_res.0.8', label = TRUE, 
                                 label.size = 3, pt.size.factor = .8, alpha = 0.8, repel = TRUE)
  
  #save plots in pdf
  pdf(paste0(PLOTS_PATH, id, 'visium_qc_plots.pdf'), height = 7, width = 7)
  print(violin_qc)
  print(spatial_qc)
  print(umap_plot)
  print(spatial_plot)
  dev.off()
  
  # Saved the processed samples
  saveRDS(srat, paste0(VISIUM_RDS_PATH, id, "_visium_processed.rds"))
  print(paste("Done with sample", id, "!!!"))
}





