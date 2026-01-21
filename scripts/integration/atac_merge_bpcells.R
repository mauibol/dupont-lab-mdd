# Maddy Peng
# Merge ATAC samples together
# Use BPCells to stream matrix from disk 

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
library(ggplot2)
library(gridExtra)

# List sample ids
sample_ids <- c()

#Read in Seurat Objects
atac_srats <- lapply(sample_ids, function(id) {
  srat <- readRDS(paste0(ATAC_SAVED_RDS_PATH, id, "_atac_combined_peaks.rds"))
  return(srat)
})

#Write each counts matrix to disk with bpcells package and create seurat object
atac_bpcells <- list()
for (srat in atac_srats) {
  id <- unique(srat$sample)
  mtx <- srat[['ATAC']]$counts
  write_matrix_dir(mat = mtx, dir = paste(MERGED_RDS_PATH, 'on_disk_mtx/', id, '/bpcells_atac'))
  data <- open_matrix_dir(dir= paste(MERGED_RDS_PATH, 'on_disk_mtx/', id, '/bpcells_atac'))
  obj <- CreateSeuratObject(counts=data, meta=srat@meta.data, assay = 'ATAC')
  print(paste('Appending sample', id))
  atac_bpcells <- append(atac_bpcells, obj)
}

# Clear memory
rm(atac_srats)
gc()

# Merge the bpcells seurat objects and join layers if v5
atac_merge <- merge(
  x = atac_bpcells[[1]],
  y = atac_bpcells[2:length(atac_bpcells)],
  add.cell.ids = sample_ids
)

#Join the layers into 1 matrix
atac_merge <- JoinLayers(atac_merge)

# Write merged matrix to disk with bpcells and read it in to normalize and reduce
write_matrix_dir(mat=atac_merge[['ATAC']]$counts, dir=paste0(MERGED_RDS_PATH, 'merged/bpcells_atac'))

# Save merged object
saveRDS(atac_merge, paste0(MERGED_RDS_PATH, 'merged_atac_bpcells.rds'))

