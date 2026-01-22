
# List sample ids
sample_ids <- c()

# Read in individual preprocessed objects
srat_objects <- lapply(sample_ids, function(i) {
  path <- paste0(VISIUM_RDS_PATH, i, "_visium_processed.rds")
  srat <- readRDS(path)
  return(srat)
})

# Merge objects
spatial <-merge(
  x = srat_objects[[1]],
  y = srat_objects[2:length(srat_objects)],
  add.cell.ids = sample_ids
)

# Remove spots with less than 500 UMI
spatial <- subset(spatial, subset = nCount_Spatial > 500)

# Join layers, Normalize, PCA, etc..
spatial <- JoinLayers(spatial)
spatial <- NormalizeData(spatial, verbose = FALSE, assay = "Spatial")
spatial <- FindVariableFeatures(spatial)
spatial <- ScaleData(spatial)
spatial <- RunPCA(spatial, assay = "Spatial", verbose = FALSE)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:30)

