suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
  library(Polychrome)
  library(scCustomize)
  library(reticulate)
})

# load seurat object
seurat <- readRDS('/data/Share/Victor_Anosike/V1_Seurat_Spatial/v1_named_clusters_11_17_2025.rds')

# update the seurat object
seurat <- UpdateSeuratObject(seurat)

Idents(seurat) <- "anatomy"

DefaultAssay(seurat) <- "Spatial"

seurat <- NormalizeData(seurat)

list_of_genes <- c('SOX2', 'PAX6', 'NES', 'ASCL1','S100B','GFAP', 'ETNPPL', 'HES5','AQP4', 'VIM',
                   'FABP7', 'NOTCH1','NOTCH2', 'MKI67', 'MCM2', 'PCNA',  'FOXO3', 'NEUROG2', 'OLIG1', 'OLIG2',
                   'NEUROD1', 'TOP2A', 'EOMES', 'NR2E1', 'MYT1L', 'CALB1', 'CALB2',
                   'RELN', 'DCX', 'STMN2', 'TUBB3', 'PROX1', 'FST', 'BHLHE22', 'POSTN',
                   'RBFOX3', 'SYT1', 'MBP', 'SOX11', 'SOX4','GAD1', 'GAD2', 'SHH', 'PTCH1', 'SMO', 'GLI1', 'GLI2',
                   'HHIP', 'HOXA9', 'HOXB4', 'HOXA10', 'HOXC6','HOXD13', 'BMP2', 'BMP4', 'BMP7','BMPR1A','BMPR1B', 'BMPR2',
                   'NOG', 'GREM1', 'WNT3A', 'WNT5A', 'WNT7A','LRP6','FZD7', 'AXIN2','DVL1','DVL2','CTNNB1', 'SEMA3D', 'FOS','HSPB8')

list_of_genes = unique(list_of_genes)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 5000)
length(VariableFeatures(seurat))
#5000
HVGs <- unique(c(VariableFeatures(seurat), list_of_genes))
length(HVGs)
#5073

# subset object to only the listed genes
seurat_subset <- subset(seurat, features = HVGs)

py_install('anndata')

# create anndata object from seurat seurat object
as.anndata(seurat_subset, file_path = "/data/Share/Tiancheng_Shi", file_name = 'v1.h5ad', main_layer = 'data', other_layers = 'counts',
           transer_dimreduc = TRUE)