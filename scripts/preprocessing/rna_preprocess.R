#RNA processing -- methods called from 'rna_methods.R' script

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

# Run config and other methods files
source('../config.R')
source('get_raw_data_paths.R')


############################# HELPER FUNCTIONS #################################

#Gets clusters runs all the normalization, dimension reduction, clustering, etc
#input: seurat obj
#output: seurat obj
get_clusters <- function(srat, res = 0.8, dims = 40) {
  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat)
  srat <- ScaleData(srat)
  srat <- RunPCA(srat)
  srat <- FindNeighbors(srat, dims = 1:dims)
  srat <- FindClusters(srat, resolution = res)
  return(srat)
  
}

#Configure metadata
#Inputs: Seurat Object
#Output: Seurat Object with metadata updated
configure_metadata <- function(srat) {
  
  sample_id <- as.character(srat$orig.ident[1])
  path <- get_raw_path_multiome(sample_id)
  meta <- read.csv(paste0(path, "/per_barcode_metrics.csv"))
  meta <- meta[meta$barcode %in% Cells(srat),]
  
  srat <- RenameCells(srat, add.cell.id = sample_id)
  srat$sample <- sample_id
  meta$sample_barcode <- Cells(srat)
  rownames(meta) <- meta$barcode
  srat$mitoPercent <- PercentageFeatureSet(srat, pattern = '^MT-')
  
  columns_needed <- c("gex_barcode","gex_raw_reads", "atac_barcode", 
                      "atac_raw_reads", "atac_unmapped_reads",
                      "atac_mitochondrial_reads", "atac_fragments", "atac_TSS_fragments",
                      "atac_peak_region_fragments", "atac_peak_region_cutsites")
  
  for (col in columns_needed) {
    srat[[col]] <- meta[[col]]
  }
  
  return(srat)
}


#Run SoupX for ambient RNA removal
#Input: sample id (string)
#Output: Seurat object strained by SoupX with configured metadata
run_soupX <- function(sample_id, auto = TRUE) {
  
  path <- get_raw_path_multiome(sample_id)
  
  if (!file.exists(paste0(path, '/raw_feature_bc_matrix.h5'))) {
    tar_file <- paste0(sample_id, "_cellranger_arc_count_outs.tar")  
    system(paste("tar -xvf", paste0(path,'/',tar_file), '-C', path))
    print(file.exists(paste0(path, '/raw_feature_bc_matrix.h5')))
  }
  
  raw <- Read10X_h5(paste0(path, "/raw_feature_bc_matrix.h5"))
  filtered <- Read10X_h5(paste0(path, "/filtered_feature_bc_matrix.h5"))
  
  sc <- SoupChannel(raw$"Gene Expression", filtered$"Gene Expression")
  srat <- CreateSeuratObject(sc$toc)
  srat <- get_clusters(srat)
  sc$metaData$clusters <- Idents(srat)
  sc <- setClusters(sc, sc$metaData$clusters)
  
  if (auto == TRUE) {
    #save soupX plots
    pdf(paste0(PLOTS_PATH, sample_id, "_soupX.pdf"), width = 8, height = 5.5)
    sc <- autoEstCont(sc)
    out <- adjustCounts(sc)
    dev.off()
  }
  
  else {
    sc <-  setContaminationFraction(sc, 0.2)
    out <- adjustCounts(sc)
    
  }
  soup_srat <- CreateSeuratObject(out, project = sample_id)
  soup_srat <- configure_metadata(soup_srat)    
  
  return(soup_srat)
}



#Detects doublets with scDblFinder, adds column in metadata for doublet/singlet
#Input: Seurat object
#Output: Seurat object
find_doublets <- function(srat) {
  
  sample_id <- as.character(srat$orig.ident[1])
  sce <- as.SingleCellExperiment(srat)
  sce <- scDblFinder(sce, clusters = 'seurat_clusters')
  srat$scDblFinder.class <- sce$scDblFinder.class
  srat$scDblFinder.score <- sce$scDblFinder.score
  srat$scDblFinder.weighted <- sce$scDblFinder.weighted
  srat$scDblFinder.cxds_score <- sce$scDblFinder.cxds_score
  
  return(srat)
}
################################################################################

# List of marker genes for plotting
genes <- c('PDGFRA','OLIG1','MBP','MOG','MAL','PTPRC','P2RY12','APBB1IP','CD163','CD247',
           'PECAM1','FLT1','NES','CEMIP',
           'TTR','HTR2C','CFAP54','AQP4','SLC1A2','ALDH1L1',
           'GLUL','GFAP','PAX6','SOX2','MYT1L','RBFOX3',
           'SOX1','DCX','GAD1','GAD2','SLC32A1',
           'RELN','CALB2','NR2E1','PROX1',
           'CCK','SLC17A8','VIP','PVALB','LHX6','SST','PENK',
           'BHLHE22','POSTN','CALB1','SEMA5A','SLC17A7',
           'CAMK2A','SV2B','SATB2')

# List the sample IDs that are being processed
sample_ids <- c('CM021', 'CM024', 'CM025', 'CM073', 'CM115', 
                'MM001', 'MM008', 'MM023','MM145','MM151')


# Loop through all the sample ids
for (id in sample_ids) {
  
  print(paste("Working on sample", id, "!!!"))
  
  if (file.exists(paste0(RNA_PROCESSED_RDS_PATH, id, "_processed.rds"))) {
    print("Finished!!")
    next
  }
  
  # Remove ambient RNA
  soupx_obj <- run_soupX(id)
  
  #QC metrics filter
  soupx_filt <- subset(soupx_obj, subset = nCount_RNA > 500 &
                         nCount_RNA < 150000 &
                         nFeature_RNA > 300 &
                         mitoPercent < 10)
  
  
  # Normalize, Dimension Reduction, Cluster, UMAP on individual sample
  soupx_filt <- get_clusters(soupx_filt)
  soupx_filt <- RunUMAP(soupx_filt, dims = 1:30)
  
  
  # Find doublets -- to be removed later
  soupx_filt <- find_doublets(soupx_filt)
  
  
  # Create QC plots -- violin, umap, gene expression dot plots
  vln_plot_unfiltered <- VlnPlot(soupx_obj, features = c("nCount_RNA", "nFeature_RNA", "mitoPercent"),
                                 ncol = 1, log = TRUE, pt.size = 0) + ggtitle(paste(id, "Violin Plot Pre QC Filtering"))
  
  vln_plot_filtered <- VlnPlot(soupx_filt, features = c("nCount_RNA", "nFeature_RNA", "mitoPercent"), 
                               ncol = 1, log = TRUE, pt.size = 0, group.by = 'sample') + ggtitle(paste(id, "Violin Plot Post QC Filtering"))
  
  umap_plot_clusters <- DimPlot(soupx_filt, group.by = 'seurat_clusters', label = T, repel = T) +
    ggtitle(paste(id, "Seurat Clusters"))
  
  dot_plot_clusters <- DotPlot(object = soupx_filt, group.by = 'seurat_clusters',
                               features = rev(genes), cols = c("yellow","blue")) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90))
  
  umap_plot_doublets <- DimPlot(soupx_filt, group.by = 'scDblFinder.class', label = F, repel = T) +
    ggtitle(paste(id, "Doublets"))
  
  
  # Print all the plots in 1 pdf
  pdf(paste0(PLOTS_PATH, id, "_qc_plots.pdf"), width=7, height=7)
  print(vln_plot_unfiltered)
  print(vln_plot_filtered)
  print(umap_plot_clusters)
  print(dot_plot_clusters)
  print(umap_plot_doublets)
  dev.off()
  
  saveRDS(soupx_filt, paste0(RNA_PROCESSED_RDS_PATH, id, "_processed.rds"))
  print(paste("Done with sample", id, "!!!"))
  
}