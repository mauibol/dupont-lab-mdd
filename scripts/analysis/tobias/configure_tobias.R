# Maddy Peng
# Calls Peaks on Individual Cell Clusters
# Isolates cell barcodes per cell cluster for each sample

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(here)
  library(tidyr)
  library(dplyr)
  library(scCustomize)
})

setwd('/groups/mb928_gp/mp4486/')
source("Multiome_process/scripts/config.r")

# Read in ATAC and Metadata
atac_merge <- LoadSeuratRds(paste0(MERGED_RDS_PATH, "merged_atac_macs3.rds"))
meta <- read.csv("Multiome_process/meta/current_working_multiome_meta (1).csv", row.names=1)
atac_merge@meta.data <- meta

# Run call peaks on a subset of cell types
cellTypes <- c('ExN.1', 'InN.5', 'GC.1')
grs <- CallPeaks(atac_merge, group.by = 'cell_annotation_leiden', idents=cellTypes, macs2.path = MACS_2_PATH, combine.peaks=FALSE)

# Configure granges peaks and save peaks for each cluster
for (i in 1:length(grs)) {
  gr <- grs[[i]]
  gr_df <- as.data.frame(gr)[,c("seqnames", "start", "end")]
  colnames(gr_df) <- c("Chromosome", "Start", "End")
  
  path = paste0("CURRENT/peak_files/", gr$ident[1], ".bed")
  write.table(gr_df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Extract Barcodes -- 1 txt file for each cell type of interest per sample
split <- SplitObject(atac_merge, split.by='sample')
for (c in cellTypes) {
  for (s in split) {
    cells <- s$cell_annotation_leiden[s$cell_annotation_leiden == c]
    if (length(cells) == 0) {
      print(paste('There are 0 cells in', c, "cluster of sample", s$sample[1], "skipping this sample..."))
      next
    }
    print(paste('There are', length(cells), 'cells in', c, "cluster of sample", s$sample[1]))
    print('Reading barcodes')
    srat <- subset(s, subset = cell_annotation_leiden == c)
    print(dim(srat))
    barcodes <- as.list(srat$gex_barcode)
    print(length(barcodes))
    barcodes <- lapply(barcodes, function(x) paste0("CB:Z:", x))
    file_path <- paste0("CURRENT/barcodes/", srat$sample[1], "_", gsub("[ /]", "_", c), ".txt")
    if (!file.exists(file_path)) {
      cat(unlist(barcodes), file = file_path, sep = "\n")
    }
    
  }
}
