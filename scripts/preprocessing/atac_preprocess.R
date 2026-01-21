# ATAC preprocessing

library(GenomicRanges)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(tidyverse)
library(biovizBase)
library(harmony)

# Genomic annotations 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
#genome(annotations) <- "hg38"


############################# HELPER FUNCTIONS #################################

# Creates a peak set from MACS call peaks function on all fragments
create_peak_set <- function(list, method = c('macs', 'cellranger')) {
  
  if (method == 'macs') {
    frag_paths <- lapply(list, function(x) {
      path <- get_raw_path_multiome(x)
      return(paste0(path, "/atac_fragments.tsv.gz"))
    })
    
    peaks <- CallPeaks(frag_paths, 
                       macs2.path = MACS_2_PATH, additional.args= '--min-len 100 --max-gap 50')
    
  } else if (method == 'cellranger') {
    
    peaks <- lapply(list, function(x) {
      path <- get_raw_path_multiome(x)
      peak <- read.table(
        file = paste0(path, "/atac_peaks.bed"),
        skip = 52,
        col.names = c("chr", "start", "end")
      )
      gr <- makeGRangesFromDataFrame(peak)
      return(gr)
    })
    
    peaks <- GRangesList(peaks)
    peaks <- unlist(peaks)
    peaks <- GenomicRanges::reduce(peaks)
  }
  
  #Filter based on peak widths for COMBINED PEAKS
  filtered_peaks <- peaks[width(peaks) < 10000 & width(peaks) > 20]
  filtered_peaks <- keepStandardChromosomes(filtered_peaks, pruning.mode = 'coarse')
  return(filtered_peaks)
}


# Creates seurat object from id and unified peak set
# Make sure to change paths accordingly
create_atac_srat <- function(id, peak_set) {
  
  if (file.exists(paste0(ATAC_PROCESSED_RDS_PATH, id, "_atac_combined_peaks.rds"))) {
    srat <- LoadSeuratRds(paste0(ATAC_PROCESSED_RDS_PATH, id, "_atac_combined_peaks.rds"))
    
  } else {
    
    print(paste0("Reading in metadata for ", id))
    
    # Read in metadata
    metadata_path <- paste0(get_raw_path_multiome(id), "/per_barcode_metrics.csv")
    meta <- read.table(
      file = metadata_path,
      stringsAsFactors = FALSE,
      sep = ",",
      header = TRUE,
      row.names = 1)
    
    # Create Fragment object
    print(paste0("Creating Fragments object for ", id))
    fragments <- CreateFragmentObject(
      path = paste0(paste0(get_raw_path_multiome(id), "/atac_fragments.tsv.gz")),
      cells = meta$gex_barcode
    )
    
    # Create peaks matrix
    print(paste0("Creating Feature Matrix for ", id))
    peak_counts <- FeatureMatrix(
      fragments = fragments,
      features = peak_set,
      cells = meta$gex_barcode
    )
    
    #Create chromatin assay
    print(paste0("Creating Seurat Object for ", id))
    sample_assay <- CreateChromatinAssay(
      peak_counts, 
      fragments = fragments,
      annotation = annotations)
    
    # Create Seurat Object from chromatin assay
    srat <- CreateSeuratObject(
      sample_assay, 
      assay = "ATAC", 
      project = id, 
      meta.data=meta)
    
    
    # QC metrics filter
    srat <- NucleosomeSignal(srat) #Low NS = high chromatin accessibility
    srat <- TSSEnrichment(srat, fast = FALSE) #High TSSE = high chromatin accessibility
    srat$pct_reads_in_peaks <- srat$atac_peak_region_fragments / srat$atac_fragments * 100
    srat$blacklist_fraction <- FractionCountsInRegion(srat, assay = 'ATAC', regions = blacklist_hg38_unified)
    srat$sample = id
    
    srat_filt <- subset(srat, subset = nCount_ATAC > 2e2 &
                     nCount_ATAC < 1e5 &
                     nucleosome_signal < 3 &
                     TSS.enrichment > 1)
    
    # save filtered object
    saveRDS(srat_filt, paste0(ATAC_SAVED_RDS_PATH, id, "_atac_combined_peaks.rds"))
    
  }
  print(paste0("Done with sample", id))
  return(srat_filt)
  
}
################################################################################

# List the IDs for processing
sample_ids <- c('CM021')
print(ids)

rna_meta$gex_barcode <- substr(rownames(rna_meta), 7,24)

# Filtered unified peak set from cellranger method
combined_peaks <- create_peak_set(ids, method = 'cellranger')

atac_srats <- lapply(ids, function(id) {
  srats <- create_atac_srat(id, combined_peaks)
  return(srats)
})


# If working with the macs3 peak set, can merge the ATAC seurat_object here
atac_merge <- merge(
  x = atac_srats[[1]],
  y = atac_srats[2:length(atac_srats)],
  add.cell.ids = ids
)

# Save the merged Atac
saveRDS(atac_merge, paste0(MERGED_RDS_PATH, "merged_atac_macs3_iter9.rds"))

# Free up space
rm(atac_srats)
gc()
options(future.globals.maxSize = 1000 * 1024^2)

# Feature selection, dim reduction, clustering etc
atac_merge <- FindTopFeatures(atac_merge, min.cutoff = 20)
atac_merge <- RunTFIDF(atac_merge)
atac_merge <- RunSVD(atac_merge)
atac_merge <- RunUMAP(atac_merge, dims = 2:50, reduction = 'lsi', verbose = F, reduction.name = "umap.atac")
atac_merge <- FindNeighbors(atac_merge, dims=2:50, reduction='lsi')
atac_merge <- FindClusters(atac_merge, algorithm=3)

# Plot pre-integration
pdf(paste0(PLOTS_PATH, "atac_merged_pre_integration.pdf"))
plt <- DimPlot(atac_merge, group.by = 'sample', pt.size = 0.5, reduction = 'umap.atac', raster = F)
print(plt)
dev.off()

# Integrate with harmony
atac_merge <- RunHarmony(atac_merge, group.by.vars = c('sample', 'ID4', 'batch'), reduction = 'lsi', dims.use = 2:50,
                         plot_convergence = T, reduction.save = "harmony.atac", assay.use = "ATAC", verbose = TRUE, project.dim = FALSE)


# Find clusters with harmony data
atac_merge <- FindNeighbors(atac_merge, reduction = "harmony.atac", dims = 1:49)
atac_merge <- FindClusters(atac_merge, resolution = 0.8, cluster.name = 'atac_harmony_clusters_0.8')
atac_merge <- RunUMAP(atac_merge, reduction = "harmony.atac", dims = 1:49, 
                      reduction.name = 'umap.atac.harmony', 
                      reduction.key = 'atacHarmonyUMAP_')


# Plot clusters post integration
pdf(paste0(PLOTS_PATH, "atac_harmony.pdf"), width = 8, height = 6)
p2<- DimPlot(atac_merge, reduction = "umap.atac.harmony", group.by =c('sample', 'atac_harmony_clusters_0.8'), 
             label = TRUE,
             label.size = 3, repel = TRUE, raster = F, ncol=2) +
  ggtitle("Post Integration")
print(p2)
dev.off()