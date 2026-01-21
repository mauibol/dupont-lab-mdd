# Merge RNA samples together
# Use BPCells to stream matrix from disk 

# List sample ids
sample_ids <- c(
  "CM021", "CM023", "CM024", "CM025", "CM026", "CM027", "CM028", "CM038",
  "CM039", "CM040", "CM041", "CM042", "CM043", "CM044", "CM053", "CM055",
  "CM056", "CM057", "CM058", "CM059", "CM060", "CM069", "CM070", "CM071",
  "CM072", "CM073", "CM074", "CM075", "CM085", "CM086", "CM087", "CM088",
  "CM089", "CM090", "CM091", "CM092", "CM109", "CM110", "CM111", "CM112",
  "CM113", "CM114", "CM115", "CM116", "MM001", "MM002", "MM006", "MM007",
  "MM008", "MM016", "MM023", "MM024", "MM080", "MM081", "MM082", "MM083",
  "MM084", "MM085", "MM086", "MM087", "MM096", "MM097", "MM100", "MM101",
  "MM102", "MM103", "MM145", "MM146", "MM147", "MM148", "MM149", "MM150",
  "MM151", "MM152"
)

#Read in Seurat Objects
rna_srats <- lapply(sample_ids, function(id) {
  srat <- readRDS(paste0(RNA_SAVED_RDS_PATH, id, "_rna_processed.rds"))
  return(srat)
})

#Write each counts matrix to disk with bpcells package and create seurat object
rna_bpcells <- list()
for (srat in rna_srats) {
  id <- unique(srat$sample)
  mtx <- srat[['RNA']]$counts
  write_matrix_dir(mat = mtx, dir = paste(MERGED_RDS_PATH, 'on_disk_mtx/', id, '/bpcells_rna'))
  data <- open_matrix_dir(dir= paste(MERGED_RDS_PATH, 'on_disk_mtx/', id, '/bpcells_rna'))
  obj <- CreateSeuratObject(counts=data, meta=srat@meta.data, assay = 'RNA')
  print(paste('Appending sample', id))
  rna_bpcells <- append(rna_bpcells, obj)
}

# Clear memory
rm(rna_srats)
gc()

# Merge the bpcells seurat objects and join layers if v5
rna <- merge(
  x = rna_bpcells[[1]],
  y = rna_bpcells[2:length(rna_bpcells)],
  add.cell.ids = sample_ids
)

#Join the layers into 1 matrix
rna <- JoinLayers(rna)

# Write merged matrix to disk with bpcells and read it in to normalize and reduce
write_matrix_dir(mat=rna[['RNA']]$counts, dir=paste0(MERGED_RDS_PATH, 'merged/bpcells_rna'))

# Save merged object
saveRDS(rna, paste0(MERGED_RDS_PATH, 'rna_bpcells.rds'))

# Preliminary processing for cleaning, visualizing and removing doublets
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna)
rna <- ScaleData(rna)
rna <- RunPCA(rna, npcs = 40)

rna <- RunHarmony(rna, group.by.vars = c('sample', 'batch', 'ID4'),
                   plot_convergence = T, reduction.save = "harmony.rna")

rna <- FindNeighbors(rna, reduction = 'harmony.rna', dims = 1:30)
rna <- FindClusters(rna, cluster.name = 'prelim_clusters_rna')
rna <- RunUMAP(rna, reduction = 'harmony.rna', dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# Plots
DimPlot_scCustom(rna9, reduction = 'harmony.rna.umap', group.by = 'prelim_clusters_rna', label = T)
DimPlot_scCustom(rna9, reduction = 'harmony.rna.umap', group.by = 'scDblFinder.class', label = T)

# Filter out the doublets both detected by doublet finder and based on clusters
rna$to_keep <- dplyr::case_when(
  rna$prelim_clusters_rna==35 ~ FALSE,
  rna$scDblFinder.class == 'doublet' ~ FALSE,
  .default=TRUE
)

rna <- subset(rna, subset = to_keep == TRUE)
saveRDS(rna, paste0(MERGED_RDS_PATH, 'rna_bpcells_filtered.rds'))





