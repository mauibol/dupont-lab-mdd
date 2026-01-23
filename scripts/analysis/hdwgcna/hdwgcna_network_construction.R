#https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html

install.packages("BiocManager")
BiocManager::install(c("WGCNA", "igraph", "devtools", "GeneOverlap", "ggrepel", "UCell"))
devtools::install_github("NightingaleHealth/ggforestplot")
BiocManager::install("GO.db")
BiocManager::install("GenomeInfoDbData")
install.packages("tidyverse")
BiocManager::install("limma")
BiocManager::install("variancePartition")

# install Seurat v5 
install.packages("ggplot2")
remotes::install_github("bnprks/BPCells/r")
devtools::install_github('smorabit/hdWGCNA', ref='dev', dependencies = T)




# single-cell analysis package
library(Seurat)
library(BPCells)
# plotting and data science packages
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(limma)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
library(variancePartition)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the Zhou et al snRNA-seq dataset
#res.file <- readRDS(rds.file)
setwd("~/hdWGCNA/4.22.25")

#whole_seurat_obj <- readRDS("/home/marimad@nyspi.local/preprocessing/rna_counts_and_wnn_umap.rds")

bpcells <- readRDS("/data/Share/pengmad/merged_CURRENT/rna_bpcells_qc_iteration9_march_19_2025.rds")
leiden_clusters <- readRDS("/data/Share/pengmad/merged_CURRENT/multimodal_bpcells_iteration9_march19_2025.rds")
colnames(leiden_clusters@meta.data)


expr_matrix <- bpcells@assays[["RNA"]]@layers[["counts"]]
expr_matrix <- expr_matrix@matrix
#dim(expr_matrix)
# 36601 495037


# Calculate fraction of cells in which each gene is expressed (>0)
expr_fraction <- Matrix::rowMeans(expr_matrix > 0)

# Filter genes: keep those expressed in at least 5% of cells
filtered_expr_matrix <- expr_matrix[expr_fraction >= 0.05, ]

#dim(filtered_expr_matrix)
# 12784 495037
filtered_genes <- rownames(filtered_expr_matrix)


rownames(expr_matrix) <- rownames(bpcells)
colnames(expr_matrix) <- colnames(bpcells)

metadata <-leiden_clusters@meta.data

# Create Seurat object from counts matrix and metadata
seurat_obj <- CreateSeuratObject(
  counts = expr_matrix,  # Raw counts data
  meta.data = metadata      # Metadata information
)

seurat_obj[["harmony.rna"]] <- leiden_clusters@reductions[["harmony.rna"]]
seurat_obj[["harmony.rna.umap"]] <- leiden_clusters@reductions[["harmony.rna.umap"]]
seurat_obj[["pca"]] <- leiden_clusters@reductions[["pca"]]

rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]@matrix)  
colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]@matrix)  

#Check if column names match between Seurat and the count matrix
all(colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]@matrix) == colnames(seurat_obj))
all(rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]@matrix) == rownames(seurat_obj))



p <- DimPlot(seurat_obj, group.by='cell_annotation_leiden', label=TRUE, reduction = "harmony.rna.umap") +
  umap_theme() + ggtitle('Cell_annotation_leiden') + NoLegend()


p
###############################################################################################
#Set WGCNA 

wgcna_name <- "tutorial_GN"
seurat_obj <- SetActiveWGCNA(seurat_obj, wgcna_name)

seurat_obj@misc[["tutorial_GN"]][["wgcna_genes"]] <- filtered_genes

seurat_obj@assays[["RNA"]]@layers[["counts"]] <- seurat_full@assays[["RNA"]]@layers[["counts"]]


#saveworkspace_pre_metacell_contrsuction

################################################################################################
#Construct Metacells with BPcells 

#!!DONT FORGET TO FIRST RUN FUNCTIONS FOR BPCELLS MTECACELL CONSTRUCTION!!

ConstructMetacells <- function(
    seurat_obj, name='agg', ident.group='broad_type_leiden', mat=seurat_obj@assays[["RNA"]]@layers[["counts"]]s, BPCells = T,
    reduction='harmony.rna', dims=1:100, 
    cells.use = NULL,
    assay='RNA', layer='counts', mode = 'sum',
    meta=NULL, return_metacell=F,
    k=50, max_shared = 15,
    target_metacells=1000,
    max_iter=5000,
    verbose=T,
    wgcna_name=NULL
){
  
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  
  # check reduction
  if(!(reduction %in% names(seurat_obj@reductions))){
    stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ", paste(names(seurat_obj@reductions), collapse=', ')))
  }
  
  # check assay
  if(!(assay %in% names(seurat_obj@assays))){
    stop(paste0("Invalid assay (", assay, "). Assays in Seurat object: ", paste(names(seurat_obj@assays), collapse=', ')))
  }
  
  # check layer
  if(!(layer %in% c('counts', 'data', 'scale.data'))){
    stop(paste0("Invalid layer (", layer, "). Valid options for layer: counts, data, scale.data "))
  }
  
  # subset seurat object by selected cells:
  if(!is.null(cells.use)){
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[,cells.use]
  }
  
  message(paste0("Aggregating cells from grouping ", name))
  
  reduced_coordinates <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1,
                   replace = F)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- k * 2
  get_shared <- function(other, this_choice) {
    k2 - length(union(cell_sample[other, ], this_choice))
  }
  while (length(good_choices) > 0 & length(chosen) < target_metacells & it < max_iter) {
    it <- it + 1
    choice <- sample(seq_len(length(good_choices)), size = 1,
                     replace = F)
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen, ]
    others <- seq_len(nrow(cell_sample) - 1)
    this_choice <- cell_sample[nrow(cell_sample), ]
    shared <- sapply(others, get_shared, this_choice = this_choice)
    
    if(max(shared) <= max_shared){
      chosen <- new_chosen
    }
  }
  
  shared_old <- shared
  cell_sample <- nn_map[chosen, , drop = F]
  
  # get a list of the cell barcodes that have been merged:
  cells_merged <- apply(cell_sample, 1, function(x){
    paste0(colnames(seurat_obj)[x], collapse=',')
  })
  
  combs <- tryCatch(
    {combn(nrow(cell_sample), 2)},
    error = function(cond){return(NA)}
  )
  if(any(is.na(combs))){
    warning('Metacell failed')
    return(NULL)
  }
  
  shared <- apply(combs, 2, function(x) {
    k2 - length(unique(as.vector(cell_sample[x, ])))
  })
  
  if(verbose){
    message(paste0("Overlap QC metrics:\nCells per bin: ",
                   k, "\nMaximum shared cells bin-bin: ", max(shared),
                   "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
                   median(shared)))
    if (mean(shared)/k > 0.1)
      warning("On average, more than 10% of cells are shared between paired bins.")
  }
  
  # get original expression matrix
  exprs_old <- mat[,colnames(seurat_obj)]
  
  # groups of cells to combine
  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
                   cell_sample[x, , drop = F])
  # mask <- mask[,which(shared_old <= max_shared)]
  # cell_sample <- cell_sample[which(shared_old <= max_shared),]
  mask <- Matrix::Matrix(mask)
  mask <- as(mask, "dMatrix")
  if(BPCells == T) {
    mask <- write_matrix_dir(mask, tempfile("Mask"))
  }
  
  # average or sum expression?
  new_exprs <- (exprs_old %*% mask)
  if(mode == 'average'){
    new_exprs <- new_exprs / k
  }
  colnames(new_exprs) <- paste0(name, '_', 1:ncol(new_exprs))
  rownames(cell_sample) <- paste0(name, '_', 1:ncol(new_exprs))
  colnames(cell_sample) <- paste0('knn_', 1:ncol(cell_sample))
  
  if(BPCells == T){
    new_exprs <- write_matrix_dir(new_exprs, tempfile("new_exprs"))
  }
  
  message(paste0("Creating Seurat object with ", ncol(new_exprs), " metacells"))
  
  # make seurat obj:
  metacell_obj <- CreateSeuratObject(
    counts = new_exprs,
    assay = assay,
    min.cells = 0,
    min.features = 0
  )
  if(layer == 'scale.data'){
    metacell_obj <- SeuratObject::SetAssayData(
      metacell_obj,
      layer=layer,
      assay=assay,
      new.data=as.matrix(new_exprs)
    )
  }
  
  # add the cells merged info to the metacell obj
  metacell_obj$cells_merged <- as.character(cells_merged)
  
  # calculate stats:
  # shared <- shared[shared <= max_shared]
  max_shared <- max(shared)
  median_shared <- median(shared)
  mean_shared <- mean(shared)
  
  # calculate matrix density:
  new_exprs_mem <- write_matrix_memory(new_exprs) %>% as.matrix()
  new_exprs_mem[new_exprs_mem > 0] <- 1
  density <- sum(Matrix::colSums(new_exprs_mem) / (nrow(new_exprs_mem)*ncol(new_exprs_mem)))
  run_stats <- data.frame(
    name = name,
    max_shared = max_shared,
    mean_shared = mean_shared,
    median_shared = median_shared,
    density = density,
    n = ncol(new_exprs)
  )
  rm(new_exprs_mem)
  
  # add to metacell seurat obj
  metacell_obj@misc$run_stats <- run_stats
  
  # add meta-data:
  if(!is.null(meta)){
    meta_names <- names(meta)
    for(x in meta_names){
      metacell_obj@meta.data[[x]] <- meta[[x]]
    }
  } else(
    warning('meta not found')
  )
  
  # add seurat metacell object to the main seurat object:
  if(return_metacell){
    out <- metacell_obj; gc()
  } else{
    
    # revert to full seurat object if we subsetted earlier
    if(!is.null(cells.use)){
      seurat_obj <- seurat_full
    }
    
    # add seurat metacell object to the main seurat object:
    seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
    
    # add other info
    seurat_obj <- SetWGCNAParams(
      seurat_obj, params = list(
        'metacell_k' = k,
        'metacell_reduction' = reduction,
        'metacell_layer' = layer,
        'metacell_assay' = assay
      ),
      wgcna_name
    )
    out <- seurat_obj; gc()
  }
  out
}

MetacellsByGroups <- function(
    seurat_obj, group.by=c('broad_type_leiden', 'sample'), ident.group='broad_type_leiden',
    matseurat_obj@assays[["RNA"]]@layers[["counts"]], BPCells = T, join = T, matrix_save_path = "BPCells/Metacells", conv_2_sparse = T, 
    reduction='harmony.rna', dims=1:100,
    cells.use = NULL, min_cells=100,
    assay="RNA", layer ='counts', mode = 'sum',
    k=50, k_max = 75, k_min = 30, max_shared_ratio = 0.25,
    target_metacells=1000, max_iter=5000, verbose=T, wgcna_name=NULL
){
  
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  
  # check group.by for invalid characters:
  if(any(grepl('#', group.by))){
    stop('Invalid character # found in group.by, please re-name the group.')
  }
  
  # check ident.group
  if(!(ident.group %in% group.by)){
    stop('ident.group must be in group.by')
  }
  
  # check mode
  if(!(mode %in% c('sum', 'average'))){
    stop('Invalid choice for mode. Mode can be either sum or average.')
  }
  
  # check reduction
  if(!(reduction %in% names(seurat_obj@reductions))){
    stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ", paste(names(seurat_obj@reductions), collapse=', ')))
  }
  
  # check assay:
  if(is.null(assay)){
    assay <- DefaultAssay(seurat_obj)
  } else if(!(assay %in% names(seurat_obj@assays))){
    stop(paste0('Assay ', assay, ' not found in seurat_obj. Select a valid assay: ', paste0(names(seurat_obj@assays), collapse = ', ')))
  }
  
  # check layer:
  if(!(layer %in% c('counts', 'data', 'scale.data'))){
    stop('Invalid input for layer. Valid choices are counts, data, scale.data.')
  } else{
    
    # check the shape of the layer
    layer_dim <- dim(GetAssayData(seurat_obj, assay=assay, layer=layer))
    if(any(layer_dim) == 0){
      stop(paste(c("Selected layer ", layer, " not found in this assay.")))
    }
  }
  
  # check that k > min_cells 
  if(min_cells < k ){
    warning("min_cells is smaller than k, this may result in downstream errors if very small groups are allowed.")
  }
  
  # subset seurat object by selected cells:
  if(!is.null(cells.use)){
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[,cells.use]
  }
  
  # setup grouping variables
  if(length(group.by) > 1){
    seurat_meta <- seurat_obj@meta.data[,group.by]
    for(col in colnames(seurat_meta)){
      seurat_meta[[col]] <- as.character(seurat_meta[[col]])
    }
    seurat_obj$metacell_grouping <- apply(seurat_meta, 1, paste, collapse='#')
  } else {
    seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
  }
  groupings <- unique(seurat_obj$metacell_grouping)
  groupings <- groupings[order(groupings)]
  
  # remove groups that are too small:
  group_counts <- table(seurat_obj$metacell_grouping) < min_cells
  if(any(group_counts)){
    warning(paste0("Removing the following groups that did not meet min_cells: ", paste(names(group_counts)[group_counts], collapse=', ')))
  }
  groupings <- groupings[table(seurat_obj$metacell_grouping) >= min_cells]
  
  if(length(groupings) == 0 ){
    stop("No groups met the min_cells requirement.")
  }
  
  # unique meta-data for each group
  meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, '#')))
  colnames(meta_df) <- group.by
  
  # list of meta-data to pass to each metacell seurat object
  meta_list <- lapply(1:nrow(meta_df), function(i){
    x <- list(as.character(meta_df[i,]))[[1]]
    names(x) <- colnames(meta_df)
    x
  })
  
  # split seurat obj by groupings
  seurat_list <- lapply(groupings, function(x){seurat_obj[,seurat_obj$metacell_grouping == x]})
  names(seurat_list) <- groupings
  
  # Adjust k and max_shared parameters for each sample
  k_list <- numeric()
  max_shared_list <- numeric()
  cells_per_sample <- numeric()
  mean_cells <- numeric()
  for(i in 1:length(seurat_list)) {
    cells_per_sample[i] <- ncol(seurat_list[[i]])
  }
  mean_cells <- mean(cells_per_sample)
  
  for(i in 1:length(cells_per_sample)) {
    k_list[i] <- ifelse(round(k*cells_per_sample[i]/mean_cells) > k_min,
                        round(k*cells_per_sample[i]/mean_cells), 
                        k_min)
    k_list[i] <- ifelse(k_list[i] > k_max, 
                        k_max, 
                        k_list[i])
    max_shared_list[i] <- round(k_list[i]*max_shared_ratio)
  }
  
  # Adjust target_metacell and max_iter parameters for each sample
  target_metacells_list <- numeric()
  max_iter_list <- numeric()
  for(i in 1:length(cells_per_sample)) {
    if(target_metacells < 1000) {
      target_metacells_list[i] <- max(target_metacells, round((target_metacells/100)*(cells_per_sample[i]/k_list[i])))
      max_iter_list[i] <- target_metacells_list[i]*5
    } else {
      target_metacells_list[i] <- target_metacells
      max_iter_list[i] <- max_iter
    }
  }
  
  # construct metacells
  metacell_list <- mapply(
    ConstructMetacells,
    seurat_obj = seurat_list,
    name = groupings,
    meta = meta_list,
    k=k_list,
    max_shared=max_shared_list,
    max_iter=max_iter_list, 
    target_metacells=target_metacells_list,
    MoreArgs = list(
      reduction=reduction, 
      dims=dims,
      assay=assay, 
      layer=layer, 
      mat=mat,
      BPCells=BPCells,
      return_metacell=F, 
      mode=mode,
      verbose=verbose, 
      wgcna_name=wgcna_name
    )
  )
  names(metacell_list) <- groupings
  
  rm(seurat_list); gc()
  
  # remove NULL
  remove <- which(sapply(metacell_list, is.null))
  if(length(remove) >= 1){
    metacell_list <- metacell_list[-remove]
  }
  
  # get the run stats:
  run_stats <- as.data.frame(do.call(rbind, lapply(metacell_list, function(x){x@misc$run_stats})))
  rownames(run_stats) <- 1:nrow(run_stats)
  for(i in 1:length(group.by)){
    run_stats[[group.by[i]]] <- do.call(rbind, strsplit(as.character(run_stats$name), '#'))[,i]
  }
  
  message(paste("Merging metacells"))
  
  # combine metacell objects
  if(length(metacell_list) > 1){
    metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])
    gc()
  } else{
    metacell_obj <- metacell_list[[1]]
    gc()
  }
  
  if(join == T & BPCells == T){
    message(paste("Converting count matrices"))
    # get vector of merged cells
    cells_merged <- strsplit(as.character(metacell_obj$cells_merged), ",")
    mat_merge <- as(matrix(
      data = sapply(cells_merged, 
                    function(x, mat) rowSums(mat[,x]), mat=mat), 
      nrow = nrow(metacell_obj), 
      ncol = ncol(metacell_obj), 
      dimnames = list(rownames(metacell_obj), 
                      colnames(metacell_obj))), 
      "sparseMatrix")
    gc()
    message(paste0("Writing BPCells Matrix to ", matrix_save_path))
    write_matrix_dir(mat_merge, matrix_save_path)
    mat_merge <- open_matrix_dir(matrix_save_path)
    metacell_obj[[assay]] <- CreateAssay5Object(counts = mat_merge)
    metacell_obj <- UpdateSeuratObject(metacell_obj)
    rm(cells_merged, mat_merge); gc()
  }
  
  if(conv_2_sparse == T){
    message(paste("Converting assay to v3 and count matrix to sparse matrix"))
    mat <- as(open_matrix_dir(matrix_save_path), Class = "sparseMatrix")
    # Make sure dimnames are okay
    rownames(mat) <- rownames(metacell_obj); colnames(mat) <- colnames(metacell_obj)
    metacell_obj[[assay]] <- CreateAssayObject(counts = mat)
    metacell_obj <- UpdateSeuratObject(metacell_obj)
  }
  
  # set idents for metacell object:
  Idents(metacell_obj) <- metacell_obj@meta.data[[ident.group]]
  
  # revert to full seurat object if we subsetted earlier
  if(!is.null(cells.use)){
    seurat_obj <- seurat_full
  }
  
  # add seurat metacell object to the main seurat object:
  seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
  gc()
  
  # add other info
  seurat_obj <- SetWGCNAParams(
    seurat_obj, params = list(
      'metacell_k' = k,
      'metacell_reduction' = reduction,
      'metacell_layer' = layer,
      'metacell_assay' = assay,
      'metacell_stats' = run_stats
    ),
    wgcna_name
  )
  seurat_obj
}



#Breakdown Function
#subset seurat_obj
#cells_use =X
seurat_full <- seurat_obj

# setup grouping variables 
group.by=c('broad_type_leiden', 'sample')


if(length(group.by) > 1){
  seurat_meta <- seurat_obj@meta.data[,group.by]
  for(col in colnames(seurat_meta)){
    seurat_meta[[col]] <- as.character(seurat_meta[[col]])
  }
  seurat_obj$metacell_grouping <- apply(seurat_meta, 1, paste, collapse='#')
} else {
  seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
}
groupings <- unique(seurat_obj$metacell_grouping)
groupings <- groupings[order(groupings)]


# remove groups that are too small:
min_cells = 100
group_counts <- table(seurat_obj$metacell_grouping) < min_cells
if(any(group_counts)){
  warning(paste0("Removing the following groups that did not meet min_cells: ", paste(names(group_counts)[group_counts], collapse=', ')))
}
groupings <- groupings[table(seurat_obj$metacell_grouping) >= min_cells]

if(length(groupings) == 0 ){
  stop("No groups met the min_cells requirement.")
}

# unique meta-data for each group
meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, '#')))
colnames(meta_df) <- group.by

# list of meta-data to pass to each metacell seurat object
meta_list <- lapply(1:nrow(meta_df), function(i){
  x <- list(as.character(meta_df[i,]))[[1]]
  names(x) <- colnames(meta_df)
  x
})


# split seurat obj by groupings
seurat_list <- lapply(groupings, function(x){seurat_obj[,seurat_obj$metacell_grouping == x]})
names(seurat_list) <- groupings

# Adjust k and max_shared parameters for each sample
k=50
k_max = 75
k_min = 30
max_shared_ratio = 0.25

k_list <- numeric()
max_shared_list <- numeric()
cells_per_sample <- numeric()
mean_cells <- numeric()
for(i in 1:length(seurat_list)) {
  cells_per_sample[i] <- ncol(seurat_list[[i]])
}
mean_cells <- mean(cells_per_sample)

for(i in 1:length(cells_per_sample)) {
  k_list[i] <- ifelse(round(k*cells_per_sample[i]/mean_cells) > k_min,
                      round(k*cells_per_sample[i]/mean_cells), 
                      k_min)
  k_list[i] <- ifelse(k_list[i] > k_max, 
                      k_max, 
                      k_list[i])
  max_shared_list[i] <- round(k_list[i]*max_shared_ratio)
}


# Adjust target_metacell and max_iter parameters for each sample

target_metacells=1000
max_iter=5000


target_metacells_list <- numeric()
max_iter_list <- numeric()
for(i in 1:length(cells_per_sample)) {
  if(target_metacells < 1000) {
    target_metacells_list[i] <- max(target_metacells, round((target_metacells/100)*(cells_per_sample[i]/k_list[i])))
    max_iter_list[i] <- target_metacells_list[i]*5
  } else {
    target_metacells_list[i] <- target_metacells
    max_iter_list[i] <- max_iter
  }
}


seurat_obj@assays[["RNA"]]@layers[["counts"]] <- seurat_obj@assays[["RNA"]]@layers[["counts"]]@matrix@matrix

# construct metacells
reduction='harmony.rna'
dims=1:100
assay = "RNA"
layer = "counts"
mat = seurat_obj@assays[["RNA"]]@layers[["counts"]]
BPCells =T
mode = "sum"
verbose = T
wgcna_name = seurat_obj@misc$active_wgcna



metacell_list <- mapply(
  ConstructMetacells,
  seurat_obj = seurat_list,
  name = groupings,
  meta = meta_list,
  k=k_list,
  max_shared=max_shared_list,
  max_iter=max_iter_list, 
  target_metacells=target_metacells_list,
  MoreArgs = list(
    reduction=reduction, 
    dims=dims,
    assay=assay, 
    layer=layer,
    BPCells=BPCells,
    return_metacell=F, 
    mat=mat,
    mode=mode,
    verbose=verbose, 
    wgcna_name=wgcna_name
  )
)
names(metacell_list) <- groupings

rm(seurat_list); gc()

# remove NULL
remove <- which(sapply(metacell_list, is.null))
if(length(remove) >= 1){
  metacell_list <- metacell_list[-remove]
}

# get the run stats:
run_stats <- as.data.frame(do.call(rbind, lapply(metacell_list, function(x){x@misc$run_stats})))
rownames(run_stats) <- 1:nrow(run_stats)
for(i in 1:length(group.by)){
  run_stats[[group.by[i]]] <- do.call(rbind, strsplit(as.character(run_stats$name), '#'))[,i]
}

message(paste("Merging metacells"))

# combine metacell objects
if(length(metacell_list) > 1){
  metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])
  gc()
} else{
  metacell_obj <- metacell_list[[1]]
  gc()
}


join = T
matrix_save_path = "/home/wamalwa@nyspi.local/hdWGCNA/4.22.25/BPCells_Matrix/return_metacell_F"
conv_2_sparse = T
ident.group='cell_group'
cells.use = NULL

if(join == T & BPCells == T){
  message(paste("Converting count matrices"))
  # get vector of merged cells
  cells_merged <- strsplit(as.character(metacell_obj$cells_merged), ",")
  mat_merge <- as(matrix(
    data = sapply(cells_merged, 
                  function(x, mat) rowSums(mat[,x]), mat=mat), 
    nrow = nrow(metacell_obj), 
    ncol = ncol(metacell_obj), 
    dimnames = list(rownames(metacell_obj), 
                    colnames(metacell_obj))), 
    "sparseMatrix")
  gc()
  message(paste0("Writing BPCells Matrix to ", matrix_save_path))
  write_matrix_dir(mat_merge, matrix_save_path)
  mat_merge <- open_matrix_dir(matrix_save_path)
  metacell_obj[[assay]] <- CreateAssay5Object(counts = mat_merge)
  metacell_obj <- UpdateSeuratObject(metacell_obj)
  rm(cells_merged, mat_merge); gc()
}

if(conv_2_sparse == T){
  message(paste("Converting assay to v3 and count matrix to sparse matrix"))
  mat <- as(open_matrix_dir(matrix_save_path), Class = "sparseMatrix")
  # Make sure dimnames are okay
  rownames(mat) <- rownames(metacell_obj); colnames(mat) <- colnames(metacell_obj)
  metacell_obj[[assay]] <- CreateAssayObject(counts = mat)
  metacell_obj <- UpdateSeuratObject(metacell_obj)
}

# set idents for metacell object:
Idents(metacell_obj) <- metacell_obj@meta.data[["broad_type_leiden"]]


# revert to full seurat object if we subsetted earlier
if(!is.null(cells.use)){
  seurat_obj <- seurat_full
}


seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)  
gc()




saveRDS(seurat_obj, "seurat_obj_k_50_4.23.25.rds")


seurat_obj_mapping <- MetacellsByGroups(seurat_obj, group.by=c('broad_type_leiden', 'sample'), ident.group='broad_type_leiden',
                                        mat=seurat_obj@assays[["RNA"]]@layers[["counts"]]@matrix@matrix, BPCells = T, join = T, matrix_save_path = "BPCells/Metacells", conv_2_sparse = T, 
                                        reduction='harmony.rna', dims=1:100,
                                        cells.use = NULL, min_cells=100,
                                        assay="RNA", layer ='counts', mode = 'sum',
                                        k=50, k_max = 75, k_min = 30, max_shared_ratio = 0.25,
                                        target_metacells=1000, max_iter=5000, verbose=T, wgcna_name=NULL)

#########################################################################################################
#seurat_obj <- readRDS("~/hdWGCNA/4.22.25/seurat_obj_k_50_4.23.25.rds")

ScaleMetacells <- function(seurat_obj, wgcna_name=NULL, ...){
  if(!exists("features")){
    features = VariableFeatures(seurat_obj)
  }
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::ScaleData(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}



seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  assay = "RNA",         # Or whatever assay you're using
  layer = "counts",      # Adjust if needed (e.g., "data" or "logcounts")
  selection.method = "vst",
  nfeatures = 2000
)

seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))

seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='sample')

seurat_obj <- RunUMAPMetacells(seurat_obj, reduction= "harmony", dims=1:15)

p1 <- DimPlotMetacells(seurat_obj, group.by='broad_type_leiden', label= T, reduction = "umap") + umap_theme() + ggtitle("Broad Type Leiden")
p2 <- DimPlotMetacells(seurat_obj, group.by='sample', reduction = "umap") + umap_theme() + ggtitle("Sample")

p1 | p2


length(unique(seurat_obj@misc[["tutorial_GN"]][["wgcna_metacell_obj"]]@meta.data[["cells_merged"]]))


library(ggplot2)
seurat_obj@misc[[wgcna_name]][["wgcna_metacell_obj"]]

# Cells in full object
head(colnames(seurat_obj)) 

# Metacell "cells"
head(colnames(seurat_obj@misc[[wgcna_name]][["wgcna_metacell_obj"]]))

########################################################################################################
#Soft Threshold Power #For GN broad leiden 

#i) 


wgcna_obj <- seurat_obj@misc[["tutorial_GN"]][["wgcna_metacell_obj"]]

head(wgcna_obj@meta.data)
table(wgcna_obj@meta.data$broad_type_leiden)  # or whatever your cell type label is called

GN <- subset(wgcna_obj, subset = broad_type_leiden == "GN")

seurat_obj@misc[["tutorial_GN"]][["GN_metacells"]] <- GN

rm(wgcna_obj)
rm(GN)

GN_expr_data <- seurat_obj@misc[["tutorial_GN"]][["GN_metacells"]]@assays[["RNA"]]@data
GC_expr_data_t <- t(GN_expr_data)

TestSoftPowers <- function(
    seurat_obj,
    powers=c(seq(1,10,by=1), seq(12,30, by=2)),
    networkType="signed",
    corFnc='bicor',
    wgcna_name = "tutorial_GN",
    ...
){
  
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  
  # get datExpr
  datExpr <-  GC_expr_data_t
  
  # Call the network topology analysis function 
  powerTable = list(
    data = WGCNA::pickSoftThreshold(
      datExpr,
      powerVector=powers,
      verbose = 100,
      networkType=networkType,
      corFnc=corFnc,
      ...
    )[[2]]
  );
  
  # set the power table in Seurat object:
  seurat_obj <- SetPowerTable(seurat_obj, powerTable$data, wgcna_name=wgcna_name)
  seurat_obj
}

# Test different soft powers:  # tart: 11:45, working in blocks of 1222 genes (30 blocks total)
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)


#length(rownames(seurat_obj))

power_table <- GetPowerTable(seurat_obj)
head(power_table) 

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)


write.csv(seurat_obj@misc[["trial 4.7.25"]][["wgcna_powerTable"]], "plot_list_4.17.25.csv")

saveRDS(seurat_obj, "seurat_obj_GC_softpower_4.17.25.rds") #11.1 GB


