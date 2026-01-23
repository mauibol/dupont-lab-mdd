# Install from CRAN or GitHub
# CRAN version (might be older)
install.packages("harmony")

# Or GitHub (latest version)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("immunogenomics/harmony")


library(Seurat)
library(harmony)
library(dplyr)

getwd()
setwd("/data/Share/Mukami/Integration/results_11.7.25_spread_parameter")

#read processed zhou (need to rerun umap)
zhou <- readRDS("/data/Share/Mukami/Integration/results_11.7.25_spread_parameter/zhou_edited_11.7.25.rds")



#zhou <- readRDS('/data/Share/marimad/preprocessing_home/lasso/merged_rpca.rds')
franjic <- readRDS('/data/Share/marimad/sc_datasets/franjic_boldrini_harmony.rds')
ayhan <- readRDS('/data/Share/marimad/sc_datasets/ayhan_boldrini_harmony.rds')
#boldrini <- readRDS('/data/Share/pengmad/merged_CURRENT/merged_rna_all_march2025.rds')
boldrini2 <- readRDS('/data/Share/pengmad/multiome_official/current_working_multimodal_bpcells.rds')

# View available metadata columns
colnames(boldrini2@meta.data)
colnames(zhou@meta.data)

# Check how cell names are structured
head(colnames(boldrini2))
head(colnames(zhou))

#add trajectory information 
trajectory <- read.csv("/data/Share/pengmad/quest_for_immature/imgns_sub_metadata_10.09.25.csv")

head(trajectory)

# make sure rownames match
head(rownames(boldrini2@meta.data))
head(trajectory$X)

# Set rownames of trajectory to the barcodes
rownames(trajectory) <- trajectory$X

# Add trajectory column to Seurat metadata
boldrini2$trajectory <- trajectory[rownames(boldrini2@meta.data), "x"]

# Check it
head(boldrini2@meta.data$trajectory)
table(boldrini2@meta.data$trajectory)
table(boldrini2@meta.data$cell_annotation_leiden_labeled)

# Create a new column: use trajectory if available, otherwise fall back to cell_annotation_leiden_labeled


# Convert to character first to avoid factor issues
trajectory_chr <- as.character(boldrini2@meta.data$trajectory)
cell_annotation_chr <- as.character(boldrini2@meta.data$cell_annotation_leiden_labeled)

# Create combined annotation: use trajectory if not NA/empty, else cell_annotation
boldrini2$combined_annotation <- coalesce(trajectory_chr, cell_annotation_chr)

# Check
head(boldrini2@meta.data$combined_annotation)

#New Boldrini Metadata
table(boldrini2@meta.data$combined_annotation)

table(boldrini2$seq_study_ID). #use this for donor 


saveRDS(boldrini2, "multiome_data_11.10.25_updated_meta.rds")


#Remove Zhong 
# 1. Check how many cells are from each investigator
table(zhou$investigator)

# Get all cells NOT from Zhong
keep_cells <- colnames(zhou)[zhou$investigator != "Zhong"]

# Subset using the vector of cells
zhou <- subset(zhou, cells = keep_cells)

# 3. Optional: check the new cell counts
table(zhou$investigator)

# 4. Normalize RNA assay data again (log-normalization)
zhou <- NormalizeData(zhou, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

# 5. Find variable features
zhou <- FindVariableFeatures(zhou, assay = "RNA", selection.method = "vst", nfeatures = 2000)

# 6. Scale the data
zhou <- ScaleData(zhou, assay = "RNA", features = rownames(zhou))

# Use integrated assay for PCA / UMAP
DefaultAssay(zhou) <- "integrated"

# Run PCA and UMAP if needed
zhou <- RunPCA(zhou, verbose = FALSE)
zhou <- RunUMAP(zhou, dims = 1:30, verbose = FALSE)

library(harmony)

# Harmony on PCA embeddings, correcting for 'investigator'
library(harmony)

zhou <- RunHarmony(
  object = zhou,          # your Seurat object
  group.by.vars = "investigator",  # variable to correct for
  reduction.use = "pca",      # which reduction to use
  assay.use = "RNA",      # which assay to use
  project.dim = FALSE
)


zhou <- RunUMAP(
  zhou,
  reduction = "harmony",
  dims = 1:15,
  reduction.name = "umap.harmony",
  reduction.key = "UMAPHarmony_",
  spread = .75
)



# 9. Check the new reductions
Reductions(zhou)



# Select metadata columns you want to transfer (add or overwrite)
cols_to_transfer <- c("combined_annotation")

# Check which of these columns actually exist in boldrini
cols_to_transfer <- cols_to_transfer[cols_to_transfer %in% colnames(boldrini2@meta.data)]

# Identify common cells between the two objects
common_cells <- intersect(colnames(zhou), colnames(boldrini2))

head(common_cells)
length(common_cells)

# Subset metadata to those cells and columns
metadata_to_add <- boldrini2@meta.data[common_cells, cols_to_transfer, drop = FALSE]

# Add to zhou (this will merge metadata columns by cell name)
zhou <- AddMetaData(zhou, metadata = metadata_to_add)

colnames(zhou@meta.data)

#
#Boldrini Old
table(zhou@meta.data$investigator, zhou@meta.data$cell_annotation)

#Boldrini New
table(zhou@meta.data$investigator, zhou@meta.data$cell_annotation_leiden_labeled)

#Zhou cell types 
table(zhou@meta.data$investigator, zhou@meta.data$cellAnnotation)



# Remove old combined column
if ("cellAnnotation_combined" %in% colnames(zhou@meta.data)) {
  zhou$cellAnnotation_combined <- NULL
}

# Initialize new combined column
zhou$cellAnnotation_combined <- NA_character_

# Boldrini cells (common)
common_cells <- intersect(colnames(zhou), colnames(boldrini2))
boldrini_labels <- setNames(as.character(boldrini2$combined_annotation), colnames(boldrini2))
zhou$cellAnnotation_combined[common_cells] <- paste0("Peng-", boldrini_labels[common_cells])

# Zhou-only cells (non-common)
zhou_only_cells <- setdiff(colnames(zhou), common_cells)

# Only add prefix to non-NA labels
non_na_zhou <- !is.na(zhou$cellAnnotation[zhou_only_cells])
zhou$cellAnnotation_combined[zhou_only_cells[non_na_zhou]] <- paste0("Zhou-", zhou$cellAnnotation[zhou_only_cells[non_na_zhou]])

# Optional: convert to factor
zhou$cellAnnotation_combined <- factor(zhou$cellAnnotation_combined)

# Check result
table(zhou$cellAnnotation_combined)


Reductions(zhou)

getwd()

saveRDS(zhou, "zhou_edited_11.7.25.rds")

# Get details for a specific one (for example, UMAP)
zhou[["umap"]]

DimPlot(zhou, reduction = "umap.harmony", group.by = "cellAnnotation_combined", label = TRUE, repel = TRUE)

# Load required libraries
library(Seurat)
library(ggplot2)
library(ggrepel)


p_umap <- DimPlot(
  zhou,
  reduction = "umap.harmony",
  group.by = "cellAnnotation_combined",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 14, family = "Arial")
  ) +
  ggtitle("Zhou Integrated Dataset: CellAnnotaion rnaHarmony_UMAP")

# Print the plot to the viewer (VS Code plot pane)
print(p_umap_2)

# Create an output directory (if it doesn’t already exist)
out_dir <- "11.10.25/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Zhou_UMAP_cellAnnotation rnaHarmonu_UMAP.pdf"),
  plot = p_umap,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Zhou_UMAP_cellAnnotation Sample rnaHarmonu_UMAP.png"),
  plot = p_umap,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)
#UMAPS for donor and investigator

#Boldrini
table(zhou@meta.data$investigator, zhou@meta.data$ID4)

#Zhou and Bolrini  cell types 
table(zhou@meta.data$investigator, zhou@meta.data$sample)



p_umap_2 <- DimPlot(
  zhou,
  reduction = "umap.harmony",
  group.by = "sample",
  label = FALSE,
  repel = TRUE,
  raster = FALSE
) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 14, family = "Arial")
  ) +
  ggtitle("Zhou Integrated Dataset: UMAP by Sample rnaHarmony_UMAP")

# Print the plot to the viewer (VS Code plot pane)
print(p_umap_2)

# Create an output directory (if it doesn’t already exist)
out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Zhou_UMAP_cell_by Sample rnaHarmonu_UMAP.pdf"),
  plot = p_umap_2,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Zhou_UMAP_cell_by Sample rnaHarmonu_UMAP.png"),
  plot = p_umap_2,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)


#By Investigator 

p_umap_3 <- DimPlot(
  zhou,
  reduction = "umap.harmony",
  group.by = "investigator",
  label = FALSE,
  repel = TRUE,
  raster = FALSE
) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, family = "Arial"),
    legend.text = element_text(size = 25, family = "Arial"),
    legend.title = element_text(size = 14, family = "Arial")
  ) +
  ggtitle("Zhou Integrated Dataset: UMAP by Investigator rnaHarmony_UMAP")
# Print the plot to the viewer (VS Code plot pane)
print(p_umap_3)

# Create an output directory (if it doesn’t already exist)
#out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Zhou_UMAP by Investigator_rnaHarmony_UMAP.pdf"),
  plot = p_umap_3,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Zhou_UMAP by Investigator_rnaHarmony_UMAP.png"),
  plot = p_umap_3,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)




#Franjic Next


# View available metadata columns
colnames(boldrini2@meta.data)
colnames(franjic@meta.data)

# Check how cell names are structured
head(colnames(boldrini2))
head(colnames(franjic))



#Boldrini New
table(boldrini2@meta.data$ID4, boldrini2@meta.data$combined_annotation )

#Franjic cls
table(franjic@meta.data$investigator, franjic@meta.data$cls )

#Franjic ct
table(franjic@meta.data$investigator, franjic@meta.data$ct )

#Boldrini
table(franjic@meta.data$Investigator)
#Franjic
table(franjic@meta.data$investigator)

table(franjic@meta.data$ct)

table(franjic@meta.data$cell_annotation_leiden_labeled)


# Select metadata columns you want to transfer (add or overwrite)
cols_to_transfer <- c("combined_annotation")

# Check which of these columns actually exist in boldrini
cols_to_transfer <- cols_to_transfer[cols_to_transfer %in% colnames(boldrini2@meta.data)]

# Identify common cells between the two objects
common_cells <- intersect(colnames(franjic), colnames(boldrini2))

head(common_cells)
length(common_cells)

# Subset metadata to those cells and columns
metadata_to_add <- boldrini2@meta.data[common_cells, cols_to_transfer, drop = FALSE]

# Add to zhou (this will merge metadata columns by cell name)
franjic <- AddMetaData(franjic, metadata = metadata_to_add)

colnames(franjic@meta.data)

table(franjic$ct)
table(franjic$combined_annotation)

# Make a new column
# assuming your object is a data frame or metadata slot of a Seurat object
# e.g., franjic@meta.data

franjic$new_annotation <- ifelse(
  grepl("^Boldrini", franjic$ct),  # detect Boldrini labels
  franjic$combined_annotation,     # replace with detailed annotation
  franjic$ct                       # otherwise keep original label
)

table(franjic$new_annotation)


franjic$new_annotation <- ifelse(
  grepl("^Franjic", franjic$new_annotation), 
  franjic$new_annotation, 
  paste0("Peng_", franjic$new_annotation)
)

table(franjic$new_annotation)




Reductions(franjic)

franjic <- RunUMAP(
  franjic,
  reduction = "rna.harmony",
  dims = 1:15,
  reduction.name = "rna.harmony.umap2",
  reduction.key = "UMAPHarmony_",
  spread = .75
)


Reductions(franjic)




# Create the UMAP plot
p_umap <- DimPlot(
  franjic,
  reduction = "rna.harmony.umap2",
  group.by = "new_annotation",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "v",
    legend.text = element_text(size = 14, family = "Arial"),
    legend.title = element_text(size = 14, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial")
  ) +
  guides(
    fill = guide_legend(ncol = 2)
  ) +
  ggtitle("Franjic Integrated Dataset: UMAP CellAnnotation rnaHarmonyUMAP")


# Print the plot to the viewer (VS Code plot pane)
print(p_umap)

# Create an output directory (if it doesn’t already exist)
out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Franjic_UMAP_cellAnnotation_rnaHarmony_UMAP.pdf"),
  plot = p_umap,
  width = 30,
  height = 30,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Franjic_UMAP_cellAnnotation_rnaHarmony_UMAP.png"),
  plot = p_umap,
  width = 30,
  height = 30,
  units = "in",
  dpi = 300
)


#UMAPS for donor and investigator

#Boldrini
table(franjic@meta.data$investigator)
table(franjic@meta.data$Investigator)

# Combine the investigator columns
franjic@meta.data$investigator_combined <- ifelse(
  !is.na(franjic@meta.data$investigator),
  as.character(franjic@meta.data$investigator),
  as.character(franjic@meta.data$Investigator)
)

# Convert to factor if needed
franjic@meta.data$investigator_combined <- factor(franjic@meta.data$investigator_combined)

# Check the result
table(franjic@meta.data$investigator_combined)

#Investigator

p_umap_2 <- DimPlot(
  franjic,
  reduction = "rna.harmony.umap2",
  group.by = "investigator_combined",
  label = F,
  repel = TRUE,
  raster = FALSE  # Set to TRUE for faster rendering, FALSE for publication quality
) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Franjic Integrated Dataset: UMAP by Investigator rnaHarmony_UMAP")

# Print the plot to the viewer (VS Code plot pane)
print(p_umap_2)

# Create an output directory (if it doesn’t already exist)
out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Franjic_UMAP_cell_by Investigator rnaHarmonu_UMAP.pdf"),
  plot = p_umap_2,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Franjic_UMAP_cell_by Investigator rnaHarmonu_UMAP.png"),
  plot = p_umap_2,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)


#By Sample

colnames(franjic@meta.data)
table(franjic@meta.data$sample)
table(franjic@meta.data$samplename)

# Combine the investigator columns
franjic@meta.data$sample_combined <- ifelse(
  !is.na(franjic@meta.data$sample),
  as.character(franjic@meta.data$sample),
  as.character(franjic@meta.data$samplename)
)

# Convert to factor if needed
franjic@meta.data$sample_combined <- factor(franjic@meta.data$sample_combined)

# Check the result
table(franjic@meta.data$sample_combined)


p_umap_3 <- DimPlot(
  franjic,
  reduction = "rna.harmony.umap2",
  group.by = "sample_combined",
  label = FALSE,
  repel = TRUE,
  raster = FALSE
) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "v",
    legend.text = element_text(size = 14, family = "Arial"),
    legend.title = element_text(size = 25, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial")
  ) +
  guides(
    fill = guide_legend(ncol = 2)
  ) +
  ggtitle("Franjic Integrated Dataset: UMAP by Sample rnaHarmonyUMAP")



# Print the plot to the viewer (VS Code plot pane)
print(p_umap_3)

# Create an output directory (if it doesn’t already exist)
#out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Franjic_UMAP_sample_rnaHarmony_UMAP.pdf"),
  plot = p_umap_3,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Franjic_UMAP_sample_rnaHarmony_UMAP.png"),
  plot = p_umap_3,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)



#Ayhan Next 



head(ayhan@meta.data)

table(ayhan@meta.data$ct)


# Select metadata columns you want to transfer (add or overwrite)
cols_to_transfer <- c("combined_annotation")

# Check which of these columns actually exist in boldrini
cols_to_transfer <- cols_to_transfer[cols_to_transfer %in% colnames(boldrini2@meta.data)]

# Identify common cells between the two objects
common_cells <- intersect(colnames(ayhan), colnames(boldrini2))

head(common_cells)
length(common_cells)

# Subset metadata to those cells and columns
metadata_to_add <- boldrini2@meta.data[common_cells, cols_to_transfer, drop = FALSE]

# Add to zhou (this will merge metadata columns by cell name)
ayhan <- AddMetaData(ayhan, metadata = metadata_to_add)

colnames(ayhan@meta.data)

table(ayhan$ct)
table(franjic$combined_annotation)

# Make a new column
# assuming your object is a data frame or metadata slot of a Seurat object
# e.g., franjic@meta.data

ayhan@meta.data$new_annotation <- NULL


ayhan$new_annotation <- ifelse(
  grepl("^Boldrini", ayhan$ct),  # detect Boldrini labels
  ayhan$combined_annotation,     # replace with detailed annotation
  ayhan$ct                       # otherwise keep original label
)

table(ayhan$new_annotation)


ayhan$new_annotation <- ifelse(
  grepl("^Ayhan", ayhan$new_annotation), 
  ayhan$new_annotation, 
  paste0("Peng_", ayhan$new_annotation)
)

table(ayhan$new_annotation)




Reductions(ayhan)

ayhan <- RunUMAP(
  ayhan,
  reduction = "rna.harmony",
  dims = 1:15,
  reduction.name = "rna.harmony.umap2",
  reduction.key = "UMAPHarmony_",
  spread = .75
)


Reductions(ayhan)


# Create the UMAP plot
p_umap <- DimPlot(
  ayhan,
  reduction = "rna.harmony.umap2",
  group.by = "new_annotation",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "v",
    legend.text = element_text(size = 14, family = "Arial"),
    legend.title = element_text(size = 14, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial")
  ) +
  guides(
    fill = guide_legend(ncol = 2)
  ) +
  ggtitle("Ayhan Integrated Dataset: UMAP CellAnnotation rnaHarmonyUMAP")


# Print the plot to the viewer (VS Code plot pane)
print(p_umap)

# Create an output directory (if it doesn’t already exist)
#out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Ayhan_UMAP_cellAnnotation_rnaHarmony_UMAP.pdf"),
  plot = p_umap,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Ayhan_UMAP_cellAnnotation_rnaHarmony_UMAP.png"),
  plot = p_umap,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)


#UMAPS for donor and investigator

#Boldrini
table(ayhan@meta.data$investigator)
table(ayhan@meta.data$Investigator)

# Combine the investigator columns
ayhan@meta.data$investigator_combined <- ifelse(
  !is.na(ayhan@meta.data$investigator),
  as.character(ayhan@meta.data$investigator),
  as.character(ayhan@meta.data$Investigator)
)

# Convert to factor if needed
ayhan@meta.data$investigator_combined <- factor(ayhan@meta.data$investigator_combined)

# Check the result
table(ayhan@meta.data$investigator_combined)

#Investigator

p_umap_2 <- DimPlot(
  ayhan,
  reduction = "rna.harmony.umap2",
  group.by = "investigator_combined",
  label = F,
  repel = TRUE,
  raster = FALSE  # Set to TRUE for faster rendering, FALSE for publication quality
) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Ayhan Integrated Dataset: UMAP by Investigator rnaHarmony_UMAP")

# Print the plot to the viewer (VS Code plot pane)
print(p_umap_2)

# Create an output directory (if it doesn’t already exist)
#out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Ayhan_UMAP_cell_by Investigator rnaHarmonu_UMAP.pdf"),
  plot = p_umap_2,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Ayhan_UMAP_cell_by Investigator rnaHarmonu_UMAP.png"),
  plot = p_umap_2,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)


#By Sample

colnames(ayhan@meta.data)
table(ayhan@meta.data$sample)
table(ayhan@meta.data$id)

# Combine the investigator columns
ayhan@meta.data$sample_combined <- ifelse(
  !is.na(ayhan@meta.data$sample),
  as.character(ayhan@meta.data$sample),
  as.character(ayhan@meta.data$id)
)

# Convert to factor if needed
ayhan@meta.data$sample_combined <- factor(ayhan@meta.data$sample_combined)

# Check the result
table(ayhan@meta.data$sample_combined)


p_umap_3 <- DimPlot(
  ayhan,
  reduction = "rna.harmony.umap2",
  group.by = "sample_combined",
  label = FALSE,
  repel = TRUE,
  raster = FALSE
) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "v",
    legend.text = element_text(size = 14, family = "Arial"),
    legend.title = element_text(size = 25, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial")
  ) +
  guides(
    fill = guide_legend(ncol = 2)
  ) +
  ggtitle("Ayhan Integrated Dataset: UMAP by Sample rnaHarmonyUMAP")



# Print the plot to the viewer (VS Code plot pane)
print(p_umap_3)

# Create an output directory (if it doesn’t already exist)
#out_dir <- "results/plots"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save as PDF (vector)
ggsave(
  filename = file.path(out_dir, "Ayhan_UMAP_sample_rnaHarmony_UMAP.pdf"),
  plot = p_umap_3,
  width = 20,
  height = 20,
  units = "in"
)

# Save as PNG (raster)
ggsave(
  filename = file.path(out_dir, "Ayhan_UMAP_sample_rnaHarmony_UMAP.png"),
  plot = p_umap_3,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300
)



saveRDS(franjic, "franjic_integrated_11.10.25.rds")
saveRDS(ayhan, "ayhan_integrated_11.10.25.rds")
saveRDS(boldrini2, "boldrini_multiome_data_11.10.25.rds")