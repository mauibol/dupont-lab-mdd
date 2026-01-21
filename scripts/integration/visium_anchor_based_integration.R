###################################
# Visium MDD paper
# Integrate V1 and V2 Visium with multiome
# Maddy Peng
###################################

# set working directory

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(here)
})

# Read in multiome and visium objs
v1 <- readRDS("/data/Share/Victor_Anosike/V1_Spatial/combined_integrated_processed_V1.rds")
v2 <- readRDS("/data/Share/Victor_Anosike/V2_Second_Attempt_Spatial/combined_integrated_processed_V2.rds")
multi <- readRDS(paste0(MERGED_RDS_PATH, 'multimodal_bpcells.rds'))

################ anchor-based integration v1 #####################
DefaultAssay(v1) <- "integrated"
v1 <- FindVariableFeatures(v1)

# Find anchors
anchors1 <- FindTransferAnchors(reference = multi, reference.assay = 'RNA', query.assay = 'integrated', query = v1)

# transfer learning with multiome
anchor.predictions1 <- TransferData(anchorset = anchors1, refdata = multi$cell_annotation_leiden, prediction.assay = T, weight.reduction = v1[['pca']], dims = 1:30)
v1[['anchor_predictions']] <- anchor.predictions1
DefaultAssay(v1) <- 'anchor_predictions'


# create a column variable for metadata that concatenates donor age and sample
v1$sampleinfo <- paste0(v1$donor, "_", v1$age, "_", v1$sample)


pdf(paste0(PLOTS_PATH, "v1_anchorbased_cell_annotation_leiden_GC.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v1$sample)) {
  p2 <- SpatialFeaturePlot(v1, features = 'GC.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v1_anchorbased_cell_annotation_leiden_Astro.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v1$sample)) {
  p2 <- SpatialFeaturePlot(v1, features = 'Astro.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v1_anchorbased_cell_annotation_leiden_Astro.2.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v1$sample)) {
  p2 <- SpatialFeaturePlot(v1, features = 'Astro.2', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v1_anchorbased_cell_annotation_leiden_Oligo.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v1$sample)) {
  p2 <- SpatialFeaturePlot(v1, features = 'Oligo.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}

dev.off()

pdf(paste0(PLOTS_PATH, "v1_anchorbased_cell_annotation_leiden_Oligo.2.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v1$sample)) {
  p2 <- SpatialFeaturePlot(v1, features = 'Oligo.2', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v1_anchorbased_cell_annotation_leiden_ExN.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v1$sample)) {
  p2 <- SpatialFeaturePlot(v1, features = 'ExN.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

################ anchor-based integration v2 #####################
DefaultAssay(v2) <- "integrated"
v2 <- FindVariableFeatures(v2)

# Find anchors
anchors2 <- FindTransferAnchors(reference = multi, reference.assay = 'RNA', query.assay = 'integrated', query = v2)

# transfer learning with multiome
anchor.predictions2 <- TransferData(anchorset = anchors2, refdata = multi$cell_annotation_leiden, prediction.assay = T, weight.reduction = v2[['pca']], dims = 1:30)
v2[['anchor_predictions']] <- anchor.predictions2
DefaultAssay(v2) <- 'anchor_predictions'


# create a column variable for metadata that concatenates donor age and sample
v2$sampleinfo <- paste0(v2$donor, "_", v2$age, "_", v2$sample)


pdf(paste0(PLOTS_PATH, "v2_anchorbased_cell_annotation_leiden_GC.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v2$sample)) {
  p2 <- SpatialFeaturePlot(v2, features = 'GC.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v2_anchorbased_cell_annotation_leiden_Astro.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v2$sample)) {
  p2 <- SpatialFeaturePlot(v2, features = 'Astro.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v2_anchorbased_cell_annotation_leiden_Astro.2.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v2$sample)) {
  p2 <- SpatialFeaturePlot(v2, features = 'Astro.2', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v2_anchorbased_cell_annotation_leiden_Oligo.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v2$sample)) {
  p2 <- SpatialFeaturePlot(v2, features = 'Oligo.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}

dev.off()

pdf(paste0(PLOTS_PATH, "v2_anchorbased_cell_annotation_leiden_Oligo.2.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v2$sample)) {
  p2 <- SpatialFeaturePlot(v2, features = 'Oligo.2', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()

pdf(paste0(PLOTS_PATH, "v2_anchorbased_cell_annotation_leiden_ExN.1.pdf"), height = 7, width = 7, onefile = T)
for (s in unique(v2$sample)) {
  p2 <- SpatialFeaturePlot(v2, features = 'ExN.1', ncol = 1, pt.size.factor = 2, alpha = c(0.1,1), images = s) +
    ggtitle(s)
  print(p2)
}
dev.off()
