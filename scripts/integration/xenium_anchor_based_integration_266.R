# Anthony Ramnauth
# multiome projection onto older Xenium


library(Seurat) # must be >= 5.1 to LoadXenium()
library(tidyverse)
library(future)
library(Polychrome)
library(grid)

reference <- readRDS("/home/ramnant@nyspi.local/Documents/rna_counts_and_wnn_umap.rds")

# load metadata
metad <- read.csv("/home/ramnant@nyspi.local/Documents/Visium/current_working_multiome_meta.csv")
reference@meta.data <- metad
rownames(reference@meta.data) <- reference@meta.data$X

xen1 <- readRDS('/home/ramnant@nyspi.local/Documents/Visium/xen1_sc_integration_2025_april_leiden.rds')
xen2 <- readRDS('/home/ramnant@nyspi.local/Documents/Visium/xen2_sc_integration_2025_april_leiden.rds')

reference[['RNA']] <- subset(reference[['RNA']], features = rownames(xen1))
reference <- subset(reference, subset = nCount_RNA > 0) # 4 cells removed

# Use SCT for normalization as that is what is used for Xenium
reference <- SCTransform(reference, assay = "RNA", verbose = FALSE)

library(future)
options(future.globals.maxSize = 20000 * 1024^2)

anchors <- FindTransferAnchors(reference = reference, query = xen1, normalization.method = 'SCT')
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = reference$cls_lasso_6b_labeled, 
                                  prediction.assay = T, 
                                  weight.reduction = xen1[["pca"]], 
                                  dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = reference$cls_lasso_6b_labeled, dims = 1:30)
xen1$predicted.id <- predictions$predicted.id

saveRDS(xen1, '/home/ramnant@nyspi.local/Documents/Visium/xen1_multiomeintegration.rds')
# xen1 <- readRDS('/home/ramnant@nyspi.local/Documents/Visium/xen1_multiomeintegration.rds')
write.csv(predictions, '/home/ramnant@nyspi.local/Documents/Visium/xen1_multiome_predictions.csv')
# predictions <- read.csv('/home/ramnant@nyspi.local/Documents/Visium/xen1_multiome_predictions.csv')

#### xen 2 ####
anchors2 <- FindTransferAnchors(reference = reference, query = xen2, normalization.method = 'SCT')
predictions.assay2 <- TransferData(anchorset = anchors2, 
                                   refdata = reference$cls_lasso_6b_labeled, 
                                   prediction.assay = T, 
                                   weight.reduction = xen2[["pca"]], 
                                   dims = 1:30)
predictions2 <- TransferData(anchorset = anchors2, refdata = reference$cls_lasso_6b_labeled, dims = 1:30)
xen2$predicted.id <- predictions2$predicted.id

saveRDS(xen2, '/home/ramnant@nyspi.local/Documents/Visium/xen2_multiomeintegration.rds')
# xen2 <- readRDS('/home/ramnant@nyspi.local/Documents/Visium/xen2_multiomeintegration.rds')
write.csv(predictions2, '/home/ramnant@nyspi.local/Documents/Visium/xen2_multiome_predictions.csv')
# predictions2 <- read.csv('/home/ramnant@nyspi.local/Documents/Visium/xen2_multiome_predictions.csv')

#############
# Plots
#############

# Spatial plots
pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/oldXenium_allMultiome_predictions.pdf', width = 15, height = 15)

cell_colors <- as.vector(Polychrome::palette36.colors(35))
#assign colors to seurat_clusters values
names(cell_colors) <- levels(xen1$predicted.id)

ImageDimPlot(xen1, group.by = "predicted.id", cols = cell_colors, border.size = NA, size = 0.75)

ImageDimPlot(xen2, group.by = "predicted.id", cols = cell_colors, border.size = NA, size = 0.75)

dev.off()

# Plot each spatial and sn cluster separately

cell_colors <- cols <- as.vector(Polychrome::palette36.colors(35))
#assign colors to seurat_clusters values
names(cols) <- levels(xen1$predicted.id)

Idents(xen1) <- xen1$predicted.id
# Plot each single cell type for xen

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/oldXenium_ImGN_predictions.pdf', width = 15, height = 60)

p1 <- ImageDimPlot(xen1, group.by = "predicted.id", cols = cell_colors, border.size = NA, size = 0.75)
p2 <- ImageDimPlot(xen1, size = 2, cols = c("red"), border.size = NA, cells = WhichCells(xen1, idents = "ImGN"))

grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

dev.off()

# plot each individual cluster

# First replace metadata with predictions for plotting prediction values
xen1@meta.data <- predictions
rownames(xen1@meta.data) <- xen1$X
xen2@meta.data <- predictions2
rownames(xen2@meta.data) <- xen2$X

celltypes <- colnames(predictions[3:37])

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/oldXenium_ImGN_predictions.pdf', width = 15, height = 60)

ImageFeaturePlot(xen1, features = 'prediction.score.ImGN', size = 0.75)
ImageFeaturePlot(xen2, features = 'prediction.score.ImGN', size = 0.75)

dev.off()

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/oldXenium_GN1_predictions.pdf', width = 10, height = 15)

ImageFeaturePlot(xen1, features = 'prediction.score.GN.1', border.size = NA)
ImageFeaturePlot(xen2, features = 'prediction.score.GN.1', border.size = NA)

dev.off()

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/xen1_lasso_wnn_predictions.pdf', width = 10, height = 10)

for (s in celltypes) {
  p1 <- ImageFeaturePlot(xen1, features = s, border.size = NA, dark.background = TRUE, axes = TRUE) +
    theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  grid.newpage()
  grid.draw(rbind(ggplotGrob(p1), size = "last"))
}

dev.off()

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/xen2_lasso_wnn_predictions.pdf', width = 10, height = 10)

for (s in celltypes) {
  p1 <- ImageFeaturePlot(xen2, features = s, border.size = NA, dark.background = TRUE, axes = TRUE) +
    theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  grid.newpage()
  grid.draw(rbind(ggplotGrob(p1), size = "last"))
}

dev.off()

