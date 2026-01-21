# Anthony Ramnauth
# multiome projection onto Xenium 5k

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

xen3 <- readRDS('/data/Share/marimad/Xenium/xen3_integration.rds')

reference[['RNA']] <- subset(reference[['RNA']], features = rownames(xen3))

# Use SCT for normalization as that is what is used for Xenium
reference <- SCTransform(reference, assay = "RNA", verbose = FALSE)

anchors <- FindTransferAnchors(reference = reference, query = xen3, normalization.method = 'SCT')
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = reference$cls_lasso_6b_labeled, 
                                  prediction.assay = T, 
                                  weight.reduction = xen3[["pca"]], 
                                  dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = reference$cls_lasso_6b_labeled, dims = 1:30)
xen3$predicted.id <- predictions$predicted.id

saveRDS(xen3, '/home/ramnant@nyspi.local/Documents/Visium/xen3_multiomeintegration.rds')
# xen3 <- readRDS('/home/ramnant@nyspi.local/Documents/Visium/xen3_multiomeintegration.rds')
write.csv(predictions, '/home/ramnant@nyspi.local/Documents/Visium/xen3_multiome_predictions.csv')
# predictions <- read.csv('/home/ramnant@nyspi.local/Documents/Visium/xen3_multiome_predictions.csv')

#### xen 4 ####
xen4 <- readRDS('/data/Share/marimad/Xenium/xen4_integration.rds')

anchors4 <- FindTransferAnchors(reference = reference, query = xen4, normalization.method = 'SCT')
predictions.assay4 <- TransferData(anchorset = anchors4, 
                                   refdata = reference$cls_lasso_6b_labeled, 
                                   prediction.assay = T, 
                                   weight.reduction = xen4[["pca"]], 
                                   dims = 1:30)
predictions4 <- TransferData(anchorset = anchors4, refdata = reference$cls_lasso_6b_labeled, dims = 1:30)
xen4$predicted.id <- predictions4$predicted.id

saveRDS(xen4, '/home/ramnant@nyspi.local/Documents/Visium/xen4_multiomeintegration.rds')
# xen4 <- readRDS('/home/ramnant@nyspi.local/Documents/Visium/xen4_multiomeintegration.rds')
write.csv(predictions4, '/home/ramnant@nyspi.local/Documents/Visium/xen4_multiome_predictions.csv')
# predictions4 <- read.csv('/home/ramnant@nyspi.local/Documents/Visium/xen4_multiome_predictions.csv')


#############
# Plots
#############

# Spatial plots
pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/allMultiome_predictions.pdf', width = 15, height = 15)

cell_colors <- as.vector(Polychrome::palette36.colors(35))
#assign colors to seurat_clusters values
names(cell_colors) <- levels(xen3$predicted.id)

ImageDimPlot(xen3, group.by = "predicted.id", cols = cell_colors, border.size = NA, size = 0.75)

ImageDimPlot(xen4, group.by = "predicted.id", cols = cell_colors, border.size = NA, size = 0.75)

dev.off()

# Plot each spatial and sn cluster separately

cell_colors <- cols <- as.vector(Polychrome::palette36.colors(35))
#assign colors to seurat_clusters values
names(cols) <- levels(xen4$predicted.id)

Idents(xen4) <- xen4$predicted.id
# Plot each single cell type for xen

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/ImGN_predictions.pdf', width = 15, height = 60)

p1 <- ImageDimPlot(xen4, group.by = "predicted.id", cols = cell_colors, border.size = NA, size = 0.75)
p2 <- ImageDimPlot(xen4, size = 2, cols = c("red"), border.size = NA, cells = WhichCells(xen4, idents = "ImGN"))

grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

dev.off()

# plot each individual cluster

# First replace metadata with predictions for plotting prediction values
xen3@meta.data <- predictions
rownames(xen3@meta.data) <- xen3$X
xen4@meta.data <- predictions4
rownames(xen4@meta.data) <- xen4$X

celltypes <- colnames(predictions[3:37])

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/ImGN_predictions.pdf', width = 15, height = 60)

ImageFeaturePlot(xen3, features = 'prediction.score.ImGN', size = 0.75)
ImageFeaturePlot(xen4, features = 'prediction.score.ImGN', size = 0.75)

dev.off()

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/GN1_predictions.pdf', width = 15, height = 60)

ImageFeaturePlot(xen3, features = 'prediction.score.GN.1', size = 0.75)
ImageFeaturePlot(xen4, features = 'prediction.score.GN.1', size = 0.75)

dev.off()

# plot each cryosection separately to adjust scale, size, etc.

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/xen3_lasso_wnn_predictions.pdf', width = 5, height = 5)

for (s in celltypes) {
  p1 <- ImageFeaturePlot(xen3, features = s, border.size = NA, dark.background = TRUE, axes = TRUE) +
    theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  grid.newpage()
  grid.draw(rbind(ggplotGrob(p1), size = "last"))
}

dev.off()

pdf('/home/ramnant@nyspi.local/Documents/Visium/latest_anchor_predictions/Xenium_Multiome/xen4_lasso_wnn_predictions.pdf', width = 7, height = 7)

for (s in celltypes) {
  p1 <- ImageFeaturePlot(xen4, features = s, border.size = NA, dark.background = TRUE, axes = TRUE) +
    theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  grid.newpage()
  grid.draw(rbind(ggplotGrob(p1), size = "last"))
}

dev.off()

