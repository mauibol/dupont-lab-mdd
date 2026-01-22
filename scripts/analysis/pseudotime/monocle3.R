# Pseudotime analysis with Monocle3
# Maddy Peng

library(monocle3)
library(ggplot2)
library(dplyr)


# Monocle3 on Astro and GC subset
expression_bold <- neuro_sub[['RNA']]@layers$counts
cell_metadata <- neuro_sub@meta.data

cds <- new_cell_data_set(expression_bold,
                         cell_metadata = cell_metadata)

cds <- preprocess_cds(cds, num_dim = 10)
cds <- align_cds(cds, alignment_group = "sample")

#Add dim reduction objects from seurat
harmony <- Embeddings(neuro_sub, 'harmony.rna')
reducedDim(cds, "PCA") <- harmony
umap <- Embeddings(imgns, 'wnn.umap.sub')
reducedDim(cds, "UMAP") <- umap

# Plot UMAP
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "trajectory_annotations")

# Monocle cluster into partitions
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

# Learn graph and pick starting point in trajectory
cds <- learn_graph(cds,use_partition = F)
cds <- order_cells(cds)


#### Plots ####
p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=F,
                 label_leaves=F,
                 label_branch_points=T,
                 graph_label_size=1,
                 trajectory_graph_color = 'black')

print(p1)

p2 <- plot_cells(cds,
                 color_cells_by = "trajectory_annotations",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=TRUE,
                 graph_label_size=1,
                 trajectory_graph_color = 'black')

print(p2)


cds_sub <- cds2[c('SOX2','PAX6', 'NES','ASCL1','NEUROD1', 'PCNA', 'MCM2','CALB2', 'RELN', 'DCX',
                  'CALB1', 'TUBB3','BHLHE22', 'FST','POSTN','PROX1','RBFOX3'),]

p3 <- plot_genes_in_pseudotime(cds_sub, label_by_short_name = F, ncol=1)

print(p3)

pdf('/data/Share/pengmad/quest_for_immature/monocle_pseudotime.pdf', width=5, height=4)
print(p2)
dev.off()

pdf('/data/Share/pengmad/quest_for_immature/monocle_by_cluster.pdf', width=5, height=4)
print(p1)
dev.off()

pdf('/data/Share/pengmad/quest_for_immature/monocle_gene_trends.pdf', width=5, height=12)
print(p3)
dev.off()






