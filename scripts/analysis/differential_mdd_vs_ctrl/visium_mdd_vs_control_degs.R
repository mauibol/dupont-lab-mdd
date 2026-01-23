######################################
# Visium v1 and v2 MDD vs Control Differential Gene Expression
# By Victor Anosike
# Adapted from code By Maddy Peng
######################################

suppressPackageStartupMessages({
  library(edgeR)
  library(scater)
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(Matrix.utils)
  library(dplyr)
  library(magrittr)
  library(Matrix)
  library(purrr)
  library(reshape2)
  library(S4Vectors)
  library(tibble)
  library(SingleCellExperiment)
  library(pheatmap)
  library(apeglm)
  library(png)
  library(DESeq2)
  library(RColorBrewer)
  library(scran)
  library(ggpubr)
  library(ggrepel)
  library(GWASTools)
  library(rio)
  library(ggbreak)
  library(ggrepel)
  library(car)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(openxlsx)
})


######################################### Visium v1 #####################################################
#########################################################################################################

#########Getting the DEG Data

# Load data and metadata
v1 <- readRDS("insert/path/to/v1.rds")

#Add sex columns to object
v1@meta.data$sex <- rep("m", 65133) 


# Define the function for cluster naming
cluster_namer <- function(i) {
  i <- as.character(i)  # ensure numeric input works
  if (i %in% names(v1_cluster_fullnames)) {
    return(v1_cluster_fullnames[i])
  } else {
    warning(paste("Cluster", i, "not found in v1 cluster list."))
    return(NA)
  }
}

#Reset Cluster Identities back to their original numbers
Idents(v1) <- v1$seurat_clusters

#Create v1_counts
v1_counts <- v1@assays[["Spatial"]]@counts
v1_meta <- v1@meta.data

# Create single cell experiment
v1_sce <- SingleCellExperiment(assays = list(counts = v1_counts), 
                               colData = v1_meta)
dim(counts(v1_sce))

# Change cls_lasso_6b_labeled to the name of the cell type column you are using
v1_sce$sample <- as.factor(v1_sce$sample)



v1_pseudo <- aggregateAcrossCells(v1_sce, 
                                  DataFrame(
                                    cluster = v1_meta[["seurat_clusters"]], #change to your col name
                                    sample = v1_sce$sample
                                  ))

# filter out ncells less than 5
v1_pseudo <- v1_pseudo[,v1_pseudo$ncells >= 5]

# winsorize PMI, RIN, & ncells
# find values that were outside the sum of the third quartile and the product of 1.5 times the inter-quartile range from the median
# Compute the median, IQR, and thresholds
Q3 <- quantile(v1_meta$PMI, 0.75)  # Third quartile
IQR_value <- IQR(v1_meta$PMI)  # Interquartile range
threshold_PMI <- Q3 + 1.5 * IQR_value  # Upper threshold

# Compute the median, IQR, and thresholds
Q1 <- quantile(v1_meta$RIN, 0.25)  # Third quartile
IQR_value <- IQR(v1_meta$RIN)  # Interquartile range
threshold_RIN <- Q1 - 1.5 * IQR_value  # Upper threshold

# Update the pseudobulk colData
colData(v1_pseudo)$RIN_w <- dplyr::case_when(
  colData(v1_pseudo)$RIN < 6.8 ~ 6.8,
  .default =colData(v1_pseudo)$RIN
)

colData(v1_pseudo)$PMI_w <- dplyr::case_when(
  colData(v1_pseudo)$PMI > 47.8 ~ 47.8,
  .default = colData(v1_pseudo)$PMI
)


colData(v1_pseudo)$diagnosis <- factor(colData(v1_pseudo)$Psych1A, levels=c('control', 'MDD'))
#saveRDS(pseudo, paste0(METADATA_PATH, 'pseudo_leiden.rds'))

# Create a data frame with the sample IDs, cluster IDs and condition. Use whatever the column name is in your object
# - All CM samples for V1 belong to one batch, and all MM samples to another batch. That means that we can use
# - batch as a covariable, but since there are only two experimenters, the experimenter column in the V1 metadata effectively acts as 
# - a batch covariable for V1. However, this is not the case for V2, since V2 has six different batches. Therefore, V1 will simply 
# - have the experimenter column's NAME replaced to "batch", but V2 will have the experimenter column removed as a covariable ENTIRELY 
# - and replaced with the variable batch. 
colData(v1_pseudo)$batch <- as.factor(as.character(substr(colData(v1_pseudo)$sample, 1, 2)))
#colData(v1_pseudo)$donor <- as.factor(colData(v1_pseudo)$sample)


# Loop through the 16 clusters (clusters 0 to 15)
cls <- unique(colData(v1_pseudo)$cluster)
cls

v1_degs <- data.frame()
for (type in cls) {
  print(type)
  
  
  
  #if (type %in% unique(v1_degs$cell_type)) next
  
  pse <- v1_pseudo[, v1_pseudo$cluster == type]
  meta <- colData(pse)[, c("sample", "diagnosis", "age", "ncells", "donor",
                           "RIN_w", "PMI_w", "cluster", "batch")]
  
  ###Add this at the start of your for loop -skip empty clusters cleanly
  if (ncol(assays(pse)$counts) == 0 || nrow(assays(pse)$counts) == 0) {
    message("Skipping cluster ", type, " (no counts present)")
    next
  }
  #####
  
  if (length(unique(meta$diagnosis)) == 1) {
    print(paste('Skipping', type, 'because only 1 diagnosis in this observation'))
    next
  }
  
  dge <- DGEList(counts = assays(pse)$counts, group = meta$diagnosis)
  print(dim(dge))
  
  keep.exprs <- filterByExpr(dge, group = 'group')
  dge <- dge[keep.exprs,, keep.lib.sizes = FALSE]
  genes <- rownames(dge)
  dge <- calcNormFactors(dge)
  print(dim(dge))
  
  # Model design with conditional batch handling
  design <- model.matrix(~diagnosis + scale(age) + RIN_w + PMI_w + batch + ncells, meta)
  
  # weights estimated without the correlation
  vobj <- voomWithQualityWeights(dge, design, plot = FALSE)
  dupcor <- duplicateCorrelation(vobj, design, block = meta$donor)
  print(dupcor$consensus.correlation)
  
  # this step uses the genome-wide average for random effect
  # Fit model
  print('fitting linear model...')
  fitDupCor <- tryCatch({
    fit <- lmFit(vobj, design, block = meta$donor, correlation = dupcor$consensus)
    fit <- eBayes(fit)
    topTable_dupcor <- topTable(fit, number = length(genes), coef = 'diagnosisMDD')
    topTable_dupcor$cell_type <- type
    topTable_dupcor$gene <- rownames(topTable_dupcor)
    v1_degs <- rbind(v1_degs, topTable_dupcor)
  }, error = function(e) {
    message("Model failed for ", type, ": ", e$message)
    return(NULL)
  })
}

# Save DEGS unfiltered by p-value to use for MDD vs Control GSEA
write.csv(v1_degs,'/path/to/visium_v1_degs_unfiltered_date.csv', row.names = T)


#########Plotting Dot Plot


# plot the dotplot 
v1_degs <- read.csv('/path/to/visium_v1_degs_unfiltered_date.csv', row.names = 1)
v1_degs_sig <- v1_degs %>% dplyr::filter(adj.P.Val < 0.05)

# Ensure cell_type is a factor with levels 0 to 15 (since there are clusters from 0 to 15
#v1_degs_sig$cell_type <- factor(v1_degs_sig$cell_type, levels = as.character(0:14))
v1_degs_sig$cell_type <- factor(v1_degs_sig$cell_type, levels = as.character(levels(Idents(v1))))


p1 <- v1_degs_sig %>%
  dplyr::mutate(cell_type = as.numeric(as.character(cell_type))) %>%
  dplyr::filter(cell_type >= 0 & cell_type <= 12) %>% #Only doing clusters 0 to 12 since those are the clusters we are measuring
  ggplot(aes(x = logFC, 
             y = factor(sapply(cell_type, cluster_namer), 
                        levels = unique(sapply(sort(unique(cell_type)), cluster_namer))), 
             color = adj.P.Val)) +
  geom_point(alpha = 0.5) +
  scale_color_distiller(palette = "Blues") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  ggtitle('DEGS MDD vs. CTRL Visium V1')

# PDF to display the dotplot
pdf('/path/to/v1_DEG_dotplot_date.pdf', width=4, height=3.5)
print(p1)
dev.off()

############ Plotting Volcano Plots

#Plotting and Posting onto PDF 
pdf('/path/to/v1_all_volcano_date.pdf', width=6, height=5)


#Changed v1 cluster names to match final decided names
v1_cluster_fullnames <- c(
  "0" = "ecm",
  "1" = "luc-rad",
  "2" = "axon",
  "3" = "sgz-ml",
  "4" = "ca1-4.1",
  "5" = "ca1-4.2",
  "6" = "dendr",
  "7" = "sgz-pl",
  "8" = "gcl",
  "9" = "ca1-2",
  "10" = "inn",
  "11" = "vasc",
  "12" = "cp"
)

for (i in 0:12) { #Had to only show 0 to 12 because those are the clusters we are measuring
  
  z <- v1_degs %>% dplyr::filter(cell_type == i)
  
  p <- EnhancedVolcano(z,
                       lab = z$gene,
                       x = 'logFC',
                       y = 'adj.P.Val',
                       FCcutoff = 1,
                       pCutoff = 0.05,
                       labSize = 3.0,
                       drawConnectors = TRUE,
                       arrowheads = FALSE,
                       max.overlaps = 20,
                       ylab = "-log10 Adjusted P-value",
                       legendPosition = 'top',
                       title = paste0("DEGs for ", cluster_namer(i)),
                       subtitle = "MDD vs. CTRL",
                       ylim = c(0,8)) +
    theme_classic() +
    theme(legend.position = "none") +
    #xlim(c(-10,10)) # in order to include more genes
    xlim(c(-5, 5)) 
  #ylim(c(0, 3))
  print(p)
  
}

dev.off()


######Saving Excel Sheets

# Save as a tab separated excel. Make sure to save it to your own path.
df_sort <- v1_degs_sig %>% arrange(adj.P.Val, desc(logFC),.by_group = T)
df_list <- split(df_sort, as.factor(df_sort$cell_type))
write.xlsx(df_list,
           file = '/path/to/v1_degs_significant_no_sex_date.xlsx',
           asTable = TRUE)


####################################### Visium v2 #####################################################
#######################################################################################################

# Define the function
cluster_namer <- function(i) {
  i <- as.character(i)  # ensure numeric input works
  if (i %in% names(v2_cluster_fullnames)) {
    return(v2_cluster_fullnames[i])
  } else {
    warning(paste("Cluster", i, "not found in v1 cluster list."))
    return(NA)
  }
}

#########Getting the DEG Data

# Load data and metadata
v2 <- readRDS("insert/path/to/v2.rds")

# Attach Row Names to Each of the Layers

rownames(v2@assays[["Spatial"]]@layers[["counts.1"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.2"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.3"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.4"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.5"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.6"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.7"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.8"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.9"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.10"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.11"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.12"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.13"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.14"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.15"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.16"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.17"]]) <- rownames(v2@assays[["Spatial"]])
rownames(v2@assays[["Spatial"]]@layers[["counts.18"]]) <- rownames(v2@assays[["Spatial"]])

#Create v2_counts
v2_counts <- cbind(v2@assays[["Spatial"]]@layers[["counts.1"]],
                   v2@assays[["Spatial"]]@layers[["counts.2"]],
                   v2@assays[["Spatial"]]@layers[["counts.3"]],
                   v2@assays[["Spatial"]]@layers[["counts.4"]],
                   v2@assays[["Spatial"]]@layers[["counts.5"]],
                   v2@assays[["Spatial"]]@layers[["counts.6"]],
                   v2@assays[["Spatial"]]@layers[["counts.7"]],
                   v2@assays[["Spatial"]]@layers[["counts.8"]],
                   v2@assays[["Spatial"]]@layers[["counts.9"]],
                   v2@assays[["Spatial"]]@layers[["counts.10"]],
                   v2@assays[["Spatial"]]@layers[["counts.11"]],
                   v2@assays[["Spatial"]]@layers[["counts.12"]],
                   v2@assays[["Spatial"]]@layers[["counts.13"]],
                   v2@assays[["Spatial"]]@layers[["counts.14"]],
                   v2@assays[["Spatial"]]@layers[["counts.15"]],
                   v2@assays[["Spatial"]]@layers[["counts.16"]],
                   v2@assays[["Spatial"]]@layers[["counts.17"]],
                   v2@assays[["Spatial"]]@layers[["counts.18"]])




v2_meta <- v2@meta.data

# Create single cell experiment
v2_sce <- SingleCellExperiment(assays = list(counts = v2_counts), 
                               colData = v2_meta)
dim(counts(v2_sce))

# Change cls_lasso_6b_labeled to the name of the cell type column you are using
#v2_sce$anatomy <- as.factor(v2_sce$anatomy)
v2_sce$sample <- as.factor(v2_sce$sample)



v2_pseudo <- aggregateAcrossCells(v2_sce, 
                                  DataFrame(
                                    cluster = v2_meta[["seurat_clusters"]], #change to your col name
                                    sample = v2_sce$sample
                                  ))

# filter out ncells less than 5
v2_pseudo <- v2_pseudo[,v2_pseudo$ncells >= 5]

# winsorize PMI, RIN, & ncells
# find values that were outside the sum of the third quartile and the product of 1.5 times the inter-quartile range from the median
# Compute the median, IQR, and thresholds
Q3 <- quantile(v2_meta$PMI, 0.75)  # Third quartile
IQR_value <- IQR(v2_meta$PMI)  # Interquartile range
threshold_PMI <- Q3 + 1.5 * IQR_value  # Upper threshold

# Compute the median, IQR, and thresholds
Q1 <- quantile(v2_meta$RIN, 0.25)  # Third quartile
IQR_value <- IQR(v2_meta$RIN)  # Interquartile range
threshold_RIN <- Q1 - 1.5 * IQR_value  # Upper threshold

# Update the pseudobulk colData
colData(v2_pseudo)$RIN_w <- dplyr::case_when(
  colData(v2_pseudo)$RIN < 6.8 ~ 6.8,
  .default =colData(v2_pseudo)$RIN
)

colData(v2_pseudo)$PMI_w <- dplyr::case_when(
  colData(v2_pseudo)$PMI > 47.8 ~ 47.8,
  .default = colData(v2_pseudo)$PMI
)


colData(v2_pseudo)$diagnosis <- factor(colData(v2_pseudo)$diagnosis, levels=c('CTRL', 'MDD'))
#saveRDS(pseudo, paste0(METADATA_PATH, 'pseudo_leiden.rds'))

# Create a data frame with the sample IDs, cluster IDs and condition. Use whatever the column name is in your object
#colData(v2_pseudo)$experimenter <- as.factor(as.character(substr(colData(v2_pseudo)$sample, 1, 2)))
#colData(v2_pseudo)$donor <- as.factor(colData(v2_pseudo)$sample)


# Loop through the 15 clusters
cls <- unique(colData(v2_pseudo)$cluster)
cls

v2_degs <- data.frame()
for (type in cls) {
  print(type)
  
  
  
  #if (type %in% unique(v2_degs$cell_type)) next
  
  pse <- v2_pseudo[, v2_pseudo$cluster == type]
  meta <- colData(pse)[, c("sample", "diagnosis", "age", "ncells", "donor",
                           "RIN_w", "PMI_w", "cluster", "batch")]
  
  ###Add this at the start of your for loop -skip empty clusters cleanly
  if (ncol(assays(pse)$counts) == 0 || nrow(assays(pse)$counts) == 0) {
    message("Skipping cluster ", type, " (no counts present)")
    next
  }
  #####
  
  if (length(unique(meta$diagnosis)) == 1) {
    print(paste('Skipping', type, 'because only 1 diagnosis in this observation'))
    next
  }
  
  dge <- DGEList(counts = assays(pse)$counts, group = meta$diagnosis)
  print(dim(dge))
  
  keep.exprs <- filterByExpr(dge, group = 'group')
  dge <- dge[keep.exprs,, keep.lib.sizes = FALSE]
  genes <- rownames(dge)
  dge <- calcNormFactors(dge)
  print(dim(dge))
  
  # Model design with conditional conditional batch handling
  design <- model.matrix(~diagnosis + scale(age) + RIN_w + PMI_w + batch + ncells, meta)
  
  # weights estimated without the correlation
  vobj <- voomWithQualityWeights(dge, design, plot = FALSE)
  dupcor <- duplicateCorrelation(vobj, design, block = meta$donor)
  print(dupcor$consensus.correlation)
  
  # this step uses the genome-wide average for random effect
  # Fit model
  print('fitting linear model...')
  fitDupCor <- tryCatch({
    fit <- lmFit(vobj, design, block = meta$donor, correlation = dupcor$consensus)
    fit <- eBayes(fit)
    topTable_dupcor <- topTable(fit, number = length(genes), coef = 'diagnosisMDD')
    topTable_dupcor$cell_type <- type
    topTable_dupcor$gene <- rownames(topTable_dupcor)
    v2_degs <- rbind(v2_degs, topTable_dupcor)
  }, error = function(e) {
    message("Model failed for ", type, ": ", e$message)
    return(NULL)
  })
}

# Save DEGS unfilitered by p-value to use for MDD vs Control GSEA
write.csv(v2_degs,'/path/to/visium_v2_degs_unfiltered_date.csv', row.names = T)


#########Plotting

# plot the dotplot 
v2_degs <- read.csv('/path/to/visium_v2_degs_unfiltered_date.csv', row.names = 1)
v2_degs_sig <- v2_degs %>% dplyr::filter(adj.P.Val < 0.05)

# Ensure cell_type is a factor with levels 0 to 18 (since there are clusters from 0 to 18, will have to change
v2_degs_sig$cell_type <- factor(v2_degs_sig$cell_type, levels = as.character(0:18))

v2_cluster_fullnames <- c(
  "0"  = "axon",
  "1"  = "ca1",
  "2"  = "dendr",
  "3"  = "ca3-4",
  "4"  = "luc",
  "5"  = "sgz-ml",
  "6"  = "inn",
  "7"  = "sub",
  "8"  = "or",
  "9"  = "rad",
  "10" = "vasc",
  "11" = "gcl",
  "12" = "cp",
  "13" = "sgz-pl",
  "14" = "ca1-rost",
  "15" = "ca2",
  "16" = "cr"
)


p1 <- v2_degs_sig %>%
  dplyr::mutate(cell_type = as.numeric(as.character(cell_type))) %>%
  dplyr::filter(cell_type >= 0 & cell_type <= 16) %>% # Only doing clusters 0 to 16, since the other clusters are not being measured
  ggplot(aes(x = logFC, 
             y = factor(sapply(cell_type, cluster_namer), 
                        levels = unique(sapply(sort(unique(cell_type)), cluster_namer))), 
             color = adj.P.Val)) +
  geom_point(alpha = 0.5) +
  scale_color_distiller(palette = "Blues") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  ggtitle('DEGS MDD vs. CTRL Visium V2')


pdf('/path/to/v2_DEG_dotplot_date.pdf', width=4, height=3.5)
print(p1)
dev.off()

######Volcano Plots

#Plotting and Posting onto PDF 
pdf('/path/to/v2_all_volcano_date.pdf', width=6, height=5)


## Changed v2 cluster names to match final decided names
v2_cluster_fullnames <- c(
  "0"  = "axon",
  "1"  = "ca1",
  "2"  = "dendr",
  "3"  = "ca3-4",
  "4"  = "luc",
  "5"  = "sgz-ml",
  "6"  = "inn",
  "7"  = "sub",
  "8"  = "or",
  "9"  = "rad",
  "10" = "vasc",
  "11" = "gcl",
  "12" = "cp",
  "13" = "sgz-pl",
  "14" = "ca1-rost",
  "15" = "ca2",
  "16" = "cr"
)

for (i in 0:15) { #11-7-2025 - Only care about clusters 0 to 15 for display purposes #I only did clusters 0 to 16 because 17 and 18 were returning an error: Error in `$<-.data.frame`(`*tmp*`, "Sig", value = "NS") : replacement has 1 row, data has 0
  z <- v2_degs %>% dplyr::filter(cell_type == i)
  
  
  p <- EnhancedVolcano(z,
                       lab = z$gene,
                       x = 'logFC',
                       y = 'adj.P.Val',
                       FCcutoff = 1,
                       pCutoff = 0.05,
                       labSize = 3.0,
                       drawConnectors = TRUE,
                       arrowheads = FALSE,
                       max.overlaps = 20,
                       ylab = "-log10 Adjusted P-value",
                       legendPosition = 'top',
                       title = paste0("DEGs for ", cluster_namer(i)),
                       subtitle = "MDD vs. CTRL",
                       ylim = c(0,8)) +
    theme_classic() +
    theme(legend.position = "none") +
    #xlim(c(-10,10)) #in order to include more genes
    xlim(c(-5, 5)) 
  #ylim(c(0, 3))
  print(p)
  
}

dev.off()


######Saving Excel Sheets
# Save as a tab separated excel. Make sure to save it to your own path.
df_sort <- v2_degs_sig %>% arrange(adj.P.Val, desc(logFC),.by_group = T)
df_list <- split(df_sort, as.factor(df_sort$cell_type))
write.xlsx(df_list,
           file = '/path/to/v2_degs_significant_no_sex_date.xlsx',
           asTable = TRUE)


