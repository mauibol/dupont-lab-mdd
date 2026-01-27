# Differential Expression with limma voom pipeline
# Maddy Peng

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
library(GWASTools)
library(rio)
library(ggbreak)
library(ggrepel)
library(car)
library(RColorBrewer)
library(EnhancedVolcano)

# Load data and metadata
multi <- LoadSeuratRds("/data/Share/pengmad/multiome_official/current_working_multimodal_bpcells.rds")

meta <- multi@meta.data
rownames(meta) <- meta$sample_barcode

# Create single cell experiment
counts <- multi[['RNA']]@layers$counts
counts <- as(counts, "dgCMatrix")
rownames(counts) <- rownames(multi)

sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = meta)
dim(counts(sce)) #36601 x 495037

# Remove genes not expressed in less than 0.1% of cells
sce <- sce[rowSums(counts(sce) > 0) >= 495, ]
dim(counts(sce)) #27510 x 495037

# Maker sure cluster column is a factor
sce$cell_annotation_leiden <- as.factor(sce$cell_annotation_leiden)
sce$sample <- as.factor(sce$sample)

# Pseudobulk by cluster and sample
pseudo <- aggregateAcrossCells(sce, 
                               DataFrame(
                                 cluster = sce$cell_annotation_leiden, #change to your col name
                                 sample = sce$sample
                               ))

# Filter out ncells less than 5
pseudo <- pseudo[,pseudo$ncells >= 5]

### Update the pseudobulk colData ###

# Winsorized RIN
colData(pseudo)$RIN_w <- dplyr::case_when(
  colData(pseudo)$RIN < 5.95 ~ 5.95,
  .default =colData(pseudo)$RIN
)

# Winsorized PMI
colData(pseudo)$PMI_w <- dplyr::case_when(
  colData(pseudo)$PMI > 32.5 ~ 32.5,
  .default = colData(pseudo)$PMI
)


colData(pseudo)$Sex <- dplyr::case_when(
  colData(pseudo)$Sex == '1' ~ 'Male',
  colData(pseudo)$Sex == '2' ~ 'Female',
  .default = colData(pseudo)$Sex
)

colData(pseudo)$Suicide <- dplyr::case_when(
  colData(pseudo)$Suicide == '1' ~ 'Suicide',
  colData(pseudo)$Suicide == '0' ~ 'NonSuicide',
  .default = colData(pseudo)$Suicide
)

# Make sure Diagnosis variable is factor in the correct order
colData(pseudo)$Diagnosis <- factor(colData(pseudo)$Diagnosis, levels=c('NONE', 'MDD'))
saveRDS(pseudo, paste0(METADATA_PATH, 'pseudo_leiden.rds'))

# Create a data frame with the sample IDs, cluster IDs and condition. Use whatever the column name is in your object
colData(pseudo)$experimenter <- as.factor(substr(colData(pseudo)$sample, 1, 2))
colData(pseudo)$batch <- as.factor(colData(pseudo)$batch)
colData(pseudo)$ID4 <- as.factor(colData(pseudo)$ID4)


# Loop through the clusters
cls <- unique(colData(pseudo)$cluster)
print(cls)

# Loop through all cell types : BATCH AS FIXED EFFECT
de_batch <- data.frame()

for (type in cls) {
  print(type)
  
  if (type %in% unique(de_batch$cell_type)) next
  
  pse <- pseudo[, pseudo$cluster == type]
  meta <- colData(pse)[, c("sample", "Diagnosis", "Age", "ncells", "ID4", "batch",
                           "Suicide", "RIN_w", "PMI_w", "Sex", "cluster", "experimenter")]
  
  if (length(unique(meta$Diagnosis)) == 1) {
    print(paste('Skipping', type, 'because only 1 diagnosis in this observation'))
    next
  }
  
  dge <- DGEList(counts = assays(pse)$counts, group = meta$Diagnosis)
  print(dim(dge))
  
  keep.exprs <- filterByExpr(dge, group = 'group')
  dge <- dge[keep.exprs,, keep.lib.sizes = FALSE]
  genes <- rownames(dge)
  dge <- calcNormFactors(dge)
  print(dim(dge))
  
  # Model design with conditional batch handling
  if (length(unique(meta$Sex)) == 1) {
    design <- model.matrix(~ Diagnosis + batch + scale(Age) + RIN_w + PMI_w + ncells, meta)
  } else {
    design <- model.matrix(~ Diagnosis + batch + Sex + scale(Age) + RIN_w + PMI_w + ncells, meta)
  }
  
  # If batch is not estimable, remove it
  nest <- nonEstimable(design)
  if (!is.null(nest)) {
    print("Not Estimable:")
    print(nest)
    if (any(grepl("^batch", nest))) {
      if (length(unique(meta$Sex)) == 1) {
        design <- model.matrix(~ Diagnosis + scale(Age) + RIN_w + PMI_w + ncells, meta)
      } else {
        design <- model.matrix(~ Diagnosis + Sex + scale(Age) + RIN_w + PMI_w + ncells, meta)
      }
    }
  }
  # weights estimated without the correlation
  vobj <- voomWithQualityWeights(dge, design, plot = FALSE)
  dupcor <- duplicateCorrelation(vobj, design, block = meta$ID4)
  print(dupcor$consensus.correlation)
  
  # this step uses the genome-wide average for random effect
  # Fit model
  print('fitting linear model...')
  fitDupCor <- tryCatch({
    fit <- lmFit(vobj, design, block = meta$ID4, correlation = dupcor$consensus)
    fit <- eBayes(fit)
    topTable_dupcor <- topTable(fit, number = length(genes), coef = 'DiagnosisMDD')
    topTable_dupcor$cell_type <- type
    topTable_dupcor$gene <- rownames(topTable_dupcor)
    de_batch <- rbind(de_batch, topTable_dupcor)
  }, error = function(e) {
    message("Model failed for ", type, ": ", e$message)
    return(NULL)
  })
}

# Save the DEGs
write.csv(de_batch,'/data/Share/pengmad/multiome_official/mdd_vs_ctrl_degs/unsupervised_degs_batch_as_fixed_062525.csv', row.names = T)

#### PLOTS ######
de_batch <- read.csv(paste0(METADATA_PATH, '/data/Share/pengmad/multiome_official/mdd_vs_ctrl_degs/unsupervised_degs_batch_as_fixed_062525.csv'), row.names = 1)
deg_sig_batch <- de_lasso_batch %>% dplyr::filter(adj.P.Val < 0.05)
deg_sig_batch$cell_type <- factor(deg_sig_batch$cell_type, levels = rev(unsupervised_cls_order))

# dotplot with significant genes
deg_sig_batch %>%
  ggplot(aes(x = logFC,y = cell_type, color  = adj.P.Val), alpha = 0.5)+
  geom_point()+
  scale_color_distiller(palette = "Blues")+
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  ggtitle(label = 'DEGS MDD vs. CTRL (Batch as Fixed Effect)')


# Save as a tab separated excel. Make sure to save it to your own path.
df_sort <- deg_sig_batch %>% arrange(adj.P.Val, desc(logFC),.by_group = T)
df_list <- split(df_sort, as.factor(df_sort$cell_type))
export(df_list, '/data/Share/pengmad/multiome_official/mdd_vs_ctrl_degs/unsupervised_degs_batch_as_fixed_062525.csv')


# Volcano Plots for all cell types
pdf(file = '/data/Share/pengmad/multiome_official/mdd_vs_ctrl_degs/degs_volcano_batch.pdf',onefile = TRUE, height = 8, width = 9)
for (c in cls) {
  ct <- de_batch %>% 
    dplyr::filter(cell_type==c)
  plot_sig <- ct %>% dplyr::filter(adj.P.Val < 0.05 & abs(logFC) >= 0.1) 
  genes <- ct$gene
  p <- ggplot(data = ct, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) + 
    geom_vline(xintercept = 0.25, linetype = 2) +
    geom_vline(xintercept = -0.25, linetype = 2) +
    geom_point(data = plot_sig, aes(x = logFC, y = -log10(adj.P.Val), color = 'black')) +
    geom_label_repel(data = plot_sig, aes(label = gene), size = 3) +
    ggtitle(paste("Cluster:", c)) +
    labs(subtitle = paste('Number of genes:', length(genes))) +
    theme(legend.position = 'none') +
    theme(text= element_text(size = 20))
  print(p)
}
dev.off()


########### FURTHER PLOTS #######################################
pdf(file = '/data/Share/pengmad/multiome_official/mdd_vs_ctrl_degs/degs_enhanced_volcano_batch.pdf',onefile = TRUE, height = 8, width = 9)
for (c in cls) {
  ct <- de_lasso_batch %>% dplyr::filter(cell_type==c)
  
  if (nrow(ct) == 0) {
    next
  }
  
  p <- EnhancedVolcano(ct,
                       lab = ct$gene,
                       x = 'logFC',
                       y = 'adj.P.Val',
                       FCcutoff = 1,
                       pCutoff = 0.05,
                       labSize = 4.0,
                       drawConnectors = TRUE,
                       arrowheads = FALSE,
                       max.overlaps = 15,
                       ylab = "-log10 Adjusted P-value",
                       legendPosition = 'top',
                       title = ct$cell_type,
                       subtitle = "MDD vs. CTRL"
  ) +
    theme_classic() +
    theme(legend.position = "none") +
    xlim(c(-4, 4)) +
    ylim(c(0, 5))
  
  print(p)
  
}
dev.off()


# Plot genes of interest
gene_of_interest <- c("ARC",'THY1', 'GABRD', 'ST8SIA3',"KCNIP4","EGR3","RPL32",
                      'TMEM47',"BHLHE40", "SCN2B"
)  # Replace with your gene of interest

# Extract logcounts for specific cell type PENK
penk <- pseudo[,pseudo$cluster =='InN.5']
penk <- logNormCounts(penk)
logcounts_matrix <- assays(penk)$logcounts

counts_matrix <- counts(penk)
logcpm_matrix <- cpm(counts_matrix, log = TRUE, prior.count = 1)


logcounts_gene <- logcpm_matrix[gene_of_interest, ]

df <- as.data.table(as.data.frame(t(logcounts_gene)))
df$Sample <- penk$sample
df$Cluster <- penk$cluster
df$Diagnosis <- penk$Diagnosis

df_long <- melt(df,
                id.vars = c("Sample", "Cluster", "Diagnosis"),
                variable.name = "Gene",
                value.name = "LogExpression")


ggplot(df_long, aes(x = Diagnosis, y = LogExpression, color = Diagnosis)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +  # individual dots
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +  # average lines
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +  # 1 row per gene
  theme_classic() +
  labs(
    x = NULL,
    y = "Log-Normalized Counts",
    title = "Log Normalized Counts for Top DEGs in InN.PENK.MME"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


p <- ggplot(df_long, aes(x = Diagnosis, y = LogExpression, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # hide outliers to reduce clutter
  geom_jitter(width = 0.2, alpha = 0.5, size = 1, color = "black") +  # optional
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Logcpm Normalized Counts",  # or "logCPM", etc.
    title = "Logcpm Normalized Counts for Top DEGs in GC.1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

pdf('/data/Share/pengmad/multiome_official/manuscript_figs/DEGs/InNPENK_DEGs_logcpm_boxplots_062625.pdf', width=8.5, height=6)
print(p)
dev.off()

# Extract logcounts for specific cell type GC1
gene_of_interest <- c('ENOX1', 'PILRB', 'SEMA6D', 'FLCN', 'ZFP37', 'CPNE8', 
                      'SERTM1', 'LITAF', 'CLPSL1', 'DSP')
library(sva)
library(edgeR)

gc <- pseudo[,pseudo$cluster =='GC.1']
gc <- logNormCounts(gc)

# Get logcounts and apply ComBat
log_expr <- logcounts(gc)
batch <- gc$ID4  # or gc$donor, gc$batch — depends on your metadata
log_expr_combat <- ComBat(dat = as.matrix(log_expr), batch = batch)

counts_matrix <- counts(gc)
logcpm_matrix <- cpm(counts_matrix, log = TRUE, prior.count = 1)

# Subset to genes of interest
logcounts_gene <- logcpm_matrix[gene_of_interest, ]

# Continue with reshaping and plotting
df <- as.data.table(as.data.frame(t(logcounts_gene)))
df$Sample <- gc$sample
df$Cluster <- gc$cluster
df$Diagnosis <- gc$Diagnosis

df_long <- melt(df,
                id.vars = c("Sample", "Cluster", "Diagnosis"),
                variable.name = "Gene",
                value.name = "LogExpression")

ggplot(df_long, aes(x = Diagnosis, y = LogExpression, color = Diagnosis)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Batch-corrected logExpr (ComBat)",
    title = "ComBat-Corrected Expression for Top DEGs in GC.1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Boxplot
p <- ggplot(df_long, aes(x = Diagnosis, y = LogExpression, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # hide outliers to reduce clutter
  geom_jitter(width = 0.2, alpha = 0.5, size = 1, color = "black") +  # optional
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Logcpm Normalized Counts",  # or "logCPM", etc.
    title = "Logcpm Normalized Counts for Top DEGs in GC.1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

pdf('/data/Share/pengmad/multiome_official/manuscript_figs/DEGs/GC1_DEGs_logcpm_boxplots_062625.pdf', width=8.5, height=6)
print(p)
dev.off()


# Extract logcounts for specific cell type CA1
gene_of_interest <- c('ANKRD6', 'FAM106C', 'ENTPD3', 'GPR26', 'SHE', "ADGRV1", 'MAP2K6',
                      'SEMA3D', 'C12orf42', 'STEAP2')

ex <- pseudo[,pseudo$cluster =='ExN.CA1-2.FIBCD1.FNDC1']
ex <- logNormCounts(ex)
logcounts_matrix <- assays(ex)$logcounts


counts_matrix <- counts(ex)
logcpm_matrix <- cpm(counts_matrix, log = TRUE, prior.count = 1)

logcounts_gene <- logcpm_matrix[gene_of_interest, ]

df <- as.data.table(as.data.frame(t(logcounts_gene)))
df$Sample <- ex$sample
df$Cluster <- ex$cluster
df$Diagnosis <- ex$Diagnosis

df_long <- melt(df,
                id.vars = c("Sample", "Cluster", "Diagnosis"),
                variable.name = "Gene",
                value.name = "LogExpression")


ggplot(df_long, aes(x = Diagnosis, y = LogExpression, color = Diagnosis)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +  # individual dots
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +  # average lines
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +  # 1 row per gene
  theme_classic() +
  labs(
    x = NULL,
    y = "Log-Normalized Counts",
    title = "Log Normalized Counts for Top DEGs in ExN.CA1-2"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Boxplot
p <- ggplot(df_long, aes(x = Diagnosis, y = LogExpression, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # hide outliers to reduce clutter
  geom_jitter(width = 0.2, alpha = 0.5, size = 1, color = "black") +  # optional
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Logcpm Normalized Counts",  # or "logCPM", etc.
    title = "Logcpm Normalized Counts for Top DEGs in ExN.CA1-2"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
pdf('/data/Share/pengmad/multiome_official/manuscript_figs/DEGs/ExNCA1-2_DEGs_logcpm_boxplots_062625.pdf', width=8.5, height=6)
print(p)
dev.off()
