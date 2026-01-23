###########################################
# Visium v1 and v2 Integration, Clustering, and Plotting
# By Victor Anosike
###########################################

# Load Relevant Libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
  library(Polychrome)
  library(scCustomize)
  library(reticulate)
  library(gridExtra)
})

####################################### Visium v1 ###################################
#####################################################################################

# Set working directory
setwd('/data/Share/Victor_Anosike/src')

# Load seurat list of individual v2 samples
list <- readRDS("v1_sample_list.rds")

# Normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  DefaultAssay(x) <- "Spatial"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)

# Save anchors
saveRDS(anchors, "anchors_v1.rds")

# Genes are relevant to neurogenesis and will be added to list of genes to investigate
list_of_genes <- c('GAD1','GAD2','CALB1','CALB2','PCNA','MKI67','PAX6','TUBB3',
                   'VIM','HMGB2','ETNPPL','STMN1','STMN2','TOP2A','NES','DCX',
                   'PROX1','ASCL1','ASCL2','NEUROG1','NEUROG2','NEUROG3','NEUROD1',
                   'SOX1','SOX2','SOX6','SOX9','SOX13','EMX1','EMX2','ASCL3',
                   'ASCL4','DLX2','DLX5','DLX6','SOX11','OTX2','TBR2','OLIG1',
                   'OLIG2','FEZF2','TRB1','NEUROG3','EOMES','FEXF1','POU3F2',
                   'POU4F1','LHX2','LHX5','LHX6','MASH1','NPAS3','CDK5R1','MYT1',
                   'MYT1L','NRSN1','SERPINF1','FGF2','NTN1','NGF','BDNF','NTRK1',
                   'ROBO1','ROBO2','BMP4','SOX10','GFAP','FABP7','WNT3','WNT7A',
                   'NR2F1','NR2F2','NR2E1','RORB','ZEB2','POU6F1','POU6F2',
                   'POU6F3','POU6F4','POU6F5','POU6F6','EZH1','CAMK2A','FLT1',
                   'RBFOX3','SLC17A7','SYT1')

# List of genes to use for integration
genes_to_integrate <- unique(c(anchors@anchor.features,list_of_genes))

# Integrating list of individual v1 samples into one v1 object
v1 <- IntegrateData(anchorset = anchors, features.to.integrate = genes_to_integrate)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(v1) <- "integrated"

# Run the standard workflow for visualization and clustering
v1 <- ScaleData(v1, verbose = FALSE)
v1 <- RunPCA(v1, npcs = 30, verbose = FALSE)
v1 <- RunUMAP(v1, reduction = "pca", dims = 1:30)
v1 <- FindNeighbors(v1, reduction = "pca", dims = 1:30)
v1 <- FindClusters(v1, resolution = 0.34)

# save combined, integrated, and processed v1 object
saveRDS(v1, 'v1.rds')

# Load v1.rds as the R object "v1". This object has been integrated and processed
#with both Seurat and Precast. However, we will be plotting the portion of the object plotted from Seurat. 
#The resolution of used to FindClusters() was set to 0.34

v1 <- readRDS('v1.rds')

# UpdateSeuratObject in case of issues due to object structure
v1 <- UpdateSeuratObject(v1)

# Get length of clusters
length <- length(unique(v1@meta.data[["seurat_clusters"]]))

# Plot UMAP of integrated data
pdf("v1_plots/v1_umap_res34p.pdf")
DimPlot(v1, reduction = "umap", group.by = 'seurat_clusters', label = TRUE)
DimPlot(v1, reduction = "umap", group.by = 'sample', label = TRUE, repel = T)
DimPlot(v1, reduction = "umap", group.by = 'donor', label = TRUE, repel = T)
dev.off()

# Seurat Clusters Number Table
table(v1$seurat_clusters, v1$sample)
#     CM013 CM014 CM015 CM016 CM017 CM018 CM019 CM020 CM117 CM118 CM119 CM120 MM063 MM064 MM067 MM068 MM069 MM070 MM071 MM072 MM073
# 0    315   193   986   881   287   326   329   445   468   309   533   516   187   270   405   412   408   446   384   372   240
# 1    276   233   341   501   243   182   693   865   412   249   366   353   299   217    99   339   722   556   419   589   393
# 2    469   338   629   670   318   313   564   449   419   288   447   445   189   270   240   540   368   469   352   320   122
# 3    200   390   318   270   304   319   151   235   833   468   413   506   206   262   210   113   215   305   341   276   271
# 4    311   492   192   337   488   400    42   261   549   405   325   461   179   289   255   289   226   208   194   190   148
# 5    153   295   207   259   304   309   171   229   261   319   408   359   183   212   166   347   446   236   298   353   264
# 6    244   720   148   162   197   165   170   156   297   146   209   288   140   172   187   220   306   492   320   290   139
# 7    237   419   146   207   318   275   141   271   406   309   210   367   179   261   168    89   113   198   177   213   141
# 8    145   207   223   163   137   153    76   103   263   191   149   162   122   152   113    75   110   146   148   166   149
# 9     29   180   169   211    69    70   234   419   167   131    34    43   148    69    23   136    75   176   292   345    60
# 10    84   116   397   361    82    77   118   141   228    86   145   178    52    61    52    73    89   190   166   164    79
# 11    43    58    62    65    38    68    70    78   100    61    57    44    34    46    52    60    59   136   133    82    49
# 12    44   143    14    17    26    20    25    12    27    23     1     3     9    10     7    17     9    22    16    21     9
# 13     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     9    33     9    26    32
# 14     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     2    14     4    15
# 15     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     2     3    17     4     1     1

# Remove clusters with low counts
v1 <- subset(v1, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12))

# Rename clusters based off of differential enrichment and gene expression information (see v1_differential_enrichment.r)
v1$seurat_clusters <- dplyr::case_when(
  v1$seurat_clusters == 0 ~ 'ecm',
  v1$seurat_clusters == 1 ~ 'luc-rad',
  v1$seurat_clusters == 2 ~ 'axon',
  v1$seurat_clusters == 3 ~ 'sgz-ml',
  v1$seurat_clusters == 4 ~ 'ca1-4.1',
  v1$seurat_clusters == 5 ~ 'ca1-4.2',
  v1$seurat_clusters == 6 ~ 'dendr',
  v1$seurat_clusters == 7 ~ 'sgz-pl',
  v1$seurat_clusters == 8 ~ 'gcl',
  v1$seurat_clusters == 9 ~ 'ca1-2',
  v1$seurat_clusters == 10 ~ 'inn',
  v1$seurat_clusters == 11 ~ 'vasc',
  v1$seurat_clusters == 12 ~ 'cp',
  .default=v1$seurat_clusters
)

# Set order of clusters to appear when doing SpatialDimPlot with scCustomize
v1_order <- c('gcl','sgz-ml','sgz-pl','ca1-2','ca1-4.2','ca1-4.1','dendr','axon','cp','luc-rad','vasc','inn','ecm')
v1$seurat_clusters <- factor(v1$seurat_clusters, levels = v1_order)

# Graph spatial dim plots of v1 samples as pdfs (placed in v1_plots)
for (s in unique(v1$sample)) {
  pdf(paste0('v1_plots/',s, '_v1_spatialdimplot.pdf'))
  p <- SpatialDimPlot_scCustom(v1, group.by = 'seurat_clusters', combine = F, label=F, repel = T,alpha=0.8, images = s, stroke = NA) 
  
  print(p)
  dev.off()
}


####################################### Visium v2 ###################################
#####################################################################################

# Set working directory
setwd('/data/Share/Victor_Anosike/src')

# Load seurat list of individual v2 samples
list <- readRDS("v2_sample_list.rds")

# Normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  DefaultAssay(x) <- "Spatial"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)

# Save anchors
saveRDS(anchors, "anchors_v2.rds")

# Genes are relevant to neurogenesis and will be added to list of genes to investigate
list_of_genes <- c('GAD1','GAD2','CALB1','CALB2','PCNA','MKI67','PAX6','TUBB3',
                   'VIM','HMGB2','ETNPPL','STMN1','STMN2','TOP2A','NES','DCX',
                   'PROX1','ASCL1','ASCL2','NEUROG1','NEUROG2','NEUROG3','NEUROD1',
                   'SOX1','SOX2','SOX6','SOX9','SOX13','EMX1','EMX2','ASCL3',
                   'ASCL4','DLX2','DLX5','DLX6','SOX11','OTX2','TBR2','OLIG1',
                   'OLIG2','FEZF2','TRB1','NEUROG3','EOMES','FEXF1','POU3F2',
                   'POU4F1','LHX2','LHX5','LHX6','MASH1','NPAS3','CDK5R1','MYT1',
                   'MYT1L','NRSN1','SERPINF1','FGF2','NTN1','NGF','BDNF','NTRK1',
                   'ROBO1','ROBO2','BMP4','SOX10','GFAP','FABP7','WNT3','WNT7A',
                   'NR2F1','NR2F2','NR2E1','RORB','ZEB2','POU6F1','POU6F2',
                   'POU6F3','POU6F4','POU6F5','POU6F6','EZH1','CAMK2A','FLT1',
                   'RBFOX3','SLC17A7','SYT1')

# List of genes to use for integration
genes_to_integrate <- unique(c(anchors@anchor.features,list_of_genes))

# Integrating list of individual v2 samples into one v2 object
v2 <- IntegrateData(anchorset = anchors, features.to.integrate = genes_to_integrate)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(v2) <- "integrated"

# Run the standard workflow for visualization and clustering
v2 <- ScaleData(v2, verbose = FALSE)
v2 <- RunPCA(v2, npcs = 30, verbose = FALSE)
v2 <- RunUMAP(v2, reduction = "pca", dims = 1:30)
v2 <- FindNeighbors(v2, reduction = "pca", dims = 1:30)
v2 <- FindClusters(v2, resolution = 0.34)

# save combined, integrated, and processed v2 object
saveRDS(v2, 'v2.rds')

# Load v2.rds as the R object "v2"
v2 <- readRDS('v2.rds')

# UpdateSeuratObject in case of issues due to object structure
v2 <- UpdateSeuratObject(v2)

# Get length of clusters
length <- length(unique(v2@meta.data[["seurat_clusters"]]))

# Plot UMAP of integrated data
pdf("v2_plots/v2_umap_res34p.pdf")
DimPlot(v2, reduction = "umap", group.by = 'seurat_clusters', label = TRUE)
DimPlot(v2, reduction = "umap", group.by = 'sample', label = TRUE, repel = T)
DimPlot(v2, reduction = "umap", group.by = 'donor', label = TRUE, repel = T)
dev.off()

# Seurat Clusters Number Table
table(v2$seurat_clusters, v2$sample)
#     MM128 MM129 MM130 MM131 MM132 MM133 MM134 MM135 MM136 MM138 MM139 MM140 MM142 MM143 RM001 RM002 RM004 RM010
# 0   1201   503  2156  1624  1284  4361  1893  1580  1921  2970  2202   960  3059  1260   193   309  3274  2067
# 1    877  1063  1454  1664  2248   223  2685  2753   735   113  2189  1256   777  1570   448   259  1586  2637
# 2   1680   406   618   554  1182   875  1583  1826  1991   692   630  1642  1939   708   202   124  1105  1080
# 3   1582  1324   312  1060  1399  1149   383   742  1112   863   318  1704   676   195   343   471   759   881
# 4    977   689   183   462  1082   594   585   501  2317   660   381  1242  2734   452   645   216   396   444
# 5   1374  1089   185   583   755   639   598   796  2239   609   191  1299  1454   214   246   280   368   523
# 6    926   402   598   536   862   415   810   756   680   177   630   841   891   484   104   134   635   624
# 7    956   586  1059   402   676   199   870   508   198     8   869   798   286   858    38    13   415   929
# 8   1368   272   918   535   382    64   735   292    86     8  1140   676   106   784    57    26   503   405
# 9    216   538   669   344   632   108   702   934   342    45   599   435   340   441   170    86   704   354
# 10   296    74   205   196   217   785   198   156   334   620   205   332   333   161    22    31   140   159
# 11   484   260   134   305   388   260   168   248   437   192   145   424   334   115    93    87   165   223
# 12   470    71   111   450    61   698    41    57   205   116   212   375   144   149    42   431   268   219
# 13   278    90    66   219   149   233   233   341   369   195   118   338   261   134    69    82   183   101
# 14   301     1     0    78     1   232   536   363    17     2   213   120    34    61    25     6   566    42
# 15    95    88     4   126   202   331    74    55    21   132     1   119     3    11    71    56    96   235
# 16    67    21    56     7    87    28    31    20    15    20    40    56    17    19     2     0     7    26
# 17     0     0     0     0     0     0    75   108     0     0     0     0     0    69    11     0    29     0
# 18     2    94     1     0     2    38     0     3     2    11     0     0     1     2     3    41     3    74

# Remove clusters with low counts
v2 <- subset(v2, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))

# Rename clusters based off of differential enrichment and gene expression information (see v2_differential_enrichment.r)
v2$seurat_clusters <- dplyr::case_when(
  v2$seurat_clusters == 0 ~ 'axon',
  v2$seurat_clusters == 1 ~ 'ca1',
  v2$seurat_clusters == 2 ~ 'dendr',
  v2$seurat_clusters == 3 ~ 'ca3-4',
  v2$seurat_clusters == 4 ~ 'luc',
  v2$seurat_clusters == 5 ~ 'sgz-ml',
  v2$seurat_clusters == 6 ~ 'inn',
  v2$seurat_clusters == 7 ~ 'sub',
  v2$seurat_clusters == 8 ~ 'or',
  v2$seurat_clusters == 9 ~ 'rad',
  v2$seurat_clusters == 10 ~ 'vasc',
  v2$seurat_clusters == 11 ~ 'gcl',
  v2$seurat_clusters == 12 ~ 'cp',
  v2$seurat_clusters == 13 ~ 'sgz-pl',
  v2$seurat_clusters == 14 ~ 'ca1-rost',
  v2$seurat_clusters == 15 ~ 'ca2',
  v2$seurat_clusters == 16 ~ 'cr',
  .default=v2$seurat_clusters
)

# Set order of clusters to appear when doing SpatialDimPlot with scCustomize
v2_order <- c('gcl','sgz-ml','sgz-pl','ca1','ca2','ca3-4','dendr','axon','cp','luc','vasc','inn','sub','rad','ca1-rost',
              'cr','or')
v2$seurat_clusters <- factor(v2$seurat_clusters, levels = v2_order)

# Cosmetic Changes -- Due to Pink Color of Some Samples, Sample Color was adjusted in Order to Better Display Spots
pink_slides <- c("RM004", "MM140", "MM139", "MM138", "MM136", "MM135", "MM134", "MM133", "MM128", "RM010")
needs_bigger_dots <- c("MM134", "MM135", "MM139", "MM128", "MM140", "MM133", "MM138", "MM136", "MM142", "RM002", "RM004")

# Graph spatial dim plots of v1 samples as pdfs (placed in v1_plots)
for (s in unique(v2$sample)) {
  pdf(paste0('v2_plots/',s, '_v2_spatialdimplot.pdf'))
  if (s %in% pink_slides && s %in% needs_bigger_dots) {
    
    p <- SpatialDimPlot_scCustom(v2, group.by = 'seurat_clusters', combine = F, label=F, repel = T,alpha=1,images = s, image.alpha = 0.25, stroke = NA, pt.size.factor = 2.15) #from 1.95 to 2.15
    
  } else if (!(s %in% pink_slides) && s %in% needs_bigger_dots) {
    
    p <- SpatialDimPlot_scCustom(v2, group.by = 'seurat_clusters', combine = F, label=F, repel = T,alpha=1,images = s, image.alpha = 0.75, stroke = NA, pt.size.factor = 2.15) #from 1.95 to 2.15
    
  } else if (s %in% pink_slides && !(s %in% needs_bigger_dots)) {
    
    p <- SpatialDimPlot_scCustom(v2, group.by = 'seurat_clusters', combine = F, label=F, repel = T,alpha=1,images = s, image.alpha = 0.25, stroke = NA, pt.size.factor = 1.95) #from 1.85 to 1.95
    
    
  } else if ( !(s %in% pink_slides) && !(s %in% needs_bigger_dots)) {
    
    p <- SpatialDimPlot_scCustom(v2, group.by = 'seurat_clusters', combine = F, label=F, repel = T,alpha=1,images = s, image.alpha = 0.75, stroke = NA, pt.size.factor = 1.95) #from 1.85 to 1.95
    
    
  } else if (s == "MM129" || s == "MM132") { #These two samples are not pink AND their dots need to be even larger than the rest
    p <- SpatialDimPlot_scCustom(v2, group.by = 'seurat_clusters', combine = F, label=F, repel = T,alpha=1,images = s, image.alpha = 0.45, stroke = NA, pt.size.factor = 2.35) #from 2.2 to 2.35
  }
  
  print(p)
  dev.off()
}


