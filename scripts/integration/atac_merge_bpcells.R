# Maddy Peng
# Merge ATAC samples together
# Use BPCells to stream matrix from disk 

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Signac)
library(tidyverse)
library(scCustomize)
library(Matrix)
library(BPCells)
library(dplyr)
library(harmony)
library(ggplot2)
library(gridExtra)

# List sample ids
sample_ids <- c()