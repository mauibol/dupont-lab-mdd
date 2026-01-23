#####################################################
# Visium v1 and v2 MDD vs Control Proportion Graphing
# By Victor Anosike
#####################################################

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(spatialLIBD)
  library(here)
  library(edgeR)
  library(scuttle)
  library(scater)
  library(scran)
  library(dplyr)
  library(gridExtra)
  library(ggforce)
  library(pheatmap)
  library(scater)
  library(scran)
  library(readxl)
  library(Matrix)
  library(rio)
  library(openxlsx)
  library(ggplot2)
  library(ggpubr)
  library(tibble)
  library(tidyverse)
})

################################# Visium v1 ###############################################
###########################################################################################

# ---- Load object ----
v1 <- readRDS("/insert/path/to/v1.rds")

# ---- Donor x cluster counts ----
counts <- as.data.frame.matrix(
  table(v1$donor, v1$seurat_clusters)
)

# ---- Normalize WITHIN donor (key change) ----
prop <- counts / rowSums(counts)

# ---- Add diagnosis metadata ----
meta <- v1@meta.data %>%
  distinct(donor, Psych1A) %>%
  arrange(donor)

prop$Diagnosis <- meta$Psych1A

# ---- Long format ----
df_long <- prop %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = matches("^\\d+$"),
    names_to = "Cluster",
    values_to = "Proportion"
  ) %>%
  filter(as.numeric(Cluster) <= 12)

df_long$Diagnosis <- factor(df_long$Diagnosis, levels = c("control", "MDD"))

# ---- Cluster name mapping ----
cluster_map <- c(
  `0` = "ecm",
  `1` = "luc-rad",
  `2` = "axon",
  `3` = "sgz-ml",
  `4` = "ca1-4.1",
  `5` = "ca1-4.2",
  `6` = "dendr",
  `7` = "sgz-pl",
  `8` = "gcl",
  `9` = "ca1-2",
  `10` = "inn",
  `11` = "vasc",
  `12` = "cp"
)

df_long$ClusterName <- cluster_map[df_long$Cluster]

# ---- Order clusters by descending MDD mean ----
cluster_order <- df_long %>%
  filter(Diagnosis == "MDD") %>%
  group_by(ClusterName) %>%
  summarise(mean_mdd = mean(Proportion, na.rm = TRUE)) %>%
  arrange(desc(mean_mdd)) %>%
  pull(ClusterName)

df_long$ClusterName <- factor(df_long$ClusterName, levels = cluster_order)

# ---- Wilcoxon tests ----
wilcox_res <- df_long %>%
  group_by(ClusterName) %>%
  summarise(
    p_value = wilcox.test(Proportion ~ Diagnosis)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", round(signif(p_value, 2), 2)))

wilcox_res$ClusterName <- factor(wilcox_res$ClusterName, levels = cluster_order)

# ---- Plot ----
ggplot(df_long, aes(x = ClusterName, y = Proportion, fill = Diagnosis)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.6,
    size = 1
  ) +
  geom_text(
    data = wilcox_res,
    aes(x = ClusterName, y = max(df_long$Proportion) * 1.05, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +

  scale_fill_manual(
    values = c("control" = "#f0697b", "MDD" = "#35c4d4"),
    labels = c("control" = "CTRL", "MDD" = "MDD")
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Cluster",
    y = "Proportion"
  ) +
  theme(
    axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right"
  )


################################################## Visium v2 ############################################
#########################################################################################################

# ---- Load object ----
v2 <- readRDS("insert/path/to/v2.rds")

# ---- Donor x cluster counts ----
counts <- as.data.frame.matrix(
  table(v2$donor, v2$seurat_clusters)
)

# ---- Normalize WITHIN donor (key change) ----
prop <- counts / rowSums(counts)

# ---- Add diagnosis metadata ----
meta <- v2@meta.data %>%
  distinct(donor, diagnosis) %>%
  arrange(donor)

prop$Diagnosis <- meta$diagnosis

# ---- Long format ----
df_2_long <- prop %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = matches("^\\d+$"),
    names_to = "Cluster",
    values_to = "Proportion"
  ) %>%
  filter(as.numeric(Cluster) <= 16)

df_2_long$Diagnosis <- factor(df_2_long$Diagnosis, levels = c("CTRL", "MDD"))

# ---- Cluster name mapping ----
cluster_2_map <- c(
  `0`  = "axon",
  `1`  = "ca1",
  `2`  = "dendr",
  `3`  = "ca3-4",
  `4`  = "luc",
  `5`  = "sgz-ml",
  `6`  = "inn",
  `7`  = "sub",
  `8`  = "or",
  `9`  = "rad",
  `10` = "vasc",
  `11` = "gcl",
  `12` = "cp",
  `13` = "sgz-pl",
  `14` = "ca1-rost",
  `15` = "ca2",
  `16` = "cr"
)

df_2_long$ClusterName <- cluster_2_map[df_2_long$Cluster]

# ---- Order clusters by descending CTRL mean ----
cluster_order <- df_2_long %>%
  filter(Diagnosis == "CTRL") %>%
  group_by(ClusterName) %>%
  summarise(mean_ctrl = mean(Proportion, na.rm = TRUE)) %>%
  arrange(desc(mean_ctrl)) %>%
  pull(ClusterName)

df_2_long$ClusterName <- factor(df_2_long$ClusterName, levels = cluster_order)

# ---- Wilcoxon tests (now statistically valid) ----
wilcox_res <- df_2_long %>%
  group_by(ClusterName) %>%
  summarise(
    p_value = wilcox.test(Proportion ~ Diagnosis)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", round(signif(p_value, 2), 2)))

wilcox_res$ClusterName <- factor(wilcox_res$ClusterName, levels = cluster_order)

# ---- Plot ----
ggplot(df_2_long, aes(x = ClusterName, y = Proportion, fill = Diagnosis)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.6,
    size = 1
  ) +
  geom_text(
    data = wilcox_res,
    aes(
      x = ClusterName,
      y = max(df_2_long$Proportion) * 1.05,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +

  scale_fill_manual(
    values = c("CTRL" = "#f0697b", "MDD" = "#35c4d4"),
    labels = c("CTRL" = "CTRL", "MDD" = "MDD")
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Cluster",
    y = "Proportion"
  ) +
  theme(
    axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right"
  )



