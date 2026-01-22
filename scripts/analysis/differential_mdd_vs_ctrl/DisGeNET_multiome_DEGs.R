#########################################
# MDD mulitome
# DisGeNET disease Gene Network analysis
# Anthony Ramnauth, June 02 2025
#########################################

# set working directory
library(devtools)
install_bitbucket("ibi_group/disgenet2r")
install_gitlab("medbio/disgenet2r")

suppressPackageStartupMessages({
  library(here)
  library(disgenet2r)
  library(readxl)
  library(knitr)
  library(tidyverse)
  library(openxlsx)
  library(ggraph)
  library(tidygraph)
  library(igraph)
  library(dplyr)
})

# here::here anchor
here::i_am("Disgenet.R")

api_key <- "d606cdd0-866c-411e-ab37-cf980ea9df3e"

Sys.setenv(DISGENET_API_KEY= api_key)

Proteomics <-read.csv("/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/proteomics_for_disgenet.csv")

GC.1 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "GC.1"))

GC.2 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "GC.2"))

GC.3 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "GC.3"))

GC.4 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "GC.4"))

ImGC <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "ImGC"))

VLMC <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "VLMC"))

Astro.1 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "Astro.1"))

Oligo.2 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "Oligo.2"))

OPC <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "OPC"))

Ependyma <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "Ependyma"))

Micro <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "Micro"))

InN.PENK <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "InN.PENK.MME"))

InN.SST.NPY <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "InN.SST.NPY"))

InN.SST.NDNF <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "InN.SST.NDNF"))

CA1 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "ExN.CA1-2.FIBCD1.FNDC1"))

CA_prox <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "ExN.proximalSUB.MUC19.GLP2R"))

CA_distal <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "ExN.distalSUB.MSLNL.ANKRD34B"))

CA_SUB_IQCF3 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "ExN.SUB.IQCF3.KRT17"))

CA_SUB_SMYD1 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/multiome_DEG_MDD_vs_CTRL_significant_by_celltype_June_2025.xlsx",
  sheet = "ExN.SUB.SMYD1.THEMIS"))

GC.1_genes <- GC.1$gene
InN.PENK_genes <- InN.PENK$gene
CA1_genes <- CA1$gene
CA_prox_genes <- CA_prox$gene
Astro.1_genes <- Astro.1$gene
Oligo.2_genes <- Oligo.2$gene
OPC_genes <- OPC$gene
Ependyma_genes <- Ependyma$gene
Micro_genes <- Micro$gene
ImGC_genes <- ImGC$gene
VLMC_genes <- VLMC$gene
InN.SST.NPY_genes <- InN.SST.NPY$gene
InN.SST.NDNF_genes <- InN.SST.NDNF$gene
CA_SUB_IQCF3_genes <- CA_SUB_IQCF3$gene
CA_SUB_SMYD1_genes <- CA_SUB_SMYD1$gene
CA_distal_genes <- CA_distal$gene

# Remove "_HUMAN" for cleaner labels
Proteomics <- Proteomics %>%
  mutate(Label = gsub("_HUMAN", "", Uniprot.ID))

Proteomics_genes <- Proteomics$Label


Proteomics_genes1 <- Proteomics_genes[1:100]
Proteomics_genes2 <- Proteomics_genes[101:200]
Proteomics_genes3 <- Proteomics_genes[201:300]
Proteomics_genes4 <- Proteomics_genes[300:308]

Proteomics_result3 <- gene2disease(
  gene     = Proteomics_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)

Proteomics_tab3 <- Proteomics_result3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

Proteomics_all_tab <- unique(c(Proteomics_tab1$gene_symbol, Proteomics_tab2$gene_symbol, Proteomics_tab3$gene_symbol, Proteomics_tab4$gene_symbol))

Proteomics_results <- gene2disease(
  gene     = Proteomics_all_tab,
  database = "PSYGENET",
  verbose  = TRUE
)


Proteomics_tab <- Proteomics_results@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

Proteomics_df <- Proteomics_results@qresult

# Save to Excel
write.xlsx(Proteomics_df, file = "PsyGeNET_results_Proteomics.xlsx")


df <- InN_results@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- Proteomics_df %>%
  filter(score >= 0.5)  # adjust threshold if needed

edges <- df_filtered %>%
  dplyr::select(gene_symbol, disease_name) %>%
  distinct() %>%
  rename(from = gene_symbol, to = disease_name)

# 2. Build igraph object
net <- graph_from_data_frame(edges, directed = FALSE)

# 3. Assign node type (Gene vs Disease)
V(net)$type <- ifelse(V(net)$name %in% Proteomics_results@qresult$gene_symbol,
                      "Gene", "Disease")

# 4. Convert to tidygraph for ggraph plotting
g <- as_tbl_graph(net)

# 5. Plot with larger text
ggraph(g, layout = "fr") + 
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 10) +
  geom_node_text(aes(label = name, color = type),
                 repel = TRUE,
                 fontface = "bold",
                 size = 10) +   # 🔥 increase text size here
  scale_color_manual(values = c(Gene = "blue", Disease = "hotpink")) +
  theme_void()


head(g)

ggraph(g, layout = 'fr') +
  geom_edge_link(color = "black", width = 0.8)+  # darker lines
  geom_node_point(aes(color = type), size = 9) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 9) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")

pdf(here::here("plots", "Proteomics"), width = 10, height = 10)

p1 <- disgenet2r::plot(Proteomics_results, type  ="Heatmap", limit  = 100, nchars = 60, verbose = T)
p1 + scale_fill_gradient(low = "white", high = "blue")

dev.off()


p1 +
  scale_fill_gradient(low = "white", high = "blue") +
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  )


proteomics_df <- Proteomics_results@qresult

write.csv(proteomics_df, "Proteomics_results_disgenet.csv")
write.xlsx(proteomics_df, file = "PsyGeNET_results_Proteomics.xlsx")


# First run for GN.1
GC.1_genes <- GC.1$gene

# Limited to 100 genes are one time, split GC.1_genes into 4 with 100, 100, 100, and 76 entries
GC.1_genes1 <- GC.1_genes[1:100]
GC.1_genes2 <- GC.1_genes[101:200]
GC.1_genes3 <- GC.1_genes[201:300]
GC.1_genes4 <- GC.1_genes[301:376]

GC1_results1 <- gene2disease(
  gene     = GC.1_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

GC1_results2 <- gene2disease(
  gene     = GC.1_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

GC1_results3 <- gene2disease(
  gene     = GC.1_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)

GC1_results4 <- gene2disease(
  gene     = GC.1_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

GC1_tab1 <- GC1_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)
GC1_tab2 <- GC1_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)
GC1_tab3 <- GC1_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)
GC1_tab4 <- GC1_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

GC1_all_tab <- unique(c(GC1_tab1$gene_symbol, GC1_tab2$gene_symbol, GC1_tab3$gene_symbol, GC1_tab4$gene_symbol))

# Rerun with just this list
GC1_results <- gene2disease(
  gene     = GC1_all_tab,
  database = "PSYGENET",
  verbose  = TRUE
)

GC1_tab <- GC1_results@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)



#enlarge the text for the network plot 

library(igraph)
library(dplyr)

# Use your DisGeNET results table
df <- GC1_results@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.1)  # adjust threshold if needed

# Create edge list: genes -> diseases
#edges <- df_filtered %>%
  #select(from = gene_symbol, to = disease_name)


# 1. Extract edges (gene–disease associations)
edges <- df_filtered %>%
  dplyr::select(gene_symbol, disease_name) %>%
  distinct() %>%
  rename(from = gene_symbol, to = disease_name)

# 2. Build igraph object
net <- graph_from_data_frame(edges, directed = FALSE)

# 3. Assign node type (Gene vs Disease)
V(net)$type <- ifelse(V(net)$name %in% GC1_results@qresult$gene_symbol,
                      "Gene", "Disease")

# 4. Convert to tidygraph for ggraph plotting
g <- as_tbl_graph(net)

# 5. Plot with larger text
ggraph(g, layout = "fr") + 
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name, color = type),
                 repel = TRUE,
                 fontface = "bold",
                 size = 7) +   # 🔥 increase text size here
  scale_color_manual(values = c(Gene = "blue", Disease = "hotpink")) +
  theme_void()



ggraph(g, layout = 'fr') +
  geom_edge_link(color = "black", width = 0.8)+  # darker lines
  geom_node_point(aes(color = type), size = 9) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 9) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")

# Set PDF output
pdf("network_plot_GC.1.pdf", width = 10, height = 10)  # adjust width/height as needed

# Plot the network
ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "darkblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired

# Close the PDF device
dev.off()


pdf(here::here("plots", "multiome", "mixedeffectsDEGs", "PsyGeNET", "PsyGeNET_heatmap_GC.1.pdf"), width = 10, height = 10)

p1 <- disgenet2r::plot(GC1_results, type  ="Heatmap", limit  = 100, nchars = 60, verbose = T)
p1 + scale_fill_gradient(low = "white", high = "blue")

dev.off()



pdf("GC1_results_Heatmap.pdf", width = 12, height = 10)  # adjust size as needed

# Generate the heatmap
p1 <- disgenet2r::plot(
  GC1_results,
  type  = "Heatmap",
  limit = 100,
  nchars = 60,
  verbose = TRUE
)

# Customize color gradient and text size
p1 +
  scale_fill_gradient(low = "white", high = "blue") +
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  )

# Close the PDF device
dev.off()





# Now run for InN.PENK

InN.PENK_genes <- InN.PENK$gene

# Limited to 100 genes are one time, split InN.PENK_genes into 6 with 100, 100, 100, 100, 100, and 65 entries
InN.PENK_genes1 <- InN.PENK_genes[1:100]
InN.PENK_genes2 <- InN.PENK_genes[101:200]
InN.PENK_genes3 <- InN.PENK_genes[201:300]
InN.PENK_genes4 <- InN.PENK_genes[301:400]
InN.PENK_genes5 <- InN.PENK_genes[401:500]
InN.PENK_genes6 <- InN.PENK_genes[501:565]

InN.PENK_results1 <- gene2disease(
  gene     = InN.PENK_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

InN.PENK_results2 <- gene2disease(
  gene     = InN.PENK_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

InN.PENK_results3 <- gene2disease(
  gene     = InN.PENK_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)

InN.PENK_results4 <- gene2disease(
  gene     = InN.PENK_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

InN.PENK_results5 <- gene2disease(
  gene     = InN.PENK_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

InN.PENK_results6 <- gene2disease(
  gene     = InN.PENK_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)

InN.PENK_tab1 <- InN.PENK_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

InN.PENK_tab2 <- InN.PENK_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

InN.PENK_tab3 <- InN.PENK_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

InN.PENK_tab4 <- InN.PENK_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

InN.PENK_tab5 <- InN.PENK_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

InN.PENK_tab6 <- InN.PENK_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

InN.PENK_all_tab <- unique(c(InN.PENK_tab1$gene_symbol, InN.PENK_tab2$gene_symbol, InN.PENK_tab3$gene_symbol, 
                              InN.PENK_tab4$gene_symbol, InN.PENK_tab5$gene_symbol, InN.PENK_tab6$gene_symbol))

# Rerun with just this list
InN.PENK_results <- gene2disease(
  gene     = InN.PENK_all_tab,
  database = "PSYGENET",
  verbose  = TRUE
)


knitr::kable(InN.PENK_tab[1:79,], caption = "GDAs for the list of genes for InN.PENK")

pdf(here::here("plots", "multiome", "mixedeffectsDEGs", "PsyGeNET", "PsyGeNET_network_InN.PENK.MME.pdf"), width = 10, height = 10)

disgenet2r::plot(InN.PENK_results, prop = 10, type = "Network",  verbose = T)

dev.off()


#enlarge the text for the network plot 

library(igraph)
library(dplyr)
library(ggraph)
library(tidygraph)

# Use your DisGeNET results table
df <- InN.PENK_results@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.65)  # adjust threshold if needed

# Create edge list: genes -> diseases
edges <- df_filtered %>%
  select(from = gene_symbol, to = disease_name)

# 1. Extract edges (gene–disease associations)
edges <- df_filtered %>%
  dplyr::select(gene_symbol, disease_name) %>%
  distinct() %>%
  rename(from = gene_symbol, to = disease_name)

# 2. Build igraph object
net <- graph_from_data_frame(edges, directed = FALSE)

# 3. Assign node type (Gene vs Disease)
V(net)$type <- ifelse(V(net)$name %in% InN.PENK_results@qresult$gene_symbol,
                      "Gene", "Disease")

# 4. Convert to tidygraph for ggraph plotting
g <- as_tbl_graph(net)

# 5. Plot with larger text
ggraph(g, layout = 'fr') +
  geom_edge_link(color = "black", width = 0.8)+  # darker lines
  geom_node_point(aes(color = type), size = 9) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 9) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")



ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "darkblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired

# Set PDF output
pdf("network_plot_InN.PENK.MME.pdf", width = 10, height = 10)  # adjust width/height as needed

# Plot the network
ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "darkblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired

# Close the PDF device
dev.off()



pdf(here::here("plots", "multiome", "mixedeffectsDEGs", "PsyGeNET", "PsyGeNET_heatmap_InN.PENK.MME.pdf"), width = 14, height = 10)

p1 <- disgenet2r::plot(InN.PENK_results, type  ="Heatmap", limit  = 100, nchars = 60, verbose = T)
p1 + scale_fill_gradient(low = "white", high = "blue")

dev.off()


pdf("InN.PENK.MME_results_Heatmap.pdf", width = 18, height = 10)  # adjust size as needed

# Generate the heatmap
p1 <- disgenet2r::plot(
  InN.PENK_results,
  type  = "Heatmap",
  limit = 100,
  nchars = 60,
  verbose = TRUE
)

# Customize color gradient and text size
p1 +
  scale_fill_gradient(low = "white", high = "blue") +
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  )

# Close the PDF device
dev.off()




# Now run for CA1
CA1_genes <- CA1$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
CA1_genes1 <- CA1_genes[1:100]
CA1_genes2 <- CA1_genes[101:145]

CA1_results1 <- gene2disease(
  gene     = CA1_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

CA1_results2 <- gene2disease(
  gene     = CA1_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

CA1_tab1 <- CA1_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)
CA1_tab2 <- CA1_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

CA1_all_tab <- unique(c(CA1_tab1$gene_symbol, CA1_tab2$gene_symbol))

# Rerun with just this list
CA1_results <- gene2disease(
  gene     = CA1_all_tab,
  database = "PSYGENET",
  verbose  = TRUE
)

CA1_tab <- CA1_results@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

knitr::kable(CA1_tab[1:15,], caption = "GDAs for the list of genes for CA1")

pdf(here::here("plots", "multiome", "mixedeffectsDEGs", "PsyGeNET", "PsyGeNET_network_ExN.CA1-2.FIBCD1.FNDC1.pdf"), width = 10, height = 10)

disgenet2r::plot(CA1_results, prop = 10, type = "Network",  verbose = T)

dev.off()


#enlarge the text for the network plot 

library(igraph)
library(dplyr)
library(ggraph)
library(tidygraph)

# Use your DisGeNET results table
df <- CA1_results@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.1)  # adjust threshold if needed

# Create edge list: genes -> diseases
edges <- df_filtered %>%
  select(from = gene_symbol, to = disease_name)


# 1. Extract edges (gene–disease associations)
edges <- df_filtered %>%
  dplyr::select(gene_symbol, disease_name) %>%
  distinct() %>%
  rename(from = gene_symbol, to = disease_name)

# 2. Build igraph object
net <- graph_from_data_frame(edges, directed = FALSE)

# 3. Assign node type (Gene vs Disease)
V(net)$type <- ifelse(V(net)$name %in% CA1_results@qresult$gene_symbol,
                      "Gene", "Disease")

# 4. Convert to tidygraph for ggraph plotting
g <- as_tbl_graph(net)

# 5. Plot with larger text
ggraph(g, layout = 'fr') +
  geom_edge_link(color = "black", width = 0.8)+  # darker lines
  geom_node_point(aes(color = type), size = 9) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 9) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")




ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "darkblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired

# Set PDF output
pdf("network_plot_CA1.pdf", width = 10, height = 10)  # adjust width/height as needed

# Plot the network
ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "darkblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired

# Close the PDF device
dev.off()


pdf(here::here("plots", "multiome", "mixedeffectsDEGs", "PsyGeNET", "PsyGeNET_heatmap_ExN.CA1-2.FIBCD1.FNDC1.pdf"), width = 10, height = 10)

p1 <- disgenet2r::plot(CA1_results, type  ="Heatmap", limit  = 100, nchars = 60, verbose = T)
p1 + scale_fill_gradient(low = "white", high = "blue")

dev.off()


pdf("CA1_results_Heatmap.pdf", width = 10, height = 10)  # adjust size as needed

# Generate the heatmap
p1 <- disgenet2r::plot(
  CA1_results,
  type  = "Heatmap",
  limit = 100,
  nchars = 60,
  verbose = TRUE
)

# Customize color gradient and text size
p1 +
  scale_fill_gradient(low = "white", high = "blue") +
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold")
  )

# Close the PDF device
dev.off()



#Save excels for all 
library(openxlsx)   # or writexl, or readr for CSV

# Extract the PsyGeNET results table
GC.1_df <- GC1_results@qresult

# Save to Excel
write.xlsx(GC.1_df, file = "PsyGeNET_results_GC.1.xlsx")

# Or save to CSV instead
write.csv(GC.1_df, file = "PsyGeNET_results_GC.1.csv", row.names = FALSE)


# Extract the PsyGeNET results table
ExN.1_df <- CA1_results@qresult

# Save to Excel
write.xlsx(ExN.1_df, file = "PsyGeNET_results_ExN.1.xlsx")

# Extract the PsyGeNET results table
InN.PENK_df <- InN.PENK_results@qresult

# Save to Excel
write.xlsx(InN.PENK_df, file = "PsyGeNET_results_InN.PENK.MME.xlsx")


rm(merged_df)

# Get all objects in environment that end with "_df"
df_names <- ls(pattern = "_df$")

# Combine them into one dataframe, adding the source name as a column
merged_df <- do.call(rbind, lapply(df_names, function(df_name) {
  df <- get(df_name)
  df$source_df <- df_name   # add a new column with dataframe name
  df
}))


merged_df_fixed <- merged_df %>%
  mutate(across(where(is.list), ~sapply(., function(x) paste(unlist(x), collapse = ";"))))

write.csv(merged_df_fixed, "merged_dfs_Disgenet.csv", row.names = FALSE)



library(dplyr)
library(httr)
library(jsonlite)

get_uniprot_function <- function(uniprot_id) {
  # remove any suffixes
  uid <- sub("_HUMAN.*$", "", uniprot_id)
  # fetch from UniProt API
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uid, ".json")
  resp <- GET(url)
  if (status_code(resp) != 200) return(NA)
  js <- content(resp, "parsed", simplifyVector = TRUE)
  # find the function/comment field
  if (!is.null(js$comments)) {
    # comments list, find type=FUNCTION
    funcoms <- js$comments[js$comments$type == "FUNCTION"]
    if (length(funcoms) >= 1) {
      return(funcoms[[1]]$text[[1]])
    }
  }
  return(NA)
}

# Example: add a function column by Uniprot ID
merged_df2 <- merged_df %>%
  mutate(
    clean_uniprot = sub("_HUMAN.*$", "", uniprotids),  # assuming your column "uniprotids"
    gene_function_uniprot = sapply(clean_uniprot, get_uniprot_function)
  )

# Inspect
head(merged_df2[, c("gene_symbol","uniprotids","gene_function_uniprot")])




#Neural trajectory 

NSC_INP_NB <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/trajectory_DEGS_mdd_vs_ctrl_10.10.25.csv.xlsx",
  sheet = "NSC_INP_NB"))

# First run 
NSC_INP_NB_genes <- NSC_INP_NB$gene
# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
NSC_INP_NB_genes1 <- NSC_INP_NB_genes[1:100]
NSC_INP_NB_genes2 <- NSC_INP_NB_genes[101:200]
NSC_INP_NB_genes3 <- NSC_INP_NB_genes[201:300]
NSC_INP_NB_genes4 <- NSC_INP_NB_genes[301:400]
NSC_INP_NB_genes5 <- NSC_INP_NB_genes[401:500]
NSC_INP_NB_genes6 <- NSC_INP_NB_genes[501:600]
NSC_INP_NB_genes7 <- NSC_INP_NB_genes[601:700]
NSC_INP_NB_genes8 <- NSC_INP_NB_genes[701:800]
NSC_INP_NB_genes9 <- NSC_INP_NB_genes[801:900]
NSC_INP_NB_genes10 <- NSC_INP_NB_genes[901:1000]
NSC_INP_NB_genes11 <- NSC_INP_NB_genes[1001:1100]
NSC_INP_NB_genes12 <- NSC_INP_NB_genes[1101:1134]

NSC_INP_NB_results1 <- gene2disease(
  gene     = NSC_INP_NB_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results2 <- gene2disease(
  gene     = NSC_INP_NB_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results3 <- gene2disease(
  gene     = NSC_INP_NB_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results3 <- gene2disease(
  gene     = NSC_INP_NB_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results3 <- gene2disease(
  gene     = NSC_INP_NB_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results4 <- gene2disease(
  gene     = NSC_INP_NB_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results5 <- gene2disease(
  gene     = NSC_INP_NB_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results6 <- gene2disease(
  gene     = NSC_INP_NB_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results7 <- gene2disease(
  gene     = NSC_INP_NB_genes7,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results8 <- gene2disease(
  gene     = NSC_INP_NB_genes8,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results9 <- gene2disease(
  gene     = NSC_INP_NB_genes9,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results10 <- gene2disease(
  gene     = NSC_INP_NB_genes10,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results11 <- gene2disease(
  gene     = NSC_INP_NB_genes11,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results12 <- gene2disease(
  gene     = NSC_INP_NB_genes12,
  database = "PSYGENET",
  verbose  = TRUE
)



NSC_INP_NB_tab1 <- NSC_INP_NB_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab2 <- NSC_INP_NB_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab3 <- NSC_INP_NB_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab4 <- NSC_INP_NB_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab5 <- NSC_INP_NB_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab6 <- NSC_INP_NB_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab7 <- NSC_INP_NB_results7@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab8 <- NSC_INP_NB_results8@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab9 <- NSC_INP_NB_results9@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab10 <- NSC_INP_NB_results10@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab11 <- NSC_INP_NB_results11@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_tab12 <- NSC_INP_NB_results12@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

NSC_INP_NB_all_tab <- unique(c(NSC_INP_NB_tab1$gene_symbol, NSC_INP_NB_tab2$gene_symbol,NSC_INP_NB_tab3$gene_symbol,NSC_INP_NB_tab4$gene_symbol,NSC_INP_NB_tab5$gene_symbol,NSC_INP_NB_tab6$gene_symbol,NSC_INP_NB_tab7$gene_symbol,NSC_INP_NB_tab8$gene_symbol,NSC_INP_NB_tab9$gene_symbol,NSC_INP_NB_tab10$gene_symbol,NSC_INP_NB_tab11$gene_symbol,NSC_INP_NB_tab12$gene_symbol))

NSC_INP_NB_all_tab1 <- NSC_INP_NB_all_tab[1:100]
NSC_INP_NB_all_tab2 <- NSC_INP_NB_all_tab[101:127]

# Rerun with just this list
NSC_INP_NB_results_all_1 <- gene2disease(
  gene     = NSC_INP_NB_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)

NSC_INP_NB_results_all_2 <- gene2disease(
  gene     = NSC_INP_NB_all_tab2,
  database = "PSYGENET",
  verbose  = TRUE
)


# Use your DisGeNET results table
df1 <- NSC_INP_NB_results_all_1@qresult
df2 <- NSC_INP_NB_results_all_2@qresult

library(dplyr)

merged_df <- bind_rows(df1, df2)

# Optional: filter by score or number of supporting PMIDs
df_filtered <- merged_df %>%
  filter(score >= 0.5)  # adjust threshold if needed

# Create edge list: genes -> diseases
#edges <- df_filtered %>%
#select(from = gene_symbol, to = disease_name)


# 1. Extract edges (gene–disease associations)
edges <- df_filtered %>%
  dplyr::select(gene_symbol, disease_name) %>%
  distinct() %>%
  rename(from = gene_symbol, to = disease_name)

# 2. Build igraph object
net <- graph_from_data_frame(edges, directed = FALSE)

# 3. Assign node type (Gene vs Disease)
V(net)$type <- ifelse(
  V(net)$name %in% NSC_INP_NB_results_all_1@qresult$gene_symbol & 
    V(net)$name %in% NSC_INP_NB_results_all_2@qresult$gene_symbol,
  "Gene", 
  "Disease"
)

V(net)$type <- ifelse(
  V(net)$name %in% df_filtered$gene_symbol,
  "Gene", 
  "Disease"
)

# 4. Convert to tidygraph for ggraph plotting
g <- as_tbl_graph(net)

# 5. Plot with larger text
ggraph(g, layout = "fr") + 
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name, color = type),
                 repel = TRUE,
                 fontface = "bold",
                 size = 7) +   # 🔥 increase text size here
  scale_color_manual(values = c(Gene = "blue", Disease = "hotpink")) +
  theme_void()

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired




V1_0 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/v1_degs_significant_no_sex_10_23_2025.xlsx",
  sheet = "0"))

ATAC <- read_csv("significant_DARs_unsupervised_batch_as_fixed_peak_annotations_10.28.25.csv")

# First run 
ATAC_genes <- ATAC$SYMBOL

ATAC_results <- gene2disease(
  gene     = ATAC_genes,
  database = "PSYGENET",
  verbose  = TRUE
)





ImGC_df <- ImGC_df_NT
NSC_INP_NB_df <- df_filtered

# Save to Excel
write.xlsx(merged_df, file = "PsyGeNET_results_NSC_INP_NB_10.15.25.xlsx")






rm(merged_df)

# Get all objects in environment that end with "_df"
df_names <- ls(pattern = "_df$")

# Combine them into one dataframe, adding the source name as a column
merged_df <- do.call(rbind, lapply(df_names, function(df_name) {
  df <- get(df_name)
  df$source_df <- df_name   # add a new column with dataframe name
  df
}))


merged_df_fixed <- merged_df %>%
  mutate(across(where(is.list), ~sapply(., function(x) paste(unlist(x), collapse = ";"))))

write.csv(merged_df_fixed, "merged_dfs_Disgenet_with_neural_trajectory_10.15.25.csv", row.names = FALSE)

head(merged_df_fixed)

head(Proteomics_df)

# Identify overlapping genes
common_genes <- intersect(merged_df_fixed$gene_symbol, Proteomics_df$gene_symbol)

# Subset each dataframe to only those genes
merged_df_sub <- merged_df_fixed[merged_df_fixed$gene_symbol %in% common_genes, ]
proteomics_df_sub <- Proteomics_df[Proteomics_df$gene_symbol %in% common_genes, ]

# Merge the two dataframes by gene_symbol
merged_combined <- merge(merged_df_sub, proteomics_df_sub, 
                         by = "gene_symbol", 
                         suffixes = c("_psygenet", "_proteomics"))

# View the merged dataframe
head(merged_combined)

sapply(merged_combined, class)

library(dplyr)

merged_combined_fixed <- merged_combined %>%
  mutate(across(where(is.list), ~sapply(., function(x) paste(unique(unlist(x)), collapse = "; "))))

write.csv(merged_combined_fixed, "overlapping_proteomics_and_multiome_psygenet.csv", row.names = FALSE)


common_genes_2 <- intersect(ImGC_genes, Proteomics_genes) #CPNE8 only overlapping 
common_genes_3 <- intersect(NSC_INP_NB_genes, Proteomics_genes) #no overlap


head(merged_df_fixed)

Gene_ala_cart_proteomics_summary <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/Genealacart_proteomics_psygenet.xlsx",
  sheet = "Summaries"))


library(dplyr)

# Merge the data frames first
all_genecard_name <- full_join(Gene_ala_cart_proteomics_gene, 
                               Gene_ala_cart_gene, 
                               by = "Symbol")

# Identify columns with .x and .y suffixes
common_cols <- intersect(
  names(all_genecard_name)[grepl("\\.x$", names(all_genecard_name))] %>% sub("\\.x$", "", .),
  names(all_genecard_name)[grepl("\\.y$", names(all_genecard_name))] %>% sub("\\.y$", "", .)
)

# Loop through each duplicated column and merge them
for (col in common_cols) {
  all_genecard_name[[col]] <- coalesce(all_genecard_name[[paste0(col, ".x")]],
                                       all_genecard_name[[paste0(col, ".y")]])
  all_genecard_name[[paste0(col, ".x")]] <- NULL
  all_genecard_name[[paste0(col, ".y")]] <- NULL
}

# Resulting data frame now has merged columns
head(all_genecard_name)

# Merge the data frames first
all_genecard_alias <- full_join(Gene_ala_cart_proteomics_alias, 
                                Gene_ala_cart_aliases, 
                                by = "Symbol")

# Identify columns with .x and .y suffixes
common_cols <- intersect(
  names(all_genecard_alias)[grepl("\\.x$", names(all_genecard_alias))] %>% sub("\\.x$", "", .),
  names(all_genecard_alias)[grepl("\\.y$", names(all_genecard_alias))] %>% sub("\\.y$", "", .)
)

# Loop through each duplicated column and merge them
for (col in common_cols) {
  all_genecard_alias[[col]] <- coalesce(all_genecard_alias[[paste0(col, ".x")]],
                                        all_genecard_alias[[paste0(col, ".y")]])
  all_genecard_alias[[paste0(col, ".x")]] <- NULL
  all_genecard_alias[[paste0(col, ".y")]] <- NULL
}


head(all_genecard_alias)


# Merge the data frames first
all_genecard_summary <- full_join(Gene_ala_cart_proteomics_summary, 
                                Gene_ala_cart_summary, 
                                by = "Symbol")

# Identify columns with .x and .y suffixes
common_cols <- intersect(
  names(all_genecard_summary)[grepl("\\.x$", names(all_genecard_summary))] %>% sub("\\.x$", "", .),
  names(all_genecard_summary)[grepl("\\.y$", names(all_genecard_summary))] %>% sub("\\.y$", "", .)
)

# Loop through each duplicated column and merge them
for (col in common_cols) {
  all_genecard_summary[[col]] <- coalesce(all_genecard_summary[[paste0(col, ".x")]],
                                          all_genecard_summary[[paste0(col, ".y")]])
  all_genecard_summary[[paste0(col, ".x")]] <- NULL
  all_genecard_summary[[paste0(col, ".y")]] <- NULL
}


head(all_genecard_summary)
head(all_genecard_alias)
head(all_genecard_name)

library(dplyr)

# Keep only the first row per Symbol in alias and name tables
alias_unique <- all_genecard_alias %>% distinct(Symbol, .keep_all = TRUE)
name_unique  <- all_genecard_name  %>% distinct(Symbol, .keep_all = TRUE)

# Merge with summary
merged_all_genecards <- all_genecard_summary %>%
  left_join(alias_unique, by = "Symbol") %>%
  left_join(name_unique, by = "Symbol")

# Check result
head(merged_all_genecards)
nrow(merged_all_genecards)


rm(merged_combined)

# Keep only one row per Symbol in the genecards table
merged_all_genecards_unique <- merged_all_genecards %>%
  group_by(Symbol) %>%
  slice(1) %>%  # take the first row per gene
  ungroup()

# Now join
merged_combined <- merged_df_fixed %>%
  left_join(merged_all_genecards_unique, by = c("gene_symbol" = "Symbol"))



write.csv(merged_combined, "All_Disgenet_with_gene_cards.csv")







