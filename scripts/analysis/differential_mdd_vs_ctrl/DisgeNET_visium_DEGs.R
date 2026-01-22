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


V1_0 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "0"))

# First run 
V1_0_genes <- V1_0$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_0_genes1 <- V1_0_genes[1:100]
V1_0_genes2 <- V1_0_genes[101:200]
V1_0_genes3 <- V1_0_genes[201:300]
V1_0_genes4 <- V1_0_genes[301:400]
V1_0_genes5 <- V1_0_genes[401:500]
V1_0_genes6 <- V1_0_genes[501:600]
V1_0_genes7 <- V1_0_genes[601:700]
V1_0_genes8 <- V1_0_genes[701:800]
V1_0_genes9 <- V1_0_genes[801:900]
V1_0_genes10 <- V1_0_genes[901:1000]
V1_0_genes11 <- V1_0_genes[1001:1100]
V1_0_genes12 <- V1_0_genes[1101:1200]
V1_0_genes13 <- V1_0_genes[1201:1300]
V1_0_genes14 <- V1_0_genes[1301:1400]
V1_0_genes15 <- V1_0_genes[1401:1500]
V1_0_genes16 <- V1_0_genes[1501:1600]
V1_0_genes17 <- V1_0_genes[1601:1700]
V1_0_genes18 <- V1_0_genes[1701:1800]
V1_0_genes19 <- V1_0_genes[1801:1900]
V1_0_genes20 <- V1_0_genes[1901:2000]
V1_0_genes21 <- V1_0_genes[2001:2100]
V1_0_genes22 <- V1_0_genes[2101:2200]
V1_0_genes23 <- V1_0_genes[2201:2300]
V1_0_genes24 <- V1_0_genes[2301:2368]

V1_0_results1 <- gene2disease(
  gene     = V1_0_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results2 <- gene2disease(
  gene     = V1_0_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results3 <- gene2disease(
  gene     = V1_0_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_0_results4 <- gene2disease(
  gene     = V1_0_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results5 <- gene2disease(
  gene     = V1_0_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results6 <- gene2disease(
  gene     = V1_0_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results7 <- gene2disease(
  gene     = V1_0_genes7,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results8 <- gene2disease(
  gene     = V1_0_genes8,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results9 <- gene2disease(
  gene     = V1_0_genes9,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results10 <- gene2disease(
  gene     = V1_0_genes10,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results11 <- gene2disease(
  gene     = V1_0_genes11,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results12 <- gene2disease(
  gene     = V1_0_genes12,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results13 <- gene2disease(
  gene     = V1_0_genes13,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results14 <- gene2disease(
  gene     = V1_0_genes14,
  database = "PSYGENET",
  verbose  = TRUE
)



V1_0_results15 <- gene2disease(
  gene     = V1_0_genes15,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results16 <- gene2disease(
  gene     = V1_0_genes16,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results17 <- gene2disease(
  gene     = V1_0_genes17,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results18 <- gene2disease(
  gene     = V1_0_genes18,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results19 <- gene2disease(
  gene     = V1_0_genes19,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results20 <- gene2disease(
  gene     = V1_0_genes20,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results21 <- gene2disease(
  gene     = V1_0_genes21,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results22 <- gene2disease(
  gene     = V1_0_genes22,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results23 <- gene2disease(
  gene     = V1_0_genes23,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results24 <- gene2disease(
  gene     = V1_0_genes24,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_0_results_tab1 <- V1_0_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab2 <- V1_0_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab3 <- V1_0_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab4 <- V1_0_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab5 <- V1_0_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab6 <- V1_0_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab7 <- V1_0_results7@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab8 <- V1_0_results8@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab9 <- V1_0_results9@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab10 <- V1_0_results10@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab11 <- V1_0_results11@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab12 <- V1_0_results12@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab13 <- V1_0_results13@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab14 <- V1_0_results14@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab15 <- V1_0_results15@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab16 <- V1_0_results16@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab17 <- V1_0_results17@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab18 <- V1_0_results18@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab19 <- V1_0_results19@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab20 <- V1_0_results20@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab21 <- V1_0_results21@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab22 <- V1_0_results22@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab23 <- V1_0_results23@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_0_results_tab24 <- V1_0_results24@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)





V1_0_all_tab <- unique(c(
  V1_0_results_tab1$gene_symbol, V1_0_results_tab2$gene_symbol, V1_0_results_tab3$gene_symbol,
  V1_0_results_tab4$gene_symbol, V1_0_results_tab5$gene_symbol, V1_0_results_tab6$gene_symbol,
  V1_0_results_tab7$gene_symbol, V1_0_results_tab8$gene_symbol, V1_0_results_tab9$gene_symbol,
  V1_0_results_tab10$gene_symbol, V1_0_results_tab11$gene_symbol, V1_0_results_tab12$gene_symbol,
  V1_0_results_tab13$gene_symbol, V1_0_results_tab14$gene_symbol, V1_0_results_tab15$gene_symbol,
  V1_0_results_tab16$gene_symbol, V1_0_results_tab17$gene_symbol, V1_0_results_tab18$gene_symbol,
  V1_0_results_tab19$gene_symbol, V1_0_results_tab20$gene_symbol, V1_0_results_tab21$gene_symbol,
  V1_0_results_tab22$gene_symbol, V1_0_results_tab23$gene_symbol, V1_0_results_tab24$gene_symbol
  
))

V1_0_all_tab1 <- V1_0_all_tab[1:100]
V1_0_all_tab2 <- V1_0_all_tab[101:200]
V1_0_all_tab3 <- V1_0_all_tab[201:205]

# Rerun with just this list
V1_0_results_all_1 <- gene2disease(
  gene     = V1_0_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results_all_2 <- gene2disease(
  gene     = V1_0_all_tab2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_0_results_all_3 <- gene2disease(
  gene     = V1_0_all_tab3,
  database = "PSYGENET",
  verbose  = TRUE
)



# Use your DisGeNET results table
df1 <- V1_0_results_all_1@qresult
df2 <- V1_0_results_all_2@qresult
df3 <- V1_0_results_all_3@qresult

library(dplyr)

merged_df <- bind_rows(df1, df2, df3)

# Optional: filter by score or number of supporting PMIDs
df_filtered <- merged_df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(merged_df, file = "PsyGeNET_results_V1_0_11.3.25.xlsx")


#Cluster 1 
###############################
readxl::excel_sheets("/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/v1_degs_significant_no_sex_10_23_2025.xlsx")

V1_1 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "1"))

# First run 
V1_1_genes <- V1_1$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_1_genes1 <- V1_1_genes[1:100]
V1_1_genes2 <- V1_1_genes[101:162]

V1_1_results1 <- gene2disease(
  gene     = V1_1_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_1_results2 <- gene2disease(
  gene     = V1_1_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)




V1_1_results_tab1 <- V1_1_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_1_results_tab2 <- V1_1_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)



V1_1_all_tab <- unique(c(
  V1_1_results_tab1$gene_symbol, V1_1_results_tab2$gene_symbol))

# Rerun with just this list
V1_1_results_all <- gene2disease(
  gene     = V1_1_all_tab,
  database = "PSYGENET",
  verbose  = TRUE
)

df <- V1_1_results_all@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.5)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_1_11.3.25.xlsx")


#Cluster 2 
######################################
V1_2 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "2"))

# First run 
V1_2_genes <- V1_2$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_2_genes1 <- V1_2_genes[1:100]
V1_2_genes2 <- V1_2_genes[101:200]
V1_2_genes3 <- V1_2_genes[201:300]
V1_2_genes4 <- V1_2_genes[301:400]
V1_2_genes5 <- V1_2_genes[401:500]
V1_2_genes6 <- V1_2_genes[501:600]
V1_2_genes7 <- V1_2_genes[601:700]
V1_2_genes8 <- V1_2_genes[701:800]
V1_2_genes9 <- V1_2_genes[801:900]
V1_2_genes10 <- V1_2_genes[901:1000]
V1_2_genes11 <- V1_2_genes[1001:1100]
V1_2_genes12 <- V1_2_genes[1101:1200]
V1_2_genes13 <- V1_2_genes[1201:1300]
V1_2_genes14 <- V1_2_genes[1301:1400]
V1_2_genes15 <- V1_2_genes[1401:1500]
V1_2_genes16 <- V1_2_genes[1501:1600]
V1_2_genes17 <- V1_2_genes[1601:1700]
V1_2_genes18 <- V1_2_genes[1701:1800]
V1_2_genes19 <- V1_2_genes[1801:1852]

V1_2_results1 <- gene2disease(
  gene     = V1_2_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results2 <- gene2disease(
  gene     = V1_2_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results3 <- gene2disease(
  gene     = V1_2_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_2_results4 <- gene2disease(
  gene     = V1_2_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results5 <- gene2disease(
  gene     = V1_2_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results6 <- gene2disease(
  gene     = V1_2_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results7 <- gene2disease(
  gene     = V1_2_genes7,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results8 <- gene2disease(
  gene     = V1_2_genes8,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results9 <- gene2disease(
  gene     = V1_2_genes9,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results10 <- gene2disease(
  gene     = V1_2_genes10,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results11 <- gene2disease(
  gene     = V1_2_genes11,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results12 <- gene2disease(
  gene     = V1_2_genes12,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results13 <- gene2disease(
  gene     = V1_2_genes13,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results14 <- gene2disease(
  gene     = V1_2_genes14,
  database = "PSYGENET",
  verbose  = TRUE
)



V1_2_results15 <- gene2disease(
  gene     = V1_2_genes15,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results16 <- gene2disease(
  gene     = V1_2_genes16,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results17 <- gene2disease(
  gene     = V1_2_genes17,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results18 <- gene2disease(
  gene     = V1_2_genes18,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results19 <- gene2disease(
  gene     = V1_2_genes19,
  database = "PSYGENET",
  verbose  = TRUE
)





V1_2_results_tab1 <- V1_2_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab2 <- V1_2_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab3 <- V1_2_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab4 <- V1_2_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab5 <- V1_2_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab6 <- V1_2_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab7 <- V1_2_results7@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab8 <- V1_2_results8@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab9 <- V1_2_results9@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab10 <- V1_2_results10@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab11 <- V1_2_results11@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab12 <- V1_2_results12@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab13 <- V1_2_results13@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab14 <- V1_2_results14@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab15 <- V1_2_results15@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab16 <- V1_2_results16@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab17 <- V1_2_results17@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab18 <- V1_2_results18@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_2_results_tab19 <- V1_2_results19@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)


V1_2_all_tab <- unique(c(
  V1_2_results_tab1$gene_symbol, V1_2_results_tab2$gene_symbol, V1_2_results_tab3$gene_symbol,
  V1_2_results_tab4$gene_symbol, V1_2_results_tab5$gene_symbol, V1_2_results_tab6$gene_symbol,
  V1_2_results_tab7$gene_symbol, V1_2_results_tab8$gene_symbol, V1_2_results_tab9$gene_symbol,
  V1_2_results_tab10$gene_symbol, V1_2_results_tab11$gene_symbol, V1_2_results_tab12$gene_symbol,
  V1_2_results_tab13$gene_symbol, V1_2_results_tab14$gene_symbol, V1_2_results_tab15$gene_symbol,
  V1_2_results_tab16$gene_symbol, V1_2_results_tab17$gene_symbol, V1_2_results_tab18$gene_symbol,
  V1_2_results_tab19$gene_symbol))

V1_2_all_tab1 <- V1_2_all_tab[1:100]
V1_2_all_tab2 <- V1_2_all_tab[101:151]

# Rerun with just this list
V1_2_results_all_1 <- gene2disease(
  gene     = V1_2_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_2_results_all_2 <- gene2disease(
  gene     = V1_2_all_tab2,
  database = "PSYGENET",
  verbose  = TRUE
)




# Use your DisGeNET results table
df1 <- V1_2_results_all_1@qresult
df2 <- V1_2_results_all_2@qresult


library(dplyr)

merged_df <- bind_rows(df1, df2)

# Optional: filter by score or number of supporting PMIDs
df_filtered <- merged_df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(merged_df, file = "PsyGeNET_results_V1_2_11.3.25.xlsx")



#Cluster 3 
######################################

V1_3 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "3"))

# First run 
V1_3_genes <- V1_3$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_3_genes1 <- V1_3_genes[1:60]

V1_3_results1 <- gene2disease(
  gene     = V1_3_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)



df <- V1_3_results1@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.5)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_3_11.3.25.xlsx")



#Cluster 4 
#############################################
V1_4 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "4"))


# First run 
V1_4_genes <- V1_4$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_4_genes1 <- V1_4_genes[1:100]
V1_4_genes2 <- V1_4_genes[101:200]
V1_4_genes3 <- V1_4_genes[201:300]
V1_4_genes4 <- V1_4_genes[301:400]
V1_4_genes5 <- V1_4_genes[401:500]
V1_4_genes6 <- V1_4_genes[501:600]
V1_4_genes7 <- V1_4_genes[601:700]
V1_4_genes8 <- V1_4_genes[701:800]
V1_4_genes9 <- V1_4_genes[801:900]
V1_4_genes10 <- V1_4_genes[901:1000]
V1_4_genes11 <- V1_4_genes[1001:1100]
V1_4_genes12 <- V1_4_genes[1101:1200]
V1_4_genes13 <- V1_4_genes[1201:1278]



V1_4_results1 <- gene2disease(
  gene     = V1_4_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results2 <- gene2disease(
  gene     = V1_4_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results3 <- gene2disease(
  gene     = V1_4_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_4_results4 <- gene2disease(
  gene     = V1_4_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results5 <- gene2disease(
  gene     = V1_4_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results6 <- gene2disease(
  gene     = V1_4_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results7 <- gene2disease(
  gene     = V1_4_genes7,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results8 <- gene2disease(
  gene     = V1_4_genes8,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results9 <- gene2disease(
  gene     = V1_4_genes9,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results10 <- gene2disease(
  gene     = V1_4_genes10,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results11 <- gene2disease(
  gene     = V1_4_genes11,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results12 <- gene2disease(
  gene     = V1_4_genes12,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results13 <- gene2disease(
  gene     = V1_4_genes13,
  database = "PSYGENET",
  verbose  = TRUE
)




V1_4_results_tab1 <- V1_4_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab2 <- V1_4_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab3 <- V1_4_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab4 <- V1_4_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab5 <- V1_4_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab6 <- V1_4_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab7 <- V1_4_results7@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab8 <- V1_4_results8@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab9 <- V1_4_results9@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab10 <- V1_4_results10@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab11 <- V1_4_results11@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab12 <- V1_4_results12@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_4_results_tab13 <- V1_4_results13@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)





V1_4_all_tab <- unique(c(
  V1_4_results_tab1$gene_symbol, V1_4_results_tab2$gene_symbol, V1_4_results_tab3$gene_symbol,
  V1_4_results_tab4$gene_symbol, V1_4_results_tab5$gene_symbol, V1_4_results_tab6$gene_symbol,
  V1_4_results_tab7$gene_symbol, V1_4_results_tab8$gene_symbol, V1_4_results_tab9$gene_symbol,
  V1_4_results_tab10$gene_symbol, V1_4_results_tab11$gene_symbol, V1_4_results_tab12$gene_symbol,
  V1_4_results_tab13$gene_symbol
))

V1_4_all_tab1 <- V1_4_all_tab[1:100]
V1_4_all_tab2 <- V1_4_all_tab[101:111]

# Rerun with just this list
V1_4_results_all_1 <- gene2disease(
  gene     = V1_4_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_4_results_all_2 <- gene2disease(
  gene     = V1_4_all_tab2,
  database = "PSYGENET",
  verbose  = TRUE
)


# Use your DisGeNET results table
df1 <- V1_4_results_all_1@qresult
df2 <- V1_4_results_all_2@qresult


library(dplyr)

merged_df <- bind_rows(df1, df2)

# Optional: filter by score or number of supporting PMIDs
df_filtered <- merged_df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(merged_df, file = "PsyGeNET_results_V1_4_11.3.25.xlsx")




#Cluster 5 
#####################################
V1_5 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "5"))

# First run 
V1_5_genes <- V1_5$gene


# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_5_genes1 <- V1_5_genes[1:100]
V1_5_genes2 <- V1_5_genes[101:200]
V1_5_genes3 <- V1_5_genes[201:300]
V1_5_genes4 <- V1_5_genes[301:397]

V1_5_results1 <- gene2disease(
  gene     = V1_5_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_5_results2 <- gene2disease(
  gene     = V1_5_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_5_results3 <- gene2disease(
  gene     = V1_5_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_5_results4 <- gene2disease(
  gene     = V1_5_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)





V1_5_results_tab1 <- V1_5_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_5_results_tab2 <- V1_5_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_5_results_tab3 <- V1_5_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_5_results_tab4 <- V1_5_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)





V1_5_all_tab <- unique(c(
  V1_5_results_tab1$gene_symbol, V1_5_results_tab2$gene_symbol, V1_5_results_tab3$gene_symbol,
  V1_5_results_tab4$gene_symbol))



V1_5_all_tab1 <- V1_5_all_tab[1:54]

# Rerun with just this list
V1_5_results_all_1 <- gene2disease(
  gene     = V1_5_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)




# Use your DisGeNET results table
df <- V1_5_results_all_1@qresult


# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_5_11.3.25.xlsx")





#Cluster 6 
##################################
V1_6 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "6"))


# First run 
V1_6_genes <- V1_6$gene


# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_6_genes1 <- V1_6_genes[1:100]
V1_6_genes2 <- V1_6_genes[101:104]


V1_6_results1 <- gene2disease(
  gene     = V1_6_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_6_results2 <- gene2disease(
  gene     = V1_6_genes2,
  database = "PSYGENET",
  verbose  = TRUE
) #no results


df <- V1_6_results1@qresult


# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.5)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_6_11.3.25.xlsx")




#Cluster 7 
#########################################
V1_7 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "7"))


# First run 
V1_7_genes <- V1_7$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_7_genes1 <- V1_7_genes[1:100]
V1_7_genes2 <- V1_7_genes[101:200]
V1_7_genes3 <- V1_7_genes[201:300]
V1_7_genes4 <- V1_7_genes[301:400]
V1_7_genes5 <- V1_7_genes[401:500]
V1_7_genes6 <- V1_7_genes[501:516]

V1_7_results1 <- gene2disease(
  gene     = V1_7_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_7_results2 <- gene2disease(
  gene     = V1_7_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_7_results3 <- gene2disease(
  gene     = V1_7_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_7_results4 <- gene2disease(
  gene     = V1_7_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_7_results5 <- gene2disease(
  gene     = V1_7_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_7_results6 <- gene2disease(
  gene     = V1_7_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)








V1_7_results_tab1 <- V1_7_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_7_results_tab2 <- V1_7_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_7_results_tab3 <- V1_7_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_7_results_tab4 <- V1_7_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_7_results_tab5 <- V1_7_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_7_results_tab6 <- V1_7_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)





V1_7_all_tab <- unique(c(
  V1_7_results_tab1$gene_symbol, V1_7_results_tab2$gene_symbol, V1_7_results_tab3$gene_symbol,
  V1_7_results_tab4$gene_symbol, V1_7_results_tab5$gene_symbol, V1_7_results_tab6$gene_symbol
))

V1_7_all_tab1 <- V1_7_all_tab[1:48]

# Rerun with just this list
V1_7_results_all_1 <- gene2disease(
  gene     = V1_7_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)



# Use your DisGeNET results table
df <- V1_7_results_all_1@qresult


# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_7_11.3.25.xlsx")




#Cluster 8 
##############################
V1_8 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "8"))

# First run 
V1_8_genes <- V1_8$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_8_genes1 <- V1_8_genes[1:100]
V1_8_genes2 <- V1_8_genes[101:200]
V1_8_genes3 <- V1_8_genes[201:300]
V1_8_genes4 <- V1_8_genes[301:400]
V1_8_genes5 <- V1_8_genes[401:500]
V1_8_genes6 <- V1_8_genes[501:600]
V1_8_genes7 <- V1_8_genes[601:666]



V1_8_results1 <- gene2disease(
  gene     = V1_8_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_8_results2 <- gene2disease(
  gene     = V1_8_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_8_results3 <- gene2disease(
  gene     = V1_8_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_8_results4 <- gene2disease(
  gene     = V1_8_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_8_results5 <- gene2disease(
  gene     = V1_8_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_8_results6 <- gene2disease(
  gene     = V1_8_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_8_results7 <- gene2disease(
  gene     = V1_8_genes7,
  database = "PSYGENET",
  verbose  = TRUE
)



V1_8_results_tab1 <- V1_8_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_8_results_tab2 <- V1_8_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_8_results_tab3 <- V1_8_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_8_results_tab4 <- V1_8_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_8_results_tab5 <- V1_8_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_8_results_tab6 <- V1_8_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_8_results_tab7 <- V1_8_results7@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)




V1_8_all_tab <- unique(c(
  V1_8_results_tab1$gene_symbol, V1_8_results_tab2$gene_symbol, V1_8_results_tab3$gene_symbol,
  V1_8_results_tab4$gene_symbol, V1_8_results_tab5$gene_symbol, V1_8_results_tab6$gene_symbol,
  V1_8_results_tab7$gene_symbol
))

V1_8_all_tab1 <- V1_8_all_tab[1:74]

# Rerun with just this list
V1_8_results_all_1 <- gene2disease(
  gene     = V1_8_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)


# Use your DisGeNET results table
df <- V1_8_results_all_1@qresult


# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_8_11.3.25.xlsx")




#Cluster 9 only 4 DEGS no psygenet result
########################################
V1_9 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "9"))


# First run 
V1_9_genes <- V1_9$gene

V1_9_results <- gene2disease(
  gene     = V1_9_genes,
  database = "PSYGENET",
  verbose  = TRUE
)



#Cluster 10
######################################
V1_10 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "10"))


# First run 
V1_10_genes <- V1_10$gene


# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
V1_10_genes1 <- V1_10_genes[1:100]
V1_10_genes2 <- V1_10_genes[101:200]
V1_10_genes3 <- V1_10_genes[201:300]
V1_10_genes4 <- V1_10_genes[301:400]
V1_10_genes5 <- V1_10_genes[401:500]
V1_10_genes6 <- V1_10_genes[501:600]
V1_10_genes7 <- V1_10_genes[601:700]
V1_10_genes8 <- V1_10_genes[701:800]
V1_10_genes9 <- V1_10_genes[801:900]
V1_10_genes10 <- V1_10_genes[901:944]


V1_10_results1 <- gene2disease(
  gene     = V1_10_genes1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results2 <- gene2disease(
  gene     = V1_10_genes2,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results3 <- gene2disease(
  gene     = V1_10_genes3,
  database = "PSYGENET",
  verbose  = TRUE
)


V1_10_results4 <- gene2disease(
  gene     = V1_10_genes4,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results5 <- gene2disease(
  gene     = V1_10_genes5,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results6 <- gene2disease(
  gene     = V1_10_genes6,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results7 <- gene2disease(
  gene     = V1_10_genes7,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results8 <- gene2disease(
  gene     = V1_10_genes8,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results9 <- gene2disease(
  gene     = V1_10_genes9,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results10 <- gene2disease(
  gene     = V1_10_genes10,
  database = "PSYGENET",
  verbose  = TRUE
)




V1_10_results_tab1 <- V1_10_results1@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab2 <- V1_10_results2@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab3 <- V1_10_results3@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab4 <- V1_10_results4@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab5 <- V1_10_results5@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab6 <- V1_10_results6@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab7 <- V1_10_results7@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab8 <- V1_10_results8@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab9 <- V1_10_results9@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)

V1_10_results_tab10 <- V1_10_results10@qresult[  ,c("gene_symbol", "disease_name","score", "yearInitial", "yearFinal")]  %>% unique()  %>%
  arrange(desc(score), yearInitial)



V1_10_all_tab <- unique(c(
  V1_10_results_tab1$gene_symbol, V1_10_results_tab2$gene_symbol, V1_10_results_tab3$gene_symbol,
  V1_10_results_tab4$gene_symbol, V1_10_results_tab5$gene_symbol, V1_10_results_tab6$gene_symbol,
  V1_10_results_tab7$gene_symbol, V1_10_results_tab8$gene_symbol, V1_10_results_tab9$gene_symbol,
  V1_10_results_tab10$gene_symbol
))

V1_10_all_tab1 <- V1_10_all_tab[1:100]
V1_10_all_tab2 <- V1_10_all_tab[101:103]

# Rerun with just this list
V1_10_results_all_1 <- gene2disease(
  gene     = V1_10_all_tab1,
  database = "PSYGENET",
  verbose  = TRUE
)

V1_10_results_all_2 <- gene2disease(
  gene     = V1_10_all_tab2,
  database = "PSYGENET",
  verbose  = TRUE
)


# Use your DisGeNET results table
df1 <- V1_10_results_all_1@qresult
df2 <- V1_10_results_all_2@qresult


library(dplyr)

merged_df <- bind_rows(df1, df2)

# Optional: filter by score or number of supporting PMIDs
df_filtered <- merged_df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(merged_df, file = "PsyGeNET_results_V1_10_11.3.25.xlsx")


#Cluster 11
######################################
V1_11 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "11"))

# First run 
V1_11_genes <- V1_11$gene



V1_11_results <- gene2disease(
  gene     = V1_11_genes,
  database = "PSYGENET",
  verbose  = TRUE
)




df <- V1_11_results@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.5)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_11_11.3.25.xlsx")



#Cluster 12
#############################################
V1_12 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "12"))


# First run 
V1_12_genes <- V1_12$gene

V1_12_results <- gene2disease(
  gene     = V1_12_genes,
  database = "PSYGENET",
  verbose  = TRUE
)



df <- V1_12_results@qresult

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.5)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V1_12_11.3.25.xlsx")


#Cluster 13 
#######################################
V1_13 <- as.data.frame(readxl::read_excel(
  "~/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v1_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "13"))

# First run 
V1_13_genes <- V1_13$gene. #0 genes



#V2 Cliuster 0 
#########################################*
V2_0 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "0"))

# First run 
V2_0_genes <- V2_0$gene


# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_0_gene_chunks <- split(
  V2_0_genes,
  ceiling(seq_along(V2_0_genes) / chunk_size)
)

# Assign names like V2_0_genes1, V2_0_genes2, ...
for (i in seq_along(V2_0_gene_chunks)) {
  assign(paste0("V2_0_genes", i), V2_0_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_0_gene_chunks)

# Example: view first few
head(V2_0_genes1)
head(V2_0_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_0_gene_chunks <- split(V2_0_genes, ceiling(seq_along(V2_0_genes) / chunk_size))

# Create an empty list to store results
V2_0_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_0_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_0_gene_chunks), "\n")
  
  V2_0_results_list[[i]] <- gene2disease(
    gene     = V2_0_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_0_qresults_list <- lapply(V2_0_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_0_results_all <- do.call(rbind, V2_0_qresults_list)

# Remove duplicates (optional)
V2_0_results_all <- unique(V2_0_results_all)

# Second run 
V2_0_genes_tab <- V2_0_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_0_gene_chunks <- split(
  V2_0_genes_tab,
  ceiling(seq_along(V2_0_genes_tab) / chunk_size)
)

# Assign names like V2_0_genes1, V2_0_genes2, ...
for (i in seq_along(V2_0_gene_chunks)) {
  assign(paste0("V2_0_genes_tab", i), V2_0_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_0_gene_chunks)

# Example: view first few
head(V2_0_genes1)
head(V2_0_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_0_gene_chunks <- split(V2_0_genes_tab, ceiling(seq_along(V2_0_genes_tab) / chunk_size))

# Create an empty list to store results
V2_0_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_0_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_0_gene_chunks), "\n")
  
  V2_0_results_list_tab[[i]] <- gene2disease(
    gene     = V2_0_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_0_qresults_list <- lapply(V2_0_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_0_results_all <- do.call(rbind, V2_0_qresults_list)

# Remove duplicates (optional)
V2_0_results_all <- unique(V2_0_results_all)


df <- V2_0_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

# 1. Extract edges (gene–disease associations)
edges <- df_filtered %>%
  dplyr::select(gene_symbol, disease_name) %>%
  distinct() %>%
  rename(from = gene_symbol, to = disease_name)

# 2. Build igraph object
net <- graph_from_data_frame(edges, directed = FALSE)

# 3. Assign node type (Gene vs Disease)
V(net)$type <- ifelse(
  V(net)$name %in% V2_0_results_all@qresult$gene_symbol & 
    V(net)$name %in% V2_0_results_all@qresult$gene_symbol,
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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_0_10.30.25.xlsx")



#V2 Cliuster 1 
#########################################*
V2_1 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "1"))

# First run 
V2_1_genes <- V2_1$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_1_gene_chunks <- split(
  V2_1_genes,
  ceiling(seq_along(V2_1_genes) / chunk_size)
)

# Assign names like V2_1_genes1, V2_1_genes2, ...
for (i in seq_along(V2_1_gene_chunks)) {
  assign(paste0("V2_1_genes", i), V2_1_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_1_gene_chunks)

# Example: view first few
head(V2_1_genes1)
head(V2_1_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_1_gene_chunks <- split(V2_1_genes, ceiling(seq_along(V2_1_genes) / chunk_size))

# Create an empty list to store results
V2_1_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_1_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_1_gene_chunks), "\n")
  
  V2_1_results_list[[i]] <- gene2disease(
    gene     = V2_1_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_1_qresults_list <- lapply(V2_1_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_1_results_all <- do.call(rbind, V2_1_qresults_list)

# Remove duplicates (optional)
V2_1_results_all <- unique(V2_1_results_all)

# Second run 
V2_1_genes_tab <- V2_1_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_1_gene_chunks <- split(
  V2_1_genes_tab,
  ceiling(seq_along(V2_1_genes_tab) / chunk_size)
)

# Assign names like V2_1_genes1, V2_1_genes2, ...
for (i in seq_along(V2_1_gene_chunks)) {
  assign(paste0("V2_1_genes_tab", i), V2_1_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_1_gene_chunks)

# Example: view first few
head(V2_1_genes1)
head(V2_1_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_1_gene_chunks <- split(V2_1_genes_tab, ceiling(seq_along(V2_1_genes_tab) / chunk_size))

# Create an empty list to store results
V2_1_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_1_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_1_gene_chunks), "\n")
  
  V2_1_results_list_tab[[i]] <- gene2disease(
    gene     = V2_1_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_1_qresults_list <- lapply(V2_1_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_1_results_all <- do.call(rbind, V2_1_qresults_list)

# Remove duplicates (optional)
V2_1_results_all <- unique(V2_1_results_all)


df <- V2_1_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_1_11.3.25.xlsx")


#V2 Cliuster 2
#########################################*
V2_2 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "2"))

# First run 
V2_2_genes <- V2_2$gene

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_2_gene_chunks <- split(
  V2_2_genes,
  ceiling(seq_along(V2_2_genes) / chunk_size)
)

# Assign names like V2_2_genes1, V2_2_genes2, ...
for (i in seq_along(V2_2_gene_chunks)) {
  assign(paste0("V2_2_genes", i), V2_2_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_2_gene_chunks)

# Example: view first few
head(V2_2_genes1)
head(V2_2_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_2_gene_chunks <- split(V2_2_genes, ceiling(seq_along(V2_2_genes) / chunk_size))

# Create an empty list to store results
V2_2_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_2_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_2_gene_chunks), "\n")
  
  V2_2_results_list[[i]] <- gene2disease(
    gene     = V2_2_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_2_qresults_list <- lapply(V2_2_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_2_results_all <- do.call(rbind, V2_2_qresults_list)

# Remove duplicates (optional)
V2_2_results_all <- unique(V2_2_results_all)

# Second run 
V2_2_genes_tab <- V2_2_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_2_gene_chunks <- split(
  V2_2_genes_tab,
  ceiling(seq_along(V2_2_genes_tab) / chunk_size)
)

# Assign names like V2_2_genes1, V2_2_genes2, ...
for (i in seq_along(V2_2_gene_chunks)) {
  assign(paste0("V2_2_genes_tab", i), V2_2_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_2_gene_chunks)

# Example: view first few
head(V2_2_genes1)
head(V2_2_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_2_gene_chunks <- split(V2_2_genes_tab, ceiling(seq_along(V2_2_genes_tab) / chunk_size))

# Create an empty list to store results
V2_2_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_2_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_2_gene_chunks), "\n")
  
  V2_2_results_list_tab[[i]] <- gene2disease(
    gene     = V2_2_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_2_qresults_list <- lapply(V2_2_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_2_results_all <- do.call(rbind, V2_2_qresults_list)

# Remove duplicates (optional)
V2_2_results_all <- unique(V2_2_results_all)


df <- V2_2_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_2_11.3.25.xlsx")



#V2 Cliuster 3
#########################################*
V2_3 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "3"))


# First run 
V2_3_genes <- V2_3$gene

length(V2_3_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_3_gene_chunks <- split(
  V2_3_genes,
  ceiling(seq_along(V2_3_genes) / chunk_size)
)

# Assign names like V2_3_genes1, V2_3_genes2, ...
for (i in seq_along(V2_3_gene_chunks)) {
  assign(paste0("V2_3_genes", i), V2_3_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_3_gene_chunks)

# Example: view first few
head(V2_3_genes1)
head(V2_3_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_3_gene_chunks <- split(V2_3_genes, ceiling(seq_along(V2_3_genes) / chunk_size))

# Create an empty list to store results
V2_3_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_3_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_3_gene_chunks), "\n")
  
  V2_3_results_list[[i]] <- gene2disease(
    gene     = V2_3_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_3_qresults_list <- lapply(V2_3_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_3_results_all <- do.call(rbind, V2_3_qresults_list)

# Remove duplicates (optional)
V2_3_results_all <- unique(V2_3_results_all)

# Second run 
V2_3_genes_tab <- V2_3_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_3_gene_chunks <- split(
  V2_3_genes_tab,
  ceiling(seq_along(V2_3_genes_tab) / chunk_size)
)

# Assign names like V2_3_genes1, V2_3_genes2, ...
for (i in seq_along(V2_3_gene_chunks)) {
  assign(paste0("V2_3_genes_tab", i), V2_3_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_3_gene_chunks)

# Example: view first few
head(V2_3_genes1)
head(V2_3_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_3_gene_chunks <- split(V2_3_genes_tab, ceiling(seq_along(V2_3_genes_tab) / chunk_size))

# Create an empty list to store results
V2_3_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_3_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_3_gene_chunks), "\n")
  
  V2_3_results_list_tab[[i]] <- gene2disease(
    gene     = V2_3_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_3_qresults_list <- lapply(V2_3_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_3_results_all <- do.call(rbind, V2_3_qresults_list)

# Remove duplicates (optional)
V2_3_results_all <- unique(V2_3_results_all)


df <- V2_3_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_3_11.3.25.xlsx")


#V2 Cliuster 4
#########################################*
V2_4 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "4"))


# First run 
V2_4_genes <- V2_4$gene

length(V2_4_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_4_gene_chunks <- split(
  V2_4_genes,
  ceiling(seq_along(V2_4_genes) / chunk_size)
)

# Assign names like V2_4_genes1, V2_4_genes2, ...
for (i in seq_along(V2_4_gene_chunks)) {
  assign(paste0("V2_4_genes", i), V2_4_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_4_gene_chunks)

# Example: view first few
head(V2_4_genes1)
head(V2_4_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_4_gene_chunks <- split(V2_4_genes, ceiling(seq_along(V2_4_genes) / chunk_size))

# Create an empty list to store results
V2_4_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_4_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_4_gene_chunks), "\n")
  
  V2_4_results_list[[i]] <- gene2disease(
    gene     = V2_4_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_4_qresults_list <- lapply(V2_4_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_4_results_all <- do.call(rbind, V2_4_qresults_list)

# Remove duplicates (optional)
V2_4_results_all <- unique(V2_4_results_all)

# Second run 
V2_4_genes_tab <- V2_4_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_4_gene_chunks <- split(
  V2_4_genes_tab,
  ceiling(seq_along(V2_4_genes_tab) / chunk_size)
)

# Assign names like V2_4_genes1, V2_4_genes2, ...
for (i in seq_along(V2_4_gene_chunks)) {
  assign(paste0("V2_4_genes_tab", i), V2_4_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_4_gene_chunks)

# Example: view first few
head(V2_4_genes1)
head(V2_4_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_4_gene_chunks <- split(V2_4_genes_tab, ceiling(seq_along(V2_4_genes_tab) / chunk_size))

# Create an empty list to store results
V2_4_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_4_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_4_gene_chunks), "\n")
  
  V2_4_results_list_tab[[i]] <- gene2disease(
    gene     = V2_4_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_4_qresults_list <- lapply(V2_4_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_4_results_all <- do.call(rbind, V2_4_qresults_list)

# Remove duplicates (optional)
V2_4_results_all <- unique(V2_4_results_all)


df <- V2_4_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_4_11.3.25.xlsx")


#V2 Cliuster 5
#########################################*
V2_5 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "5"))


# First run 
V2_5_genes <- V2_5$gene

length(V2_5_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_5_gene_chunks <- split(
  V2_5_genes,
  ceiling(seq_along(V2_5_genes) / chunk_size)
)

# Assign names like V2_5_genes1, V2_5_genes2, ...
for (i in seq_along(V2_5_gene_chunks)) {
  assign(paste0("V2_5_genes", i), V2_5_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_5_gene_chunks)

# Example: view first few
head(V2_5_genes1)
head(V2_5_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_5_gene_chunks <- split(V2_5_genes, ceiling(seq_along(V2_5_genes) / chunk_size))

# Create an empty list to store results
V2_5_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_5_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_5_gene_chunks), "\n")
  
  V2_5_results_list[[i]] <- gene2disease(
    gene     = V2_5_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_5_qresults_list <- lapply(V2_5_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_5_results_all <- do.call(rbind, V2_5_qresults_list)

# Remove duplicates (optional)
V2_5_results_all <- unique(V2_5_results_all)

# Second run 
V2_5_genes_tab <- V2_5_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_5_gene_chunks <- split(
  V2_5_genes_tab,
  ceiling(seq_along(V2_5_genes_tab) / chunk_size)
)

# Assign names like V2_5_genes1, V2_5_genes2, ...
for (i in seq_along(V2_5_gene_chunks)) {
  assign(paste0("V2_5_genes_tab", i), V2_5_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_5_gene_chunks)

# Example: view first few
head(V2_5_genes1)
head(V2_5_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_5_gene_chunks <- split(V2_5_genes_tab, ceiling(seq_along(V2_5_genes_tab) / chunk_size))

# Create an empty list to store results
V2_5_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_5_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_5_gene_chunks), "\n")
  
  V2_5_results_list_tab[[i]] <- gene2disease(
    gene     = V2_5_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_5_qresults_list <- lapply(V2_5_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_5_results_all <- do.call(rbind, V2_5_qresults_list)

# Remove duplicates (optional)
V2_5_results_all <- unique(V2_5_results_all)


df <- V2_5_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_5_11.3.25.xlsx")


#V2 Cliuster 6
#########################################*
V2_6 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "6"))


# First run 
V2_6_genes <- V2_6$gene

length(V2_6_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_6_gene_chunks <- split(
  V2_6_genes,
  ceiling(seq_along(V2_6_genes) / chunk_size)
)

# Assign names like V2_6_genes1, V2_6_genes2, ...
for (i in seq_along(V2_6_gene_chunks)) {
  assign(paste0("V2_6_genes", i), V2_6_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_6_gene_chunks)

# Example: view first few
head(V2_6_genes1)
head(V2_6_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_6_gene_chunks <- split(V2_6_genes, ceiling(seq_along(V2_6_genes) / chunk_size))

# Create an empty list to store results
V2_6_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_6_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_6_gene_chunks), "\n")
  
  V2_6_results_list[[i]] <- gene2disease(
    gene     = V2_6_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_6_qresults_list <- lapply(V2_6_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_6_results_all <- do.call(rbind, V2_6_qresults_list)

# Remove duplicates (optional)
V2_6_results_all <- unique(V2_6_results_all)

# Second run 
V2_6_genes_tab <- V2_6_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_6_gene_chunks <- split(
  V2_6_genes_tab,
  ceiling(seq_along(V2_6_genes_tab) / chunk_size)
)

# Assign names like V2_6_genes1, V2_6_genes2, ...
for (i in seq_along(V2_6_gene_chunks)) {
  assign(paste0("V2_6_genes_tab", i), V2_6_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_6_gene_chunks)

# Example: view first few
head(V2_6_genes1)
head(V2_6_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_6_gene_chunks <- split(V2_6_genes_tab, ceiling(seq_along(V2_6_genes_tab) / chunk_size))

# Create an empty list to store results
V2_6_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_6_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_6_gene_chunks), "\n")
  
  V2_6_results_list_tab[[i]] <- gene2disease(
    gene     = V2_6_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_6_qresults_list <- lapply(V2_6_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_6_results_all <- do.call(rbind, V2_6_qresults_list)

# Remove duplicates (optional)
V2_6_results_all <- unique(V2_6_results_all)


df <- V2_6_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_6_11.3.25.xlsx")


#V2 Cliuster 7 #only 2 genes, no result
#########################################*
V2_7 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "7"))


# First run 
V2_7_genes <- V2_7$gene

length(V2_7_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_7_gene_chunks <- split(
  V2_7_genes,
  ceiling(seq_along(V2_7_genes) / chunk_size)
)

# Assign names like V2_7_genes1, V2_7_genes2, ...
for (i in seq_along(V2_7_gene_chunks)) {
  assign(paste0("V2_7_genes", i), V2_7_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_7_gene_chunks)

# Example: view first few
head(V2_7_genes1)
head(V2_7_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_7_gene_chunks <- split(V2_7_genes, ceiling(seq_along(V2_7_genes) / chunk_size))

# Create an empty list to store results
V2_7_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_7_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_7_gene_chunks), "\n")
  
  V2_7_results_list[[i]] <- gene2disease(
    gene     = V2_7_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


#V2 Cliuster 8
#########################################*
V2_8 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "8"))


# First run 
V2_8_genes <- V2_8$gene

length(V2_8_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_8_gene_chunks <- split(
  V2_8_genes,
  ceiling(seq_along(V2_8_genes) / chunk_size)
)

# Assign names like V2_8_genes1, V2_8_genes2, ...
for (i in seq_along(V2_8_gene_chunks)) {
  assign(paste0("V2_8_genes", i), V2_8_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_8_gene_chunks)

# Example: view first few
head(V2_8_genes1)
head(V2_8_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_8_gene_chunks <- split(V2_8_genes, ceiling(seq_along(V2_8_genes) / chunk_size))

# Create an empty list to store results
V2_8_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_8_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_8_gene_chunks), "\n")
  
  V2_8_results_list[[i]] <- gene2disease(
    gene     = V2_8_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_8_qresults_list <- lapply(V2_8_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_8_results_all <- do.call(rbind, V2_8_qresults_list)

# Remove duplicates (optional)
V2_8_results_all <- unique(V2_8_results_all)

# Second run 
V2_8_genes_tab <- V2_8_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_8_gene_chunks <- split(
  V2_8_genes_tab,
  ceiling(seq_along(V2_8_genes_tab) / chunk_size)
)

# Assign names like V2_8_genes1, V2_8_genes2, ...
for (i in seq_along(V2_8_gene_chunks)) {
  assign(paste0("V2_8_genes_tab", i), V2_8_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_8_gene_chunks)

# Example: view first few
head(V2_8_genes1)
head(V2_8_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_8_gene_chunks <- split(V2_8_genes_tab, ceiling(seq_along(V2_8_genes_tab) / chunk_size))

# Create an empty list to store results
V2_8_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_8_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_8_gene_chunks), "\n")
  
  V2_8_results_list_tab[[i]] <- gene2disease(
    gene     = V2_8_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_8_qresults_list <- lapply(V2_8_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_8_results_all <- do.call(rbind, V2_8_qresults_list)

# Remove duplicates (optional)
V2_8_results_all <- unique(V2_8_results_all)


df <- V2_8_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


library(openxlsx)

# Explicitly set how NAs are written
write.xlsx(df, file = "PsyGeNET_results_V2_8_11.3.25.xlsx", na.string = "NA")



#V2 Cliuster 9
#########################################*
V2_9 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "9"))


# First run 
V2_9_genes <- V2_9$gene

length(V2_9_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_9_gene_chunks <- split(
  V2_9_genes,
  ceiling(seq_along(V2_9_genes) / chunk_size)
)

# Assign names like V2_9_genes1, V2_9_genes2, ...
for (i in seq_along(V2_9_gene_chunks)) {
  assign(paste0("V2_9_genes", i), V2_9_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_9_gene_chunks)

# Example: view first few
head(V2_9_genes1)
head(V2_9_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_9_gene_chunks <- split(V2_9_genes, ceiling(seq_along(V2_9_genes) / chunk_size))

# Create an empty list to store results
V2_9_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_9_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_9_gene_chunks), "\n")
  
  V2_9_results_list[[i]] <- gene2disease(
    gene     = V2_9_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_9_qresults_list <- lapply(V2_9_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_9_results_all <- do.call(rbind, V2_9_qresults_list)

# Remove duplicates (optional)
V2_9_results_all <- unique(V2_9_results_all)

# Second run 
V2_9_genes_tab <- V2_9_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_9_gene_chunks <- split(
  V2_9_genes_tab,
  ceiling(seq_along(V2_9_genes_tab) / chunk_size)
)

# Assign names like V2_9_genes1, V2_9_genes2, ...
for (i in seq_along(V2_9_gene_chunks)) {
  assign(paste0("V2_9_genes_tab", i), V2_9_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_9_gene_chunks)

# Example: view first few
head(V2_9_genes1)
head(V2_9_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_9_gene_chunks <- split(V2_9_genes_tab, ceiling(seq_along(V2_9_genes_tab) / chunk_size))

# Create an empty list to store results
V2_9_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_9_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_9_gene_chunks), "\n")
  
  V2_9_results_list_tab[[i]] <- gene2disease(
    gene     = V2_9_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_9_qresults_list <- lapply(V2_9_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_9_results_all <- do.call(rbind, V2_9_qresults_list)

# Remove duplicates (optional)
V2_9_results_all <- unique(V2_9_results_all)


df <- V2_9_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_9_11.3.25.xlsx")



#V2 Cliuster 10
#########################################*
V2_10 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "10"))


# First run 
V2_10_genes <- V2_10$gene

length(V2_10_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_10_gene_chunks <- split(
  V2_10_genes,
  ceiling(seq_along(V2_10_genes) / chunk_size)
)

# Assign names like V2_10_genes1, V2_10_genes2, ...
for (i in seq_along(V2_10_gene_chunks)) {
  assign(paste0("V2_10_genes", i), V2_10_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_10_gene_chunks)

# Example: view first few
head(V2_10_genes1)
head(V2_10_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_10_gene_chunks <- split(V2_10_genes, ceiling(seq_along(V2_10_genes) / chunk_size))

# Create an empty list to store results
V2_10_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_10_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_10_gene_chunks), "\n")
  
  V2_10_results_list[[i]] <- gene2disease(
    gene     = V2_10_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract qresult only from valid S4 results, skip those that are characters
V2_10_qresults_list <- lapply(V2_10_results_list, function(x) {
  if (isS4(x) && "qresult" %in% slotNames(x)) {
    return(x@qresult)
  } else {
    return(NULL)
  }
})

# Drop NULL entries
V2_10_qresults_list <- Filter(Negate(is.null), V2_10_qresults_list)

# Combine all qresult data frames into one big data frame
V2_10_results_all <- do.call(rbind, V2_10_qresults_list)

# Remove duplicates (optional)
V2_10_results_all <- unique(V2_10_results_all)

# Second run 
V2_10_genes_tab <- V2_10_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_10_gene_chunks <- split(
  V2_10_genes_tab,
  ceiling(seq_along(V2_10_genes_tab) / chunk_size)
)

# Assign names like V2_10_genes1, V2_10_genes2, ...
for (i in seq_along(V2_10_gene_chunks)) {
  assign(paste0("V2_10_genes_tab", i), V2_10_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_10_gene_chunks)

# Example: view first few
head(V2_10_genes1)
head(V2_10_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_10_gene_chunks <- split(V2_10_genes_tab, ceiling(seq_along(V2_10_genes_tab) / chunk_size))

# Create an empty list to store results
V2_10_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_10_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_10_gene_chunks), "\n")
  
  V2_10_results_list_tab[[i]] <- gene2disease(
    gene     = V2_10_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_10_qresults_list <- lapply(V2_10_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_10_results_all <- do.call(rbind, V2_10_qresults_list)

# Remove duplicates (optional)
V2_10_results_all <- unique(V2_10_results_all)


df <- V2_10_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_10_11.3.25.xlsx")



#V2 Cliuster 11
#########################################*
V2_11 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "11"))


# First run 
V2_11_genes <- V2_11$gene

length(V2_11_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_11_gene_chunks <- split(
  V2_11_genes,
  ceiling(seq_along(V2_11_genes) / chunk_size)
)

# Assign names like V2_11_genes1, V2_11_genes2, ...
for (i in seq_along(V2_11_gene_chunks)) {
  assign(paste0("V2_11_genes", i), V2_11_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_11_gene_chunks)

# Example: view first few
head(V2_11_genes1)
head(V2_11_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_11_gene_chunks <- split(V2_11_genes, ceiling(seq_along(V2_11_genes) / chunk_size))

# Create an empty list to store results
V2_11_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_11_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_11_gene_chunks), "\n")
  
  V2_11_results_list[[i]] <- gene2disease(
    gene     = V2_11_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_11_qresults_list <- lapply(V2_11_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_11_results_all <- do.call(rbind, V2_11_qresults_list)

# Remove duplicates (optional)
V2_11_results_all <- unique(V2_11_results_all)

# Second run 
V2_11_genes_tab <- V2_11_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_11_gene_chunks <- split(
  V2_11_genes_tab,
  ceiling(seq_along(V2_11_genes_tab) / chunk_size)
)

# Assign names like V2_11_genes1, V2_11_genes2, ...
for (i in seq_along(V2_11_gene_chunks)) {
  assign(paste0("V2_11_genes_tab", i), V2_11_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_11_gene_chunks)

# Example: view first few
head(V2_11_genes1)
head(V2_11_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_11_gene_chunks <- split(V2_11_genes_tab, ceiling(seq_along(V2_11_genes_tab) / chunk_size))

# Create an empty list to store results
V2_11_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_11_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_11_gene_chunks), "\n")
  
  V2_11_results_list_tab[[i]] <- gene2disease(
    gene     = V2_11_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_11_qresults_list <- lapply(V2_11_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_11_results_all <- do.call(rbind, V2_11_qresults_list)

# Remove duplicates (optional)
V2_11_results_all <- unique(V2_11_results_all)


df <- V2_11_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_11_11.3.25.xlsx")



#V2 Cliuster 12
#########################################*
V2_12 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "12"))


# First run 
V2_12_genes <- V2_12$gene

length(V2_12_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_12_gene_chunks <- split(
  V2_12_genes,
  ceiling(seq_along(V2_12_genes) / chunk_size)
)

# Assign names like V2_12_genes1, V2_12_genes2, ...
for (i in seq_along(V2_12_gene_chunks)) {
  assign(paste0("V2_12_genes", i), V2_12_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_12_gene_chunks)

# Example: view first few
head(V2_12_genes1)
head(V2_12_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_12_gene_chunks <- split(V2_12_genes, ceiling(seq_along(V2_12_genes) / chunk_size))

# Create an empty list to store results
V2_12_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_12_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_12_gene_chunks), "\n")
  
  V2_12_results_list[[i]] <- gene2disease(
    gene     = V2_12_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_12_qresults_list <- lapply(V2_12_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_12_results_all <- do.call(rbind, V2_12_qresults_list)

# Remove duplicates (optional)
V2_12_results_all <- unique(V2_12_results_all)

# Second run 
V2_12_genes_tab <- V2_12_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_12_gene_chunks <- split(
  V2_12_genes_tab,
  ceiling(seq_along(V2_12_genes_tab) / chunk_size)
)

# Assign names like V2_12_genes1, V2_12_genes2, ...
for (i in seq_along(V2_12_gene_chunks)) {
  assign(paste0("V2_12_genes_tab", i), V2_12_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_12_gene_chunks)

# Example: view first few
head(V2_12_genes1)
head(V2_12_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_12_gene_chunks <- split(V2_12_genes_tab, ceiling(seq_along(V2_12_genes_tab) / chunk_size))

# Create an empty list to store results
V2_12_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_12_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_12_gene_chunks), "\n")
  
  V2_12_results_list_tab[[i]] <- gene2disease(
    gene     = V2_12_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_12_qresults_list <- lapply(V2_12_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_12_results_all <- do.call(rbind, V2_12_qresults_list)

# Remove duplicates (optional)
V2_12_results_all <- unique(V2_12_results_all)


df <- V2_12_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_12_11.3.25.xlsx")



#V2 Cliuster 13
#########################################*
V2_13 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "13"))


# First run 
V2_13_genes <- V2_13$gene

length(V2_13_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_13_gene_chunks <- split(
  V2_13_genes,
  ceiling(seq_along(V2_13_genes) / chunk_size)
)

# Assign names like V2_13_genes1, V2_13_genes2, ...
for (i in seq_along(V2_13_gene_chunks)) {
  assign(paste0("V2_13_genes", i), V2_13_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_13_gene_chunks)

# Example: view first few
head(V2_13_genes1)
head(V2_13_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_13_gene_chunks <- split(V2_13_genes, ceiling(seq_along(V2_13_genes) / chunk_size))

# Create an empty list to store results
V2_13_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_13_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_13_gene_chunks), "\n")
  
  V2_13_results_list[[i]] <- gene2disease(
    gene     = V2_13_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_13_qresults_list <- lapply(V2_13_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_13_results_all <- do.call(rbind, V2_13_qresults_list)

# Remove duplicates (optional)
V2_13_results_all <- unique(V2_13_results_all)

# Second run 
V2_13_genes_tab <- V2_13_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_13_gene_chunks <- split(
  V2_13_genes_tab,
  ceiling(seq_along(V2_13_genes_tab) / chunk_size)
)

# Assign names like V2_13_genes1, V2_13_genes2, ...
for (i in seq_along(V2_13_gene_chunks)) {
  assign(paste0("V2_13_genes_tab", i), V2_13_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_13_gene_chunks)

# Example: view first few
head(V2_13_genes1)
head(V2_13_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_13_gene_chunks <- split(V2_13_genes_tab, ceiling(seq_along(V2_13_genes_tab) / chunk_size))

# Create an empty list to store results
V2_13_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_13_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_13_gene_chunks), "\n")
  
  V2_13_results_list_tab[[i]] <- gene2disease(
    gene     = V2_13_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_13_qresults_list <- lapply(V2_13_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_13_results_all <- do.call(rbind, V2_13_qresults_list)

# Remove duplicates (optional)
V2_13_results_all <- unique(V2_13_results_all)


df <- V2_13_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_13_11.3.25.xlsx")



#V2 Cliuster 14
#########################################*
V2_14 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "14"))


# First run 
V2_14_genes <- V2_14$gene

length(V2_14_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_14_gene_chunks <- split(
  V2_14_genes,
  ceiling(seq_along(V2_14_genes) / chunk_size)
)

# Assign names like V2_14_genes1, V2_14_genes2, ...
for (i in seq_along(V2_14_gene_chunks)) {
  assign(paste0("V2_14_genes", i), V2_14_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_14_gene_chunks)

# Example: view first few
head(V2_14_genes1)
head(V2_14_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_14_gene_chunks <- split(V2_14_genes, ceiling(seq_along(V2_14_genes) / chunk_size))

# Create an empty list to store results
V2_14_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_14_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_14_gene_chunks), "\n")
  
  V2_14_results_list[[i]] <- gene2disease(
    gene     = V2_14_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


#query had 0 pages


#V2 Cliuster 15
#########################################*
V2_15 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "15"))


# First run 
V2_15_genes <- V2_15$gene

length(V2_15_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_15_gene_chunks <- split(
  V2_15_genes,
  ceiling(seq_along(V2_15_genes) / chunk_size)
)

# Assign names like V2_15_genes1, V2_15_genes2, ...
for (i in seq_along(V2_15_gene_chunks)) {
  assign(paste0("V2_15_genes", i), V2_15_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_15_gene_chunks)

# Example: view first few
head(V2_15_genes1)
head(V2_15_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_15_gene_chunks <- split(V2_15_genes, ceiling(seq_along(V2_15_genes) / chunk_size))

# Create an empty list to store results
V2_15_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_15_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_15_gene_chunks), "\n")
  
  V2_15_results_list[[i]] <- gene2disease(
    gene     = V2_15_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_15_qresults_list <- lapply(V2_15_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_15_results_all <- do.call(rbind, V2_15_qresults_list)

# Remove duplicates (optional)
V2_15_results_all <- unique(V2_15_results_all)

# Second run 
V2_15_genes_tab <- V2_15_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_15_gene_chunks <- split(
  V2_15_genes_tab,
  ceiling(seq_along(V2_15_genes_tab) / chunk_size)
)

# Assign names like V2_15_genes1, V2_15_genes2, ...
for (i in seq_along(V2_15_gene_chunks)) {
  assign(paste0("V2_15_genes_tab", i), V2_15_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_15_gene_chunks)

# Example: view first few
head(V2_15_genes1)
head(V2_15_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_15_gene_chunks <- split(V2_15_genes_tab, ceiling(seq_along(V2_15_genes_tab) / chunk_size))

# Create an empty list to store results
V2_15_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_15_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_15_gene_chunks), "\n")
  
  V2_15_results_list_tab[[i]] <- gene2disease(
    gene     = V2_15_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_15_qresults_list <- lapply(V2_15_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_15_results_all <- do.call(rbind, V2_15_qresults_list)

# Remove duplicates (optional)
V2_15_results_all <- unique(V2_15_results_all)


df <- V2_15_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_15_11.3.25.xlsx")
library(openxlsx)
write.xlsx(df, file = "PsyGeNET_results_V2_15_11.3.25.xlsx", na.string = "NA")



#V2 Cliuster 16
#########################################*
V2_16 <- as.data.frame(readxl::read_excel(
  "/Users/mukamiwamalwa/Library/CloudStorage/GoogleDrive-alexandra.wamalwa@nyspi.columbia.edu/Shared drives/Boldrini Lab/Mukami Wamalwa/Disgenet/visiums_11.3/v2_second_attempt_mdd_vs_control_degs_with_ncbi_uniprot_11_3_2025.xlsx",
  sheet = "16"))


# First run 
V2_16_genes <- V2_16$gene

length(V2_16_genes)

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_16_gene_chunks <- split(
  V2_16_genes,
  ceiling(seq_along(V2_16_genes) / chunk_size)
)

# Assign names like V2_16_genes1, V2_16_genes2, ...
for (i in seq_along(V2_16_gene_chunks)) {
  assign(paste0("V2_16_genes", i), V2_16_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_16_gene_chunks)

# Example: view first few
head(V2_16_genes1)
head(V2_16_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_16_gene_chunks <- split(V2_16_genes, ceiling(seq_along(V2_16_genes) / chunk_size))

# Create an empty list to store results
V2_16_results_list <- list()

# Loop over each chunk
for (i in seq_along(V2_16_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_16_gene_chunks), "\n")
  
  V2_16_results_list[[i]] <- gene2disease(
    gene     = V2_16_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}

# Extract the qresult slot (the actual table) from each S4 result
V2_16_qresults_list <- lapply(V2_16_results_list, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_16_results_all <- do.call(rbind, V2_16_qresults_list)

# Remove duplicates (optional)
V2_16_results_all <- unique(V2_16_results_all)

# Second run 
V2_16_genes_tab <- V2_16_results_all$gene_symbol

# Limited to 100 genes are one time, split CA1_genes into 2 with 100 and 45 entries
# Number of genes per chunk
chunk_size <- 100

# Split into groups of 100
V2_16_gene_chunks <- split(
  V2_16_genes_tab,
  ceiling(seq_along(V2_16_genes_tab) / chunk_size)
)

# Assign names like V2_16_genes1, V2_16_genes2, ...
for (i in seq_along(V2_16_gene_chunks)) {
  assign(paste0("V2_16_genes_tab", i), V2_16_gene_chunks[[i]])
}

# Check how many chunks were created
length(V2_16_gene_chunks)

# Example: view first few
head(V2_16_genes1)
head(V2_16_genes2)

# Ensure you have all gene chunks
chunk_size <- 100
V2_16_gene_chunks <- split(V2_16_genes_tab, ceiling(seq_along(V2_16_genes_tab) / chunk_size))

# Create an empty list to store results
V2_16_results_list_tab <- list()

# Loop over each chunk
for (i in seq_along(V2_16_gene_chunks)) {
  cat("Processing chunk", i, "of", length(V2_16_gene_chunks), "\n")
  
  V2_16_results_list_tab[[i]] <- gene2disease(
    gene     = V2_16_gene_chunks[[i]],
    database = "PSYGENET",
    verbose  = TRUE
  )
  
  # Optional short delay to avoid overwhelming the API
  Sys.sleep(2)
}


# Extract the qresult slot (the actual table) from each S4 result
V2_16_qresults_list <- lapply(V2_16_results_list_tab, function(x) x@qresult)

# Combine all qresult data frames into one big data frame
V2_16_results_all <- do.call(rbind, V2_16_qresults_list)

# Remove duplicates (optional)
V2_16_results_all <- unique(V2_16_results_all)


df <- V2_16_results_all

# Optional: filter by score or number of supporting PMIDs
df_filtered <- df %>%
  filter(score >= 0.8)  # adjust threshold if needed

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

ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(color = "grey60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 6) +
  scale_color_manual(values = c(Gene = "lightblue", Disease = "hotpink")) +
  theme_void() +
  theme(legend.position = "none")  # hide legend if desired


write.xlsx(df, file = "PsyGeNET_results_V2_16_11.3.25.xlsx")

