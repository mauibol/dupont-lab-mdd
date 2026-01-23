#################
#Annotate Gene Lists for Visium v1 and v2 Differential  Enrichment Tables with NCBI and UNIPROT
# Requires v1 and v2 differential_enrichment_batch  .xlsx sheet from visium_differential_gene_enrichment.R
#By Victor Anosike
#Adapted from code by Lucia Polizzi
#################

# Important Libraries

start_time <- proc.time()

pacman::p_load('mgsub', 'gsubfn', 'readxl', 'openxlsx', 'gridExtra', 'cowplot', 'tidyr',
               'dplyr', 'tidyverse', 'reshape2', 'rstatix', 'HybridMTest', 'ggrepel',
               'plyr', 'foreach', 'doParallel', 'plotly', 'pheatmap', 'gplots',
               'FactoMineR', 'factoextra', 'limma', 'leukemiasEset', 'svglite', 'patchwork',
               'ggplot2', 'UpSetR', 'VennDiagram', 'stringr', 'ggpubr', 'pdftools',
               'httr', 'readr', 'tibble', 'data.table', 'ggvenn', 'ComplexUpset',
               'eulerr', 'circlize', 'ComplexHeatmap', 'purrr', 'arrow', 'Rfit', 'scales',
               'rentrez', 'progress', 'jsonlite', 'UniProt.ws', 'rvest', 'readxl')


############################# Visium v1 ##############################################################################
######################################################################################################################


# Get the names of all the sheets in the gene_enrichment excel sheet generated from visium_differential_gene_enrichment.R
sheet_names <- unlist(as.list(excel_sheets('/path/to/v1_differential_enrichment_date_batch_effect.xlsx')))


#Create List That will hold all Cluster_Sheet S4 Objects, ie all the sheets within the loaded workbook
book_df <- list()

#Create S4 object Cluster_Sheet that will contain both the sheet name and a seperate dataframe that contains all the information from that sheet
setClass("Cluster_Sheet", 
         slots = c(sheet_name = "character",
                   sheet_data = "data.frame",
                   ncbi_summary = "data.frame",
                   uniprot_summary = "data.frame"))

#Loop through sheet_names to add each sheet from the excel file as a Cluster_Sheet S4 object
for (i in sheet_names) {
  
  data <- read_excel(
    "/path/to/v1_differential_enrichment_date_batch_effect.xlsx",
    sheet = i
  )
  
  sheet_df <- new("Cluster_Sheet",
                  sheet_name = i,
                  sheet_data = data,
                  ncbi_summary = data.frame(),
                  uniprot_summary = data.frame())
  
  book_df <- append(book_df,sheet_df)
}


#################################
# Add NCBI summaries via gene_info
#################################

# ---- STEP 1: Download Homo sapiens gene_info for mapping ----
url_info  <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
dest_info <- "/path/to/save/Homo_sapiens.gene_info.gz"

if (!file.exists(dest_info)) {
  download.file(url_info, dest_info, mode = "wb")
}

gene_info <- read.delim(
  gzfile(dest_info),
  header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE
)

# ---- STEP 2: Map Symbol -> GeneID ----
gene_map <- gene_info[, c("Symbol", "GeneID")]


for (i in 1:16) {
  
  # ---- STEP 3: Prepare GeneIDs to query ----
  ids <- na.omit(gene_map$GeneID[gene_map$Symbol %in% book_df[[i]]@sheet_data$gene])
  ids <- unique(ids)
  
  # ---- STEP 4: Function to fetch summaries in batch ----
  fetch_summaries_batch <- function(id_vec) {
    summs <- entrez_summary(db = "gene", id = id_vec)
    tibble(
      GeneID = names(summs),
      NCBI_Summary = vapply(summs, function(x) {
        if (is.null(x$summary)) NA_character_ else x$summary
      }, character(1))
    )
  }
  
  # ---- STEP 5: Split IDs into chunks of 200 and fetch ----
  id_chunks <- split(ids, ceiling(seq_along(ids)/200))
  summaries_list <- lapply(id_chunks, fetch_summaries_batch)
  gene_summaries <- bind_rows(summaries_list)
  gene_map <- gene_map %>% mutate(GeneID = as.character(GeneID))
  gene_summaries <- gene_summaries %>% mutate(GeneID = as.character(GeneID))
  gene_map_clean <- gene_info %>%
    dplyr::filter(Synonyms == "-" | Symbol %in% book_df[[i]]@sheet_data$gene) %>%
    dplyr::select(Symbol, GeneID) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE)
  
  # ---- STEP 6: Merge summaries back into book_df[[i]]@sheet_data ----
  book_df[[i]]@ncbi_summary <- book_df[[i]]@sheet_data %>%
    left_join(gene_map_clean %>% mutate(GeneID = as.character(GeneID)),
              by = c("gene" = "Symbol")) %>%
    left_join(gene_summaries %>% mutate(GeneID = as.character(GeneID)),
              by = "GeneID")
  
}


#################################
# Add UniProt functions via UniProt REST API
#################################

url <- "https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606%20AND%20reviewed:true)&format=tsv&fields=accession,gene_primary,cc_function"
dest <- "/path/to/save/uniprot_human_function.tsv"
if (!file.exists(dest)) { download.file(url, dest, mode = "wb")}
uniprot_func <- read.delim(dest, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)

for (i in 1:16) {
  
  # 1. Build a gene → function mapping
  uniprot_map <- uniprot_func %>%
    dplyr::select(gene = 'Gene.Names..primary.', Uniprot_Function = 'Function..CC.')
  uniprot_map_clean <- uniprot_map %>%
    group_by(gene) %>%
    summarise(Uniprot_Function = first(Uniprot_Function), .groups = "drop") %>%
    mutate(Uniprot_Function = sub("^FUNCTION: ", "", Uniprot_Function))
  
  # 2. Merge into book_df[[i]]@sheet_data_clean
  book_df[[i]]@uniprot_summary <- book_df[[i]]@ncbi_summary %>%
    left_join(uniprot_map_clean, by = "gene")
  
  
}



#################################
# Save data
#################################

#3) Save as xlsx

wb <- createWorkbook()
for (i in seq_along(book_df)) {
  addWorksheet(wb, book_df[[i]]@sheet_name)
  writeData(wb, sheet = book_df[[i]]@sheet_name, x = as.data.frame(book_df[[i]]@uniprot_summary))
}

saveWorkbook(wb, "/path/to/v1_gene_differential_enrichment_batch_effect_with_ncbi_uniprot_date.xlsx", overwrite = TRUE)

#################################
# Record the end time
#################################
end_time <- proc.time()
time_taken <- end_time - start_time
print(time_taken)

############################# Visium v2 ##############################################################################
######################################################################################################################

# Get the names of all the sheets in the gene_enrichment excel sheet generated from visium_differential_gene_enrichment.R
sheet_names <- unlist(as.list(excel_sheets("/path/to/v2_differential_enrichment_date_batch_effect.xlsx")))


#Create List That will hold all Cluster_Sheet S4 Objects, ie all the sheets within the loaded workbook
book_df <- list()

#Create S4 object Cluster_Sheet that will contain both the sheet name and a seperate dataframe that contains all the information from that sheet
setClass("Cluster_Sheet", 
         slots = c(sheet_name = "character",
                   sheet_data = "data.frame",
                   ncbi_summary = "data.frame",
                   uniprot_summary = "data.frame"))

#Loop through sheet_names to add each sheet from the excel file as a Cluster_Sheet S4 object
for (i in sheet_names) {
  
  data <- read_excel(
    "/path/to/v2_differential_enrichment_date_batch_effect.xlsx",
    sheet = i
  )
  
  sheet_df <- new("Cluster_Sheet",
                  sheet_name = i,
                  sheet_data = data,
                  ncbi_summary = data.frame(),
                  uniprot_summary = data.frame())
  
  book_df <- append(book_df,sheet_df)
}


#################################
# Add NCBI summaries via gene_info
#################################

# ---- STEP 1: Download Homo sapiens gene_info for mapping ----
url_info  <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
dest_info <- "/path/to/save/Homo_sapiens.gene_info.gz"

if (!file.exists(dest_info)) {
  download.file(url_info, dest_info, mode = "wb")
}

gene_info <- read.delim(
  gzfile(dest_info),
  header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE
)

# ---- STEP 2: Map Symbol -> GeneID ----
gene_map <- gene_info[, c("Symbol", "GeneID")]


for (i in 1:19) {
  
  # ---- STEP 3: Prepare GeneIDs to query ----
  ids <- na.omit(gene_map$GeneID[gene_map$Symbol %in% book_df[[i]]@sheet_data$gene])
  ids <- unique(ids)
  
  # ---- STEP 4: Function to fetch summaries in batch ----
  fetch_summaries_batch <- function(id_vec) {
    summs <- entrez_summary(db = "gene", id = id_vec)
    tibble(
      GeneID = names(summs),
      NCBI_Summary = vapply(summs, function(x) {
        if (is.null(x$summary)) NA_character_ else x$summary
      }, character(1))
    )
  }
  
  # ---- STEP 5: Split IDs into chunks of 200 and fetch ----
  id_chunks <- split(ids, ceiling(seq_along(ids)/200))
  summaries_list <- lapply(id_chunks, fetch_summaries_batch)
  gene_summaries <- bind_rows(summaries_list)
  gene_map <- gene_map %>% mutate(GeneID = as.character(GeneID))
  gene_summaries <- gene_summaries %>% mutate(GeneID = as.character(GeneID))
  gene_map_clean <- gene_info %>%
    dplyr::filter(Synonyms == "-" | Symbol %in% book_df[[i]]@sheet_data$gene) %>%
    dplyr::select(Symbol, GeneID) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE)
  
  # ---- STEP 6: Merge summaries back into book_df[[i]]@sheet_data ----
  book_df[[i]]@ncbi_summary <- book_df[[i]]@sheet_data %>%
    left_join(gene_map_clean %>% mutate(GeneID = as.character(GeneID)),
              by = c("gene" = "Symbol")) %>%
    left_join(gene_summaries %>% mutate(GeneID = as.character(GeneID)),
              by = "GeneID")
  
}



#################################
# Add UniProt functions via UniProt REST API
#################################

url <- "https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606%20AND%20reviewed:true)&format=tsv&fields=accession,gene_primary,cc_function"
dest <- "/path/to/save/uniprot_human_function.tsv"
if (!file.exists(dest)) { download.file(url, dest, mode = "wb")}
uniprot_func <- read.delim(dest, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)

for (i in 1:19) {
  
  # 1. Build a gene → function mapping
  uniprot_map <- uniprot_func %>%
    dplyr::select(gene = 'Gene.Names..primary.', Uniprot_Function = 'Function..CC.')
  uniprot_map_clean <- uniprot_map %>%
    group_by(gene) %>%
    summarise(Uniprot_Function = first(Uniprot_Function), .groups = "drop") %>%
    mutate(Uniprot_Function = sub("^FUNCTION: ", "", Uniprot_Function))
  
  # 2. Merge into book_df[[i]]@sheet_data_clean
  book_df[[i]]@uniprot_summary <- book_df[[i]]@ncbi_summary %>%
    left_join(uniprot_map_clean, by = "gene")
  
  
}



#################################
# Save data
#################################

#3) Save as xlsx

wb <- createWorkbook()
for (i in seq_along(book_df)) {
  addWorksheet(wb, book_df[[i]]@sheet_name)
  writeData(wb, sheet = book_df[[i]]@sheet_name, x = as.data.frame(book_df[[i]]@uniprot_summary))
}

saveWorkbook(wb, "/path/to/v2_gene_differential_enrichment_batch_effect_with_ncbi_uniprot_date.xlsx", overwrite = TRUE)

#################################
# Record the end time
#################################
end_time <- proc.time()
time_taken <- end_time - start_time
print(time_taken)




