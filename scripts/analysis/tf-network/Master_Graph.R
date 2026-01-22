library(tidyverse)
library(dplyr)
library(tidygraph)
library(stringr)

.build_master_graph <- function(df) {
  
  # --- 1) Make a 'motif' column from TFBS_name --------------------------------
  df <- df %>%
    mutate(
      motif = sub(".*(MA[^_\\s]+)$", "\\1", TFBS_name)
    )
  
  # --- 2) Create a motif index per TF -----------------------------------------
  motif_index <- df %>%
    distinct(tf_name, motif) %>%
    group_by(tf_name) %>%
    arrange(motif, .by_group = TRUE) %>%
    mutate(motif_num = row_number()) %>%
    ungroup()
  
  # --- 3) Join back; build tf_detail id ---------------------------------------
  df <- df %>%
    left_join(motif_index, by = c("tf_name","motif")) %>%
    mutate(
      tf_detail = paste0(tf_name, "@", TFBS_start, "#", motif_num)
    )
  
  # --- 4) EDGES ---------------------------------------------------------------
  
  # 4a) tf_detail -> TF
  edges_detail_tf <- df %>%
    transmute(
      from = tf_detail,
      to   = TFBS_name
    ) %>%
    distinct()
  
  # 4b) TF -> gene (minimal: connect both ends, no aggregation)
  edges_tf_gene <- df %>%
    group_by(TFBS_name, gene) %>%
    summarise(
      from = first(TFBS_name),
      to   = first(gene),
      mean_log2fc = mean(footprint_log2fc, na.rm = TRUE),
      .groups = "drop" 
    ) %>%
    distinct()
  
  # 4c) gene -> protein
  edges_gene_protein <- df %>%
    filter(!is.na(Protein) & Protein != "") %>%
    transmute(
      from = gene,
      to   = Protein,
      direction = direction_match
    ) %>%
    distinct()
  
  # 4d) all edges
  edges_all <- bind_rows(edges_detail_tf, edges_tf_gene, edges_gene_protein)
  
  # --- 5) NODES ---------------------------------------------------------------
  
  # 5a) tf_detail nodes (site-level)
  tf_detail_nodes <- df %>%
    transmute(
      name         = tf_detail,
      node_type    = "tf_detail",
      # label available for later use; keep it compact
      label        = paste0(tf_name, " #", motif_num),
      tf_name      = tf_name,
      tfbs_score   = coalesce(TFBS_score, NA_real_),
      tfbs_start   = TFBS_start,
      tfbs_end     = TFBS_end,
      tfbs_chr     = TFBS_chr,
      motif        = motif,
      motif_num    = motif_num,
      binding_size = as.numeric(TFBS_end) - as.numeric(TFBS_start)
    ) %>%
    distinct(name, .keep_all = TRUE)
  
  # 5b) TF parent nodes (one per TFBS_name)
  tf_parent_nodes <- df %>%
    distinct(TFBS_name, tf_name) %>%
    transmute(
      name      = TFBS_name,
      node_type = "tf",
      label     = tf_name
    )
  
  # 5c) gene nodes (collapse duplicates; carry gene-level metrics once)
  gene_nodes <- df %>%
    distinct(gene, .keep_all = TRUE) %>%
    transmute(
      name      = gene,
      node_type = "gene",
      label     = gene,
      logfc     = coalesce(logFC, NA_real_),
      avgexp    = coalesce(AveExpr, NA_real_)
    )
  
  # 5d) protein nodes (significance: boolean + labeled factor)
  protein_nodes <- df %>%
    filter(!is.na(Protein) & Protein != "") %>%
    distinct(Protein, .keep_all = TRUE) %>%
    transmute(
      name        = Protein,
      node_type   = "protein",
      label       = Protein,
      p_value     = suppressWarnings(as.numeric(p_value_protein)),
      protein_sig = !is.na(p_value_protein) & p_value_protein < 0.05,
      protein_sig_label = case_when(
        is.na(p_value_protein) ~ "NA",
        p_value_protein < 0.05 ~ "p<0.05",
        TRUE                   ~ "p≥0.05"
      )
    ) %>%
    mutate(
      protein_sig_label = factor(protein_sig_label,
                                 levels = c("p<0.05","p≥0.05","NA"))
    )
  
  # 5e) Combine nodes; keep unique names and drop NAs
  nodes_all <- bind_rows(tf_detail_nodes, tf_parent_nodes, gene_nodes, protein_nodes) %>%
    filter(!is.na(name) & name != "") %>%
    distinct(name, .keep_all = TRUE)
  
  # --- 6) Ensure edges have valid endpoints -----------------------------------
  edges_all <- edges_all %>%
    inner_join(nodes_all %>% select(name) %>% rename(from = name), by = "from") %>%
    inner_join(nodes_all %>% select(name) %>% rename(to   = name), by = "to")
  
  # --- 7) Build graph (directed) ---------------------------------------------
  g <- tidygraph::tbl_graph(nodes = nodes_all, edges = edges_all, directed = TRUE)
  
  # --- 8) Return --------------------------------------------------------------
  list(
    nodes = nodes_all,
    edges = edges_all,
    g     = g
  )
}


summary(.build_master_graph(ExN.1_network)$edges)
