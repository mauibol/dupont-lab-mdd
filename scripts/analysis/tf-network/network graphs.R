source("Master_Graph.R")

library(tidyverse)
library(tidygraph)
library(dplyr)
library(ggraph)
library(grid)

Full_graph <- function(df, sig_protein_only = TRUE) {
  set.seed(42)
  
  g <- .build_master_graph(df)$g
  
  # ---- keep only TF / gene / protein; drop tf_detail and its edges -----------
  g <- g %>%
    activate(nodes) %>%
    mutate(
      type        = tolower(node_type),
      degree      = centrality_degree(),
      node_size   = 1,
      sig_protein = if ("sig_protein" %in% names(cur_data_all())) {
        as.logical(sig_protein)
      } else if ("protein_sig" %in% names(cur_data_all())) {
        as.logical(protein_sig)
      } else if ("p_value" %in% names(cur_data_all())) {
        !is.na(p_value) & p_value < 0.05
      } else {NA}) %>%
    filter(type %in% c("tf","gene","protein")) %>%
    activate(edges) %>%
    mutate(
      from_type = .N()$type[from],
      to_type   = .N()$type[to]
    ) %>%
    filter(from_type %in% c("tf","gene","protein") & to_type %in% c("tf","gene","protein"))
  
  # ---- if sig-only: keep significant proteins and their upstreams ------------
  if (isTRUE(sig_protein_only)) {
    nd <- as_tibble(g, active = "nodes")
    ed <- as_tibble(g, active = "edges")
    
    prot_idx <- which(nd$type == "protein" & nd$sig_protein %in% TRUE)
    if (length(prot_idx) == 0L) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No significant proteins (p < 0.05)", size = 5) +
          theme_void()
      )
    }
    
    # genes upstream of those proteins (gene -> protein)
    gene_idx <- unique(ed$from[ed$to %in% prot_idx & nd$type[ed$from] == "gene"])
    
    # TFs upstream of those genes (tf -> gene)
    tf_idx <- unique(ed$from[ed$to %in% gene_idx & nd$type[ed$from] == "tf"])
    
    keep_idx   <- sort(unique(c(prot_idx, gene_idx, tf_idx)))
    keep_names <- nd$name[keep_idx]
    
    # In tidygraph, filtering nodes keeps only incident edges automatically
    g <- g %>%
      activate(nodes) %>%
      filter(name %in% keep_names) %>%
      activate(edges)
  }
  
  n_nodes <- g %>% as_tibble(active = "nodes") %>% nrow()
  n_edges <- g %>% as_tibble(active = "edges") %>% nrow()
  
  # ---- plot ------------------------------------------------------------------
  ggraph(g, layout = "stress") +
    geom_edge_link(linetype = "solid", colour = "grey60", width = 0.5, alpha = 0.7,
                   arrow = arrow(length = grid::unit(1, "mm")),
                   start_cap = circle(3.5, "mm"), end_cap = circle(3.5, "mm")) +
    
    # 1) TF (parent) first (no fill scale used here)
    geom_node_point(
      data = ~ dplyr::filter(., type == "tf"),
      aes(size = node_size),
      shape = 22, fill = NA, color = "black", stroke = 0.5, show.legend = FALSE
    ) +
    
    # 2) PROTEINS (discrete fill = sig_protein) + its OWN scale
    geom_node_point(
      data = ~ dplyr::filter(., type == "protein"),
      aes(size = node_size, fill = sig_protein),
      shape = 24, color = "black", stroke = 0.25
    ) +
    scale_fill_manual(
      name = if (isTRUE(sig_protein_only)) "Protein (p<0.05 only)" else "Protein p<0.05",
      values = c(`TRUE` = "#4DAF4A", `FALSE` = "#FFFFFF"),
      na.value = "#BDBDBD"
    ) +
    
    # >>> reset fill scale so next layer can use a continuous gradient <<<
    ggnewscale::new_scale_fill() +
    
    # 3) GENES (continuous fill = logfc) + its OWN scale
    geom_node_point(
      data = ~ dplyr::filter(., type == "gene"),
      aes(size = node_size, fill = logfc),
      shape = 21, color = "black", stroke = 0.25
    ) +
    scale_fill_gradient2(
      name = "Gene log2FC",
      low = "blue", mid = "white", high = "red",
      midpoint = 0, 
      limits = c(-1.5, 1.5),na.value = "#E0E0E0"
    ) +
    
    ggrepel::geom_text_repel(
      data = ~ dplyr::filter(., !is.na(name)),
      aes(x = x, y = y, label = name),
      size = 3, max.overlaps = Inf, segment.size = 0.1,
      box.padding = 0.3, point.padding = 0.2, min.segment.length = 0, force = 2,
      bg.color = "white", bg.r = 0.15
    ) +
    
    guides(size = "none") +
    labs(title = paste0("Full Network (N=", n_nodes, ", E=", n_edges, ")")) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 12))
}

.match_tf_nodes <- function(g, TF) {
  # Normalize TF input into a clean character vector
  TF <- unique(TF[nzchar(TF)])              # drop empties, dedupe
  TF_lc <- tolower(TF)
  
  # Work on TF nodes only
  nd_all <- g %>% activate(nodes) %>% as_tibble()
  is_tf  <- tolower(nd_all$node_type) == "tf"
  nd_tf  <- nd_all[is_tf, , drop = FALSE]
  
  # Columns we’ll try to match against, in this order
  cand_cols <- c("name", "label", "tf_name", "output_prefix")
  cand_cols <- cand_cols[cand_cols %in% names(nd_tf)]  # keep those that exist
  
  # Lowercase lookup vectors for each candidate column
  lookup <- lapply(cand_cols, function(col) tolower(as.character(nd_tf[[col]])))
  names(lookup) <- cand_cols
  
  # For each requested TF token, mark TF rows that match in ANY candidate column
  match_any <- rep(FALSE, nrow(nd_tf))
  for (t in TF_lc) {
    hit_one <- Reduce(`|`, lapply(lookup, function(v) ifelse(is.na(v), FALSE, v == t)))
    match_any <- match_any | hit_one
  }
  
  # Indices back to the full node table
  idx <- which(is_tf)[match_any]
  
  # Pretty display names for your title: prefer 'label' if present, else 'name'
  display <- if ("label" %in% names(nd_tf)) nd_tf$label[match_any] else nd_tf$name[match_any]
  display <- unique(display)
  
  # Optional: warn about tokens that didn’t match anything
  if (length(idx) == 0L) {
    stop(sprintf(
      "Requested TFs not found (tried columns: %s): %s",
      paste(cand_cols, collapse = ", "),
      paste(TF, collapse = ", ")
    ))
  } else {
    # find non-matching tokens for debugging
    found_vals <- unique(unlist(lookup))
    miss <- TF_lc[!TF_lc %in% found_vals]
    if (length(miss)) message("No TF match for: ", paste(miss, collapse = ", "))
  }
  
  list(idx = idx, display = display)
}

Partial_graph <- function(df, TF = "", cell = "", sig_protein_only = TRUE) {
  set.seed(42)
  
  # 1) Build graph (must contain mean_log2fc on TF->gene edges, logfc on genes)
  g <- .build_master_graph(df)$g
  
  # 2) Node attributes (safe & consistent)
  g <- g %>%
    activate(nodes) %>%
    mutate(
      type      = tolower(node_type),
      degree    = centrality_degree(),
      node_size = case_when(
        type == "tf"      ~ 4,
        type == "gene"    ~ 4,
        type == "protein" ~ 4,
        TRUE              ~ 4
      )
    )
  
  # Derive/normalize sig_protein
  
  if (!"logfc" %in% names(as_tibble(g, active = "nodes"))) {
    g <- g %>% activate(nodes) %>% mutate(logfc = NA_real_)
  }
  
  # Keep only TF/Gene/Protein nodes
  g <- g %>% activate(nodes) %>% filter(type %in% c("tf","gene","protein"))
  
  # --- Multi-TF radius (≤2 hops) using tidygraph only --------------------------
  nd_tbl <- g %>% activate(nodes) %>% as_tibble()
  
  m <- .match_tf_nodes(g, TF)
  tf_idx <- m$idx
  if (length(tf_idx) == 0L) stop(sprintf("Requested TFs not found: %s", paste(TF, collapse = ", ")))
  
  # Compute min distance from any TF seed using tidygraph::node_distance_from
  dist_min <- rep(Inf, nrow(nd_tbl))
  for (i in tf_idx) {
    d <- g %>%
      activate(nodes) %>%
      mutate(d = node_distance_from(node = i, mode = "out")) %>%
      pull(d)
    # take elementwise min, ignoring NA
    dist_min <- pmin(dist_min, ifelse(is.na(d), Inf, d))
  }
  
  g <- g %>%
    activate(nodes) %>%
    mutate(dist_from_tf = dist_min) %>%
    filter(is.finite(dist_from_tf) & dist_from_tf <= 2) %>%
    activate(edges) %>%
    mutate(
      from_type = .N()$type[from],
      to_type   = .N()$type[to]
    ) %>%
    filter(from_type %in% c("tf","gene","protein") & to_type %in% c("tf","gene","protein"))
  
  # --- STRICT subgraph: only significant proteins & their upstreams ------------
  nd <- as_tibble(g, active = "nodes")
  ed <- as_tibble(g, active = "edges")
  
  prot_idx <- which(nd$type == "protein" & nd$protein_sig %in% TRUE)
  if (length(prot_idx) == 0L) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No significant proteins (p < 0.05)", size = 5) +
             theme_void())
  }
  gene_idx  <- unique(ed$from[ed$to %in% prot_idx & nd$type[ed$from] == "gene"])
  tf_up_idx <- unique(ed$from[ed$to %in% gene_idx  & nd$type[ed$from] == "tf"])
  
  keep_idx   <- sort(unique(c(prot_idx, gene_idx, tf_up_idx)))
  keep_names <- nd$name[keep_idx]
  
  g <- g %>% activate(nodes) %>% filter(name %in% keep_names) %>% activate(edges)
  
  # --- Plot --------------------------------------------------------------------
  n_nodes <- g %>% as_tibble(active = "nodes") %>% nrow()
  n_edges <- g %>% as_tibble(active = "edges") %>% nrow()
  
  p <- ggraph(g, layout = "stress") +
    # Edges colored by mean_log2fc (red=+1, blue=-1, NA=grey)
    geom_edge_link(
      aes(colour = mean_log2fc, linetype = direction),
      width = 1.5, alpha = 0.7,
      arrow = arrow(length = grid::unit(1, "mm")),
      start_cap = circle(3.5, "mm"), end_cap = circle(3.5, "mm")
    ) +
    
    scale_edge_linetype_manual(
      name   = "Direction",
      values = c("opposite" = "dashed",
                 "same" = "solid"),
      na.value = "solid"
    ) +
    
    guides(
      edge_linetype = guide_legend(
        keywidth = unit(0.8, "cm"),   # ✅ only way to widen the dashed line
        override.aes = list(
          linetype = c("dashed", "solid", "None"),
          color = "black",
          linewidth = 1.5
        )
      )
    ) +
    
    scale_edge_colour_gradient2(
      low = "blue", mid = "grey90", high = "red",
      midpoint = 0, limits = c(-0.7, 0.7), na.value = "black",
      name = "mean log2FC"
    ) +
    
    # TFs
    geom_node_point(
      data = ~ dplyr::filter(., type == "tf"),
      aes(size = node_size, shape = type), color = "black", stroke = 0.5, fill = NA
    ) +
    
    # Proteins (significant = green fill; NA grey)
    geom_node_point(
      data = ~ dplyr::filter(., type == "protein"),
      aes(size = node_size, shape = type),
      color = "black", stroke = 0.25, fill = "#4DAF4A"
    ) +
    
    ggnewscale::new_scale_fill() +
    
    # Genes (filled by gene log2FC)
    geom_node_point(
      data = ~ dplyr::filter(., type == "gene"),
      aes(size = node_size, fill = logfc, shape = type),color = "black", stroke = 0.25
    ) +
    
    # 4️⃣ Shape legend (manual)
    scale_shape_manual(
      name = "Node type",
      values = c(
        "gene" = 21,      # circle
        "tf" = 22,    # square
        "protein" = 24  # triangle
      )
    ) +
    
    scale_fill_gradient2(
      name = "Gene log2FC",
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1.5, 1.5), na.value = "#E0E0E0"
    ) +
    
    
    
    scale_size_identity(guide = "none") +
    ggrepel::geom_text_repel(
      data = ~ dplyr::filter(., !is.na(name)),
      aes(x = x, y = y, label = label),
      size = 4, max.overlaps = Inf, segment.size = 0.1,
      box.padding = 0.3, point.padding = 0.2, min.segment.length = 0, force = 2,
      bg.color = "white", bg.r = 0.15, fontface = "bold"
    ) +
    labs(title = paste0(cell, " Network")) +
    theme_void() +
    theme(
      text = element_text(family = "Helvetica", face = "bold", size = 6),
      plot.title = element_text(face = "bold", size = 6),
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 6),
      legend.key.size = unit(4, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.box.spacing = unit(1, "mm")
    )
  
  return(p)
}

Partial_graph_same_direction <- function(df, TF = "", cell = "", sig_protein_only = TRUE) {
  set.seed(42)
  
  # filter direction that is opposite
  df <- df %>% 
    filter(is.na(direction_match) | direction_match != "opposite")
  
  # 1) Build graph (must contain mean_log2fc on TF->gene edges, logfc on genes)
  g <- .build_master_graph(df)$g
  
  # 2) Node attributes (safe & consistent)
  g <- g %>%
    activate(nodes) %>%
    mutate(
      type      = tolower(node_type),
      degree    = centrality_degree(),
      node_size = case_when(
        type == "tf"      ~ 4,
        type == "gene"    ~ 4,
        type == "protein" ~ 4,
        TRUE              ~ 4
      )
    )
  
  # Derive/normalize sig_protein
  
  if (!"logfc" %in% names(as_tibble(g, active = "nodes"))) {
    g <- g %>% activate(nodes) %>% mutate(logfc = NA_real_)
  }
  
  # Keep only TF/Gene/Protein nodes
  g <- g %>% activate(nodes) %>% filter(type %in% c("tf","gene","protein"))
  
  # --- Multi-TF radius (≤2 hops) using tidygraph only --------------------------
  nd_tbl <- g %>% activate(nodes) %>% as_tibble()
  
  m <- .match_tf_nodes(g, TF)
  tf_idx <- m$idx
  if (length(tf_idx) == 0L) stop(sprintf("Requested TFs not found: %s", paste(TF, collapse = ", ")))
  
  # Compute min distance from any TF seed using tidygraph::node_distance_from
  dist_min <- rep(Inf, nrow(nd_tbl))
  for (i in tf_idx) {
    d <- g %>%
      activate(nodes) %>%
      mutate(d = node_distance_from(node = i, mode = "out")) %>%
      pull(d)
    # take elementwise min, ignoring NA
    dist_min <- pmin(dist_min, ifelse(is.na(d), Inf, d))
  }
  
  g <- g %>%
    activate(nodes) %>%
    mutate(dist_from_tf = dist_min) %>%
    filter(is.finite(dist_from_tf) & dist_from_tf <= 2) %>%
    activate(edges) %>%
    mutate(
      from_type = .N()$type[from],
      to_type   = .N()$type[to]
    ) %>%
    filter(from_type %in% c("tf","gene","protein") & to_type %in% c("tf","gene","protein"))
  
  # --- STRICT subgraph: only significant proteins & their upstreams ------------
  nd <- as_tibble(g, active = "nodes")
  ed <- as_tibble(g, active = "edges")
  
  prot_idx <- which(nd$type == "protein" & nd$protein_sig %in% TRUE)
  if (length(prot_idx) == 0L) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No significant proteins (p < 0.05)", size = 5) +
             theme_void())
  }
  gene_idx  <- unique(ed$from[ed$to %in% prot_idx & nd$type[ed$from] == "gene"])
  tf_up_idx <- unique(ed$from[ed$to %in% gene_idx  & nd$type[ed$from] == "tf"])
  
  keep_idx   <- sort(unique(c(prot_idx, gene_idx, tf_up_idx)))
  keep_names <- nd$name[keep_idx]
  
  g <- g %>% activate(nodes) %>% filter(name %in% keep_names) %>% activate(edges)
  
  # --- Plot --------------------------------------------------------------------
  n_nodes <- g %>% as_tibble(active = "nodes") %>% nrow()
  n_edges <- g %>% as_tibble(active = "edges") %>% nrow()
  
  p <- ggraph(g, layout = "stress") +
    # Edges colored by mean_log2fc (red=+1, blue=-1, NA=grey)
    geom_edge_link(
      aes(colour = mean_log2fc),
      width = 1.5, alpha = 0.7,
      arrow = arrow(length = grid::unit(1, "mm")),
      start_cap = circle(3.5, "mm"), end_cap = circle(3.5, "mm")
    ) +
    
    # scale_edge_linetype_manual(
    #   name   = "Direction",
    #   values = c("same" = "solid"),
    #   na.value = "solid"
    # ) +
    
    guides(
      edge_linetype = guide_legend(
        keywidth = unit(0.8, "cm"),   # ✅ only way to widen the dashed line
        override.aes = list(
          linetype = c("solid", "None"),
          color = "black",
          linewidth = 1.5
        )
      )
    ) +
    
    scale_edge_colour_gradient2(
      low = "blue", mid = "grey90", high = "red",
      midpoint = 0, limits = c(-0.7, 0.7), na.value = "black",
      name = "mean log2FC"
    ) +
    
    # TFs
    geom_node_point(
      data = ~ dplyr::filter(., type == "tf"),
      aes(size = node_size, shape = type), color = "black", stroke = 0.5, fill = NA
    ) +
    
    # Proteins (significant = green fill; NA grey)
    geom_node_point(
      data = ~ dplyr::filter(., type == "protein"),
      aes(size = node_size, shape = type),
      color = "black", stroke = 0.25, fill = "#4DAF4A"
    ) +
    
    ggnewscale::new_scale_fill() +
    
    # Genes (filled by gene log2FC)
    geom_node_point(
      data = ~ dplyr::filter(., type == "gene"),
      aes(size = node_size, fill = logfc, shape = type),color = "black", stroke = 0.25
    ) +
    
    # 4️⃣ Shape legend (manual)
    scale_shape_manual(
      name = "Node type",
      values = c(
        "gene" = 21,      # circle
        "tf" = 22,    # square
        "protein" = 24  # triangle
      )
    ) +
    
    scale_fill_gradient2(
      name = "Gene log2FC",
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1.5, 1.5), na.value = "#E0E0E0"
    ) +
    
    
    
    scale_size_identity(guide = "none") +
    ggrepel::geom_text_repel(
      data = ~ dplyr::filter(., !is.na(name)),
      aes(x = x, y = y, label = label),
      size = 4, max.overlaps = Inf, segment.size = 0.1,
      box.padding = 0.3, point.padding = 0.2, min.segment.length = 0, force = 2,
      bg.color = "white", bg.r = 0.15, fontface = "bold"
    ) +
    labs(title = paste0(cell, " Network")) +
    theme_void() +
    theme(
      text = element_text(family = "Helvetica", face = "bold", size = 6),
      plot.title = element_text(face = "bold", size = 6),
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 6),
      legend.key.size = unit(4, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.box.spacing = unit(1, "mm")
    )
  
  return(p)
}