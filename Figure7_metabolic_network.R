# Figure7_metabolic_network.R
# Script to generate Figure7 of the manuscript.
# # =============================================================================
# Figure 7 — Co-expression networks of hydrocarbon-degrading MAGs
# -----------------------------------------------------------------------------
# For each focal taxon (Burkholderiales, Pseudomonadales), builds a bipartite
# MAG x metabolic-category network from summed, log2-transformed expression.
# Nodes are MAGs (green) or metabolic categories (gold); an edge connects a MAG
# to a category it expresses, edge width = normalized mean expression of that
# category in the MAG, and category node size = number of connected MAGs (degree).
# (Panel A of Figure 7, the marker-gene presence/absence heatmap, is generated
#  separately.)
#
# Input (per taxon): an Excel file with columns
#   mags, order, gene.name, category, pathway, + one column per sample.
# Output: a co-expression network figure per taxon.
# =============================================================================

library(readxl)
library(dplyr)
library(stringr)
library(igraph)
library(ggraph)
library(ggplot2)

# ---- input files (edit paths) -----------------------------------------------
taxon_files <- c(
  Burkholderiales = "data/Burkholderiales_gene_summed.xlsx",
  Pseudomonadales = "data/Pseudomonadales_gene_summed.xlsx"
)

# five metabolic categories, fixed order and colors
category_levels <- c("Hydrocarbon degradation", "C1 metabolism",
                     "Complex carbon degradation", "Fermentation",
                     "Central carbon metabolism")

# =============================================================================
# Build a co-expression network for one taxon
# =============================================================================
build_network <- function(path, taxon_name) {

  df <- read_excel(path)
  names(df) <- tolower(gsub(" ", "_", names(df)))

  meta_cols  <- c("mags", "order", "gene.name", "category", "pathway")
  count_cols <- setdiff(names(df), meta_cols)

  # numeric counts; stop if any non-numeric
  df[count_cols] <- lapply(df[count_cols], as.numeric)
  if (any(is.na(df[count_cols]))) stop("NA values in count columns of ", path)

  # recode category into the five standard groups
  cat_clean <- df$category %>% str_trim() %>% str_to_lower() %>% str_replace_all("\\s+", " ")
  df$category <- factor(case_when(
    str_detect(cat_clean, "hydrocarbon") ~ "Hydrocarbon degradation",
    str_detect(cat_clean, "^c1")         ~ "C1 metabolism",
    str_detect(cat_clean, "complex")     ~ "Complex carbon degradation",
    str_detect(cat_clean, "ferment")     ~ "Fermentation",
    str_detect(cat_clean, "central")     ~ "Central carbon metabolism",
    TRUE                                 ~ NA_character_
  ), levels = category_levels)

  # log2 expression, MAG number
  mat_log <- log2(as.matrix(df[, count_cols]) + 1)
  df$mag_num <- as.numeric(str_extract(tolower(df$mags), "(?<=_mag_)\\d+"))

  # mean expression per MAG x category  ->  bipartite edges
  df$row_mean <- rowMeans(mat_log)
  agg <- df %>%
    group_by(mag_num, category) %>%
    summarise(mean_expr = mean(row_mean), .groups = "drop")

  edges <- data.frame(
    from   = as.character(agg$mag_num),
    to     = as.character(agg$category),
    weight = agg$mean_expr / max(agg$mean_expr)
  )
  nodes <- data.frame(name = unique(c(edges$from, edges$to)))
  nodes$type <- ifelse(nodes$name %in% edges$from, "MAG", "Category")

  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  V(g)$degree <- igraph::degree(g)

  # ---- plot: category labels bold, MAG labels small; node size = degree ----
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = weight, alpha = weight), color = "grey40") +
    geom_node_point(aes(color = type, size = degree)) +
    geom_node_text(aes(label = ifelse(type == "MAG", name, NA)),
                   repel = TRUE, size = 3) +
    geom_node_text(aes(label = ifelse(type == "Category", name, NA)),
                   repel = FALSE, fontface = "bold", size = 4.5, color = "black") +
    scale_size(range = c(4, 16)) +
    scale_edge_width(range = c(0.2, 2)) +
    scale_color_manual(values = c(MAG = "darkgreen", Category = "goldenrod")) +
    theme_void() +
    ggtitle(paste0(taxon_name, " co-expression network"))

  ggsave(paste0("Figure7_network_", taxon_name, ".pdf"), p, width = 9, height = 8)
  cat("Saved Figure7_network_", taxon_name, ".pdf\n", sep = "")
  invisible(g)
}

# =============================================================================
# Run for each taxon
# =============================================================================
for (taxon in names(taxon_files)) {
  build_network(taxon_files[[taxon]], taxon)
}
