# Figure2_pathway_abundance_order.R
# Script to generate Figure2 of the manuscript.
# # =============================================================================
# Figure 2 — Hydrocarbon-degradation pathway abundance at the Order level
# -----------------------------------------------------------------------------
# Heatmap of mean pathway abundance (log2 RPKG) summarized by taxonomic Order.
# Abundance is summed across genes within each pathway, averaged across samples,
# and aggregated across MAGs within each Order. Rows = pathways, columns = Orders
# (both clustered, Ward.D2 on Euclidean distance).
#
# Input:
#   data/pathway_RPKG.csv - columns: MAG, Order, Pathway, and one column per
#                           sample (sample columns start with "DE" or "CP")
# Output:
#   Figure2_pathway_abundance.pdf
# =============================================================================

library(tidyverse)
library(pheatmap)

# ---- input file -------------------------------------------------------------
input_file <- "data/pathway_RPKG.csv"

# ---- 1. Load & reshape to long ----------------------------------------------
df <- read.csv(input_file, check.names = FALSE)

df_long <- df %>%
  pivot_longer(cols = matches("^(DE|CP)"),
               names_to = "Sample", values_to = "Abundance")

# ---- 2. Mean abundance per Order x Pathway ----------------------------------
# average across samples, then aggregate across MAGs within each Order
df_summary <- df_long %>%
  group_by(Order, Pathway, Sample) %>%
  summarise(sample_ab = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Order, Pathway) %>%
  summarise(mean_abundance = mean(sample_ab, na.rm = TRUE) * 100, .groups = "drop")

# ---- 3. Matrix (pathways x Orders), log2-transformed ------------------------
heatmap_mat <- df_summary %>%
  pivot_wider(names_from = Order, values_from = mean_abundance, values_fill = 0) %>%
  column_to_rownames("Pathway")

heatmap_mat_log <- log2(as.matrix(heatmap_mat) + 1)

# ---- 4. Colors --------------------------------------------------------------
my_colors <- colorRampPalette(
  c("white", "#abd9e9", "#74add1", "#2171b5", "#08306b",
    "#fdae61", "#f16913", "#cb181d"))(256)

# ---- 5. Heatmap -------------------------------------------------------------
pdf("Figure2_pathway_abundance.pdf", width = 12, height = 9)
pheatmap(
  heatmap_mat_log,
  scale                    = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  cluster_rows             = TRUE,
  cluster_cols             = TRUE,
  color                    = my_colors,
  border_color             = "grey97",
  fontsize                 = 10,
  treeheight_row           = 50,
  treeheight_col           = 50,
  cellwidth                = 22,
  cellheight               = 26,
  main = "Pathway abundance by Order (log2 RPKG)"
)
dev.off()

cat("Figure2_pathway_abundance.pdf saved.\n")
