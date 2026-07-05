# Figure5_gene_environment_correlation.R
# Script to generate Figure5 of the manuscript.
# # =============================================================================
# Figure 5 — Correlation of HC-degradation gene abundance with environmental data
# -----------------------------------------------------------------------------
# Spearman correlations between per-gene abundance (RPKG) and environmental
# variables across samples. FDR-adjusted significance is shown as asterisks;
# genes are grouped by degradation pathway (row annotation).
#
# Inputs:
#   data/gene_RPKG.csv        - gene abundance (Gene column + one col per sample)
#   data/metadata.csv         - sample metadata (samples as rows; numeric env vars)
#   data/HC_genes_pathway_map.csv - columns: gene, pathway
# Output:
#   Figure5_gene_environment_correlation.pdf
#   gene_metadata_correlation_spearman.csv  (correlation/FDR table)
# =============================================================================

library(tidyverse)
library(pheatmap)
library(paletteer)

# ---- input files ------------------------------------------------------------
gene_file     <- "data/gene_RPKG.csv"
metadata_file <- "data/metadata.csv"
pathway_file  <- "data/HC_genes_pathway_map.csv"

# ---- 1. Load & align --------------------------------------------------------
feature_df <- read.csv(gene_file, header = TRUE, check.names = FALSE)
rownames(feature_df) <- feature_df$Gene
feature_df <- feature_df[, -1]

metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
common   <- intersect(colnames(feature_df), rownames(metadata))
feature_df <- feature_df[, common]
metadata   <- metadata[common, ]
metadata_num <- metadata %>% select(where(is.numeric))

# ---- 2. Spearman correlation + FDR ------------------------------------------
cor_mat <- matrix(NA, nrow(feature_df), ncol(metadata_num),
                  dimnames = list(rownames(feature_df), colnames(metadata_num)))
p_mat <- cor_mat
for (i in seq_len(nrow(feature_df))) {
  for (j in seq_len(ncol(metadata_num))) {
    test <- suppressWarnings(
      cor.test(as.numeric(feature_df[i, ]), metadata_num[[j]], method = "spearman"))
    cor_mat[i, j] <- test$estimate
    p_mat[i, j]   <- test$p.value
  }
}
p_adj <- matrix(p.adjust(p_mat, method = "fdr"), nrow(p_mat), ncol(p_mat),
                dimnames = dimnames(p_mat))

# ---- 3. Pathway row annotation ----------------------------------------------
map <- read.csv(pathway_file, check.names = FALSE)
annotation_row <- data.frame(
  Pathway = map$pathway[match(rownames(cor_mat), map$gene)],
  row.names = rownames(cor_mat))

pathways <- unique(annotation_row$Pathway)
pal <- if (length(pathways) <= 12)
  as.character(paletteer_d("tidyquant::tq_light", n = length(pathways))) else
  colorRampPalette(as.character(paletteer_d("RColorBrewer::Set3")))(length(pathways))
pathway_colors <- setNames(pal, pathways)

# ---- 4. Significance labels (FDR) -------------------------------------------
sig <- ifelse(p_adj < 0.001, "***",
       ifelse(p_adj < 0.01,  "**",
       ifelse(p_adj < 0.05,  "*", "")))

# ---- 5. Order genes by pathway ----------------------------------------------
ord <- order(annotation_row$Pathway)
cor_ord <- cor_mat[ord, ]
sig_ord <- sig[ord, ]
ann_ord <- annotation_row[ord, , drop = FALSE]

# ---- 6. Heatmap -------------------------------------------------------------
pdf("Figure5_gene_environment_correlation.pdf", width = 10, height = 12)
pheatmap(
  cor_ord,
  cluster_rows    = FALSE,
  cluster_cols    = TRUE,
  main            = "Gene vs environmental metadata correlation (Spearman)",
  fontsize_row    = 10,
  fontsize_col    = 10,
  display_numbers = sig_ord,
  number_color    = "black",
  color           = rev(as.character(paletteer_c("ggthemes::Classic Red-Blue", 10))),  # reversed: blue = positive
  breaks          = seq(-1, 1, length.out = 11),
  na_col          = "white",
  border_color    = "white",
  annotation_row    = ann_ord,
  annotation_colors = list(Pathway = pathway_colors)
)
dev.off()

# ---- 7. Save correlation table ----------------------------------------------
tidy <- expand.grid(Gene = rownames(cor_mat), Metadata = colnames(cor_mat)) %>%
  mutate(Correlation = as.vector(cor_mat),
         P_value     = as.vector(p_mat),
         FDR         = as.vector(p_adj)) %>%
  arrange(FDR)
write.csv(tidy, "gene_metadata_correlation_spearman.csv", row.names = FALSE)

cat("Figure5_gene_environment_correlation.pdf saved.\n")
