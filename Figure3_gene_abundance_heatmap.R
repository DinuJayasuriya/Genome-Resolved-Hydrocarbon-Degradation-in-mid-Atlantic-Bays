# Figure3_gene_abundance_heatmap.R
# Script to generate Figure3 of the manuscript.
# # =============================================================================
# Figure 3 — Abundance of hydrocarbon-degradation genes within MAGs
# -----------------------------------------------------------------------------
# Heatmap of per-gene abundance (log10 RPKG) across samples. Genes are grouped
# by degradation pathway (row annotation); samples are ordered Chesapeake (CP)
# then Delaware (DE).
#
# Input:
#   data/gene_pathway_RPKG.csv - columns: Genes, Pathway, and one column per
#                                sample (sample columns start with "CP" or "DE")
# Output:
#   Figure3_gene_abundance_heatmap.pdf
# =============================================================================

library(tidyverse)
library(pheatmap)
library(paletteer)

# ---- input file -------------------------------------------------------------
input_file <- "data/gene_pathway_RPKG.csv"

# ---- 1. Load ----------------------------------------------------------------
rpkg <- read.csv(input_file, check.names = FALSE)
rownames(rpkg) <- rpkg$Genes

pathway_annot <- data.frame(Pathway = rpkg$Pathway, row.names = rownames(rpkg))
rpkg_data <- rpkg %>% select(-Genes, -Pathway)

# ---- 2. Transform: RPKG x100, then log10 ------------------------------------
rpkg_log10 <- log10(rpkg_data * 100 + 1)

# ---- 3. Order rows by pathway, columns CP then DE ---------------------------
ord <- order(pathway_annot$Pathway)
rpkg_log10    <- rpkg_log10[ord, ]
pathway_annot <- pathway_annot[ord, , drop = FALSE]

cp_cols <- grep("^CP", colnames(rpkg_log10), value = TRUE)
de_cols <- grep("^DE", colnames(rpkg_log10), value = TRUE)
rpkg_log10 <- rpkg_log10[, c(cp_cols, de_cols), drop = FALSE]

# ---- 4. Colors --------------------------------------------------------------
cols <- colorRampPalette(
  c("white", "#abd9e9", "#74add1", "#2171b5", "#08306b",
    "#fdae61", "#f16913", "#cb181d"))(256)

pathway_colors <- setNames(
  as.character(paletteer_d("tidyquant::tq_light", length(unique(pathway_annot$Pathway)))),
  unique(pathway_annot$Pathway))

# ---- 5. Heatmap -------------------------------------------------------------
pdf("Figure3_gene_abundance_heatmap.pdf", width = 10, height = 12)
pheatmap(
  as.matrix(rpkg_log10),
  color             = cols,
  cluster_rows      = FALSE,      # keep pathway order
  cluster_cols      = FALSE,      # keep CP-then-DE column order
  annotation_row    = pathway_annot,
  annotation_colors = list(Pathway = pathway_colors),
  fontsize_row      = 8,
  fontsize_col      = 10,
  border_color      = "white",
  cellwidth         = 10,
  cellheight        = 10,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  main              = ""
)
dev.off()

cat("Figure3_gene_abundance_heatmap.pdf saved.\n")
