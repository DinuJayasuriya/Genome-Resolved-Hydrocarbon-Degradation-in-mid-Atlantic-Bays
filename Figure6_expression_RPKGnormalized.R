# ======================================================
# VERSION A - FINAL (with SEPARATE pathway file, as in original)
# RPKG-normalized HC gene expression heatmap
#   - zeros kept (offset 0.01) -> relative expression ACROSS ALL conditions
#   - row Z-score
#   - palette recentered: white = 0 (gene average); blue below, warm above
#   - colors mapped to real Z range via breaks; genuine NA shown grey
#   - PATHWAYS loaded from a SEPARATE pathway file (gene + pathway columns)
#   - your colors, grid, pathway grouping, no gaps, no clustering
# ======================================================

library(pheatmap)
library(RColorBrewer)
library(paletteer)

# ---- knobs ----
OFFSET         <- 0.01
CONTRAST_LIMIT <- 2

# ======================================================
# 1. Load expression data (RPKG-normalized)
# ======================================================
cat("Select RPKG-normalized file (HC_RPKGnormalized_GENELEVEL.csv):\n")
expr_in <- read.csv(file.choose(), row.names = 1, check.names = FALSE)
# NOTE: row.names = 1 -> Gene column becomes rownames (matches your original style).
# If the file's first column is 'Gene', this makes genes the rownames.
# Drop the 'pathway' column if present (we'll get pathways from the separate file).
if ("pathway" %in% colnames(expr_in)) expr_in$pathway <- NULL
expr_data <- as.matrix(expr_in)
cat("Genes:", nrow(expr_data), "| samples:", ncol(expr_data), "\n")

# ======================================================
# 2. Metadata
# ======================================================
cat("\nSelect metadata file:\n")
metadata <- read.csv(file.choose(), row.names = 1, check.names = FALSE)
metadata$Bay      <- factor(metadata$Bay,      levels = c("CP","DE"))
metadata$Season   <- factor(metadata$Season,   levels = c("Spr","Sum","Fall"))
metadata$Salinity <- factor(metadata$Salinity, levels = c("Low","Med","High"))
metadata <- metadata[!is.na(metadata$Salinity), ]
common <- intersect(colnames(expr_data), rownames(metadata))
expr_data <- expr_data[, common, drop = FALSE]; metadata <- metadata[common, ]
cat("Samples after alignment:", ncol(expr_data), "\n")

# ======================================================
# 3. Pathway annotations  (SEPARATE FILE - as in your original)
#    File must have a 'gene' column and a 'pathway' column.
# ======================================================
cat("\nSelect pathway file (gene + pathway columns):\n")
pathway_lookup <- NULL
tryCatch({
  pathway_info <- read.csv(file.choose(), check.names = FALSE)
  if ("gene" %in% colnames(pathway_info)) {
    pathway_info <- pathway_info[!duplicated(pathway_info$gene), ]
    pathway_lookup <- setNames(pathway_info$pathway, pathway_info$gene)
    cat("Genes with pathway annotations:", length(pathway_lookup), "\n")
  } else {
    cat("ERROR: no 'gene' column in pathway file. Columns:",
        paste(colnames(pathway_info), collapse=", "), "\n")
  }
}, error = function(e) cat("Skipping pathway annotations\n"))

# ======================================================
# 4. Transform (zeros KEPT via offset; NAs stay NA)
# ======================================================
expr_log <- log10(expr_data + OFFSET)
all_na <- apply(expr_log, 1, function(x) all(is.na(x)))
if (any(all_na)) { cat("Dropping all-NA genes:", sum(all_na), "\n"); expr_log <- expr_log[!all_na,,drop=FALSE] }

# ======================================================
# 5. Column order (Bay > Season > Salinity), col annotation
# ======================================================
annotation_col <- metadata[, c("Bay","Season","Salinity")]
col_order <- rownames(annotation_col[order(annotation_col$Bay, annotation_col$Season, annotation_col$Salinity), ])
expr_log <- expr_log[, col_order, drop = FALSE]; annotation_col <- annotation_col[col_order, ]

# ======================================================
# 6. Row (pathway) annotation - matched from separate file
# ======================================================
annotation_row <- NULL
if (!is.null(pathway_lookup)) {
  annotation_row <- data.frame(
    Pathway = pathway_lookup[rownames(expr_log)],
    row.names = rownames(expr_log))
  miss <- is.na(annotation_row$Pathway)
  if (any(miss)) {
    cat("Genes with NO pathway match:", sum(miss), "\n")
    cat("  -> ", paste(rownames(annotation_row)[miss], collapse=", "), "\n")
    annotation_row$Pathway[miss] <- "Unknown"
  }
  cat("\nPathway distribution:\n"); print(sort(table(annotation_row$Pathway), decreasing=TRUE))
}

# ======================================================
# 7. Colors
# ======================================================
ann_colors <- list(
  Bay      = c(CP="midnightblue", DE="#F08080"),
  Season   = c(Spr="#98D8C8", Sum="#F7DC6F", Fall="#E67E22"),
  Salinity = c(Low="#D6EAF8", Med="#5DADE2", High="#1B4F72"))
if (!is.null(annotation_row))
  ann_colors$Pathway <- setNames(
    paletteer::paletteer_d("tidyquant::tq_light", length(unique(annotation_row$Pathway))),
    unique(annotation_row$Pathway))

# recentered palette: white = 0 (gene average)
cols <- colorRampPalette(c("#08306b","#2171b5","#74add1","#abd9e9",
                           "white",
                           "#fdae61","#f16913","#cb181d"))(256)
breaks <- seq(-CONTRAST_LIMIT, CONTRAST_LIMIT, length.out = length(cols) + 1)

# ======================================================
# 8. Order genes by pathway (grouped, NO gap)
# ======================================================
if (!is.null(annotation_row)) {
  gd <- data.frame(gene = rownames(expr_log),
                   pathway = annotation_row$Pathway,
                   me = rowMeans(expr_log, na.rm = TRUE))
  gd <- gd[order(gd$pathway, -gd$me), ]
  expr_log       <- expr_log[gd$gene, , drop = FALSE]
  annotation_row <- annotation_row[gd$gene, , drop = FALSE]
}

# ======================================================
# 9. Heatmap
# ======================================================
hm <- function() pheatmap(
  expr_log,
  annotation_col = annotation_col, annotation_row = annotation_row,
  annotation_colors = ann_colors,
  scale = "row", breaks = breaks, color = cols, na_col = "grey50",
  show_rownames = TRUE, show_colnames = FALSE,
  main = "HC gene expression (TPM / carrier RPKG) - row Z-score",
  cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = FALSE,
  fontsize_row = 7, fontsize_col = 8, fontsize = 11,
  border_color = "grey20", cellwidth = 12, cellheight = 10, legend = TRUE)
hm()

# ======================================================
# 10. Save (optional)
# ======================================================
cat("\nSave figure? (y/n): ")
if (tolower(readline()) %in% c("y","yes")) {
  png("Heatmap_A_RPKGnormalized.png", width = 12, height = 11, units = "in", res = 300); hm(); dev.off()
  pdf("Heatmap_A_RPKGnormalized.pdf", width = 12, height = 11); hm(); dev.off()
  cat("Saved PNG + PDF\n")
}
cat("\nDone.\n")
