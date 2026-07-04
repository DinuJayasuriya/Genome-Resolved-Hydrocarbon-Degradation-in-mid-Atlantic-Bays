# =============================================================================
# DNA-abundance normalization of hydrocarbon-degradation gene expression
# GENE-LEVEL output, RPKG denominator
# -----------------------------------------------------------------------------
# For each GENE (collapsed across all its MAG copies) in each SAMPLE:
#
#   Normalized(gene, sample) =        sum of TPM of that gene's transcripts
#                              ---------------------------------------------------
#                              sum of RPKG of all MAGs carrying that gene, in sample
#
# Why RPKG (not CPM): RPKG is a per-genome-equivalent abundance measure suited
# to normalizing expression to population size, and it resolves low-abundance
# carrier MAGs that CPM rounds to zero (reduces undefined/NA cells from ~41% to <1%).
#
# Inputs:
#   1. average_TPM_geneid_MAG.csv  (HC genes; per-sample TPM; cols Geneid,kegg_hit,pathway,Gene,...)
#   2. RPKG_MAG.csv                (MAG x sample RPKG; col 'MAG' + sample cols)
# Output:
#   HC_RPKGnormalized_GENELEVEL.csv
#     - one row per Gene
#     - columns: Gene, pathway, + 32 shared sample columns
#     - 0  = gene silent (carriers present, no transcripts)
#     - NA = carriers absent in that sample (undefined)
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)

# ---- 1. Load -----------------------------------------------------------------
tpm  <- read.csv("average_TPM_geneid_MAG.csv", check.names = FALSE)
rpkg <- read.csv("RPKG_MAG.csv",              check.names = FALSE)

meta_cols       <- c("Geneid", "kegg_hit", "pathway", "Gene")
tpm_sample_cols <- setdiff(colnames(tpm), meta_cols)

# ---- 2. Match sample names ---------------------------------------------------
# RPKG uses DE_/CP_ (same as TPM). TPM has two columns ending in _RNA -> strip it.
strip_rna <- function(s) str_replace(s, "_RNA$", "")
new_names <- sapply(tpm_sample_cols, strip_rna)
colnames(tpm)[match(tpm_sample_cols, colnames(tpm))] <- new_names
tpm_sample_cols <- unname(new_names)

rpkg_sample_cols <- setdiff(colnames(rpkg), "MAG")
shared_samples   <- intersect(tpm_sample_cols, rpkg_sample_cols)
cat("Shared samples:", length(shared_samples), "\n")
if (length(shared_samples) == 0) stop("No shared samples.")

# ---- 3. Extract MAG from Geneid ---------------------------------------------
tpm <- tpm %>% mutate(MAG = str_extract(Geneid, "^[^_]+_MAG_\\d+"))

missing_mags <- setdiff(unique(tpm$MAG), rpkg$MAG)
if (length(missing_mags) > 0)
  warning("HC carrier MAGs missing from RPKG (treated as 0): ",
          paste(missing_mags, collapse = ", "))

# ---- 4. Gene -> pathway (1:1) ------------------------------------------------
gene_path <- tpm %>% select(Gene, pathway) %>% distinct()

# ---- 5. Denominator: summed RPKG of carrier MAGs, per Gene x Sample ----------
gene_to_mags <- tpm %>% select(Gene, MAG) %>% distinct()

rpkg_long <- rpkg %>%
  pivot_longer(all_of(shared_samples), names_to = "Sample", values_to = "RPKG")

denom <- gene_to_mags %>%
  left_join(rpkg_long, by = "MAG", relationship = "many-to-many") %>%
  group_by(Gene, Sample) %>%
  summarise(denom = sum(RPKG, na.rm = TRUE), .groups = "drop")

# ---- 6. Numerator: summed TPM per Gene x Sample -----------------------------
num <- tpm %>%
  select(Gene, all_of(shared_samples)) %>%
  pivot_longer(all_of(shared_samples), names_to = "Sample", values_to = "TPM") %>%
  group_by(Gene, Sample) %>%
  summarise(sumTPM = sum(TPM, na.rm = TRUE), .groups = "drop")

# ---- 7. Divide once ----------------------------------------------------------
norm <- num %>%
  left_join(denom, by = c("Gene", "Sample")) %>%
  mutate(norm = if_else(denom > 0, sumTPM / denom, NA_real_))

# ---- 8. Wide: one row per Gene, pathway attached ----------------------------
out <- norm %>%
  select(Gene, Sample, norm) %>%
  pivot_wider(names_from = Sample, values_from = norm) %>%
  left_join(gene_path, by = "Gene") %>%
  relocate(pathway, .after = Gene)

# ---- 9. Save -----------------------------------------------------------------
write.csv(out, "HC_RPKGnormalized_GENELEVEL.csv", row.names = FALSE)
cat("Saved: HC_RPKGnormalized_GENELEVEL.csv\n")
cat("Genes:", nrow(out), "| Sample cols:", length(shared_samples), "\n")
na_frac <- mean(is.na(as.matrix(out[, shared_samples])))
cat(sprintf("NA cells (carriers absent): %.2f%%\n", 100 * na_frac))

# =============================================================================
# HEATMAP (optional) — log + per-gene Z-score, NA shown as grey, rows by pathway
# =============================================================================
# library(pheatmap)
# mat <- as.matrix(out[, shared_samples]); rownames(mat) <- out$Gene
# mat_log <- log10(mat + 1e-6)                 # small offset; values are small
# mat_z   <- t(scale(t(mat_log)))              # per-gene Z-score across samples
# row_ann <- data.frame(Pathway = out$pathway); rownames(row_ann) <- out$Gene
# pheatmap(mat_z,
#          cluster_rows = TRUE, cluster_cols = TRUE,
#          annotation_row = row_ann,
#          na_col = "grey90",                  # NA = carriers absent
#          main = "HC gene expression (TPM / carrier RPKG), per-gene Z-score")
