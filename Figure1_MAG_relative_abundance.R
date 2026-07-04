# Figure1_MAG_relative_abundance.R
# Script to generate Figure1 of the manuscript.
# # =============================================================================
# Figure 1 — Relative abundance of hydrocarbon-degrading MAGs
# -----------------------------------------------------------------------------
# Stacked bar chart of MAG relative abundance (from RPKG/CPM), restricted to
# hydrocarbon-degrading MAGs in the top 15 taxonomic Classes, colored by Order,
# faceted by season. MAGs outside the top 15 Classes are pooled as "Others".
#
# Inputs:
#   data/MAG_abundance.csv   - abundance table (MAG column + one column per sample)
#   data/taxonomy.csv        - Genome, Domain..Species
#   data/metadata.csv        - Genome (= Sample ID), Bay, Season, Salinity_class
#   data/Hydrocarbon_MAGs.csv- Genome column listing HC-degrading MAG IDs
# Output:
#   Figure1_MAG_relative_abundance.pdf
# =============================================================================

library(tidyverse)
library(paletteer)
library(ggh4x)     # nested faceting (Season x Salinity)

# ---- input files ------------------------------------------------------------
abundance_file <- "data/MAG_abundance.csv"
taxonomy_file  <- "data/taxonomy.csv"
metadata_file  <- "data/metadata.csv"
hc_mags_file   <- "data/Hydrocarbon_MAGs.csv"

N_TOP_CLASSES  <- 15

# ---- 1. Load ----------------------------------------------------------------
abund <- read.csv(abundance_file, check.names = FALSE)
rownames(abund) <- abund$MAG
abund <- abund %>% select(-MAG)

taxonomy <- read.csv(taxonomy_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
hc_mags  <- read.csv(hc_mags_file,  check.names = FALSE)$Genome

# ---- 2. Relative abundance (% per sample) -----------------------------------
rel <- abund %>% mutate(across(everything(), ~ . / sum(.) * 100))

# ---- 3. Keep only hydrocarbon-degrading MAGs --------------------------------
rel_hc <- rel[rownames(rel) %in% hc_mags, ]

# ---- 4. Top 15 Classes; pool everything else as "Others" --------------------
top_classes <- rel_hc %>%
  rownames_to_column("Genome") %>%
  left_join(taxonomy %>% select(Genome, Class), by = "Genome") %>%
  pivot_longer(-c(Genome, Class), names_to = "Sample", values_to = "Abundance") %>%
  group_by(Class, Sample) %>%
  summarise(tot = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Class) %>%
  summarise(mean_ab = mean(tot), .groups = "drop") %>%
  slice_max(mean_ab, n = N_TOP_CLASSES) %>%
  pull(Class)

in_top  <- rel_hc[rownames(rel_hc) %in%
                    (taxonomy %>% filter(Class %in% top_classes) %>% pull(Genome)), ]
out_top <- rel_hc[!rownames(rel_hc) %in% rownames(in_top), ]

plot_mat <- in_top
if (nrow(out_top) > 0) plot_mat <- rbind(plot_mat, Others = colSums(out_top))

# ---- 5. Long format + taxonomy + metadata -----------------------------------
tax <- taxonomy %>% select(Genome, Class, Order) %>%
  bind_rows(data.frame(Genome = "Others", Class = "Others", Order = "Others"))

plot_data <- plot_mat %>%
  rownames_to_column("Genome") %>%
  pivot_longer(-Genome, names_to = "Sample", values_to = "Abundance") %>%
  left_join(tax, by = "Genome") %>%
  left_join(metadata %>% select(Genome, Bay, Season, Salinity_class),
            by = c("Sample" = "Genome")) %>%
  mutate(
    Bay            = factor(Bay, levels = c("CP", "DE")),
    Season         = factor(Season, levels = c("Spring", "Summer", "Fall")),
    Salinity_class = factor(Salinity_class, levels = c("Low", "Medium", "High")),
    Class          = replace_na(Class, "Others"),
    Order          = replace_na(Order, "Others")
  )

# order samples by season then bay; order MAGs by class/order/abundance
sample_order <- plot_data %>% distinct(Sample, Season, Bay) %>%
  arrange(Season, Bay) %>% pull(Sample)
plot_data$Sample <- factor(plot_data$Sample, levels = sample_order)

mag_order <- plot_data %>%
  group_by(Genome, Class, Order) %>%
  summarise(m = mean(Abundance), .groups = "drop") %>%
  arrange(Class, Order, desc(m)) %>% pull(Genome)
plot_data$Genome <- factor(plot_data$Genome, levels = rev(mag_order))

# ---- 6. Colors: one per Order (nested in Class), Others = grey ---------------
orders <- plot_data %>% filter(Order != "Others") %>%
  distinct(Class, Order) %>% arrange(Class, Order) %>% pull(Order)
pal <- colorRampPalette(as.character(paletteer_d("ggsci::nrc_npg")))(length(orders))
order_colors <- c(setNames(pal, orders), Others = "grey80")
plot_data$Order <- factor(plot_data$Order, levels = c(orders, "Others"))

# ---- 7. Plot ----------------------------------------------------------------
p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_grid(. ~ Season, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = order_colors, name = "Order") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40),
                     breaks = seq(0, 40, 10)) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.y       = element_text(face = "bold"),
    strip.text         = element_text(face = "bold"),
    strip.background    = element_rect(fill = "grey90", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA),
    legend.key.size    = unit(0.4, "cm"),
    legend.text        = element_text(size = 8)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave("Figure1_MAG_relative_abundance.pdf", p, width = 16, height = 8)
cat("Figure1_MAG_relative_abundance.pdf saved.\n")
