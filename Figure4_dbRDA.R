# Figure4_dbRDA.R
# Script to generate Figure4 of the manuscript.
# # =============================================================================
# Figure 4 — Distance-based redundancy analysis (db-RDA) of HC-degrading MAGs
# -----------------------------------------------------------------------------
# Constrained ordination of MAG community structure against environmental
# variables (robust Aitchison / rclr distance). Points are colored by season
# and shaped by salinity class; significant environmental predictors are drawn
# as biplot vectors.
#
# Inputs:
#   data/metadata.csv   - sample metadata (Season, Salinity_class, Bay + env vars)
#   data/MAG_CPM.csv    - MAG CPM table (MAGs as rows, samples as columns)
# Output:
#   Figure4_dbRDA.pdf
#
# Requires vegan >= 2.6 (robust Aitchison distance).
# =============================================================================

library(vegan)

# ---- input files ------------------------------------------------------------
metadata_file <- "data/metadata.csv"
cpm_file      <- "data/MAG_CPM.csv"

# ---- 1. Load & prepare ------------------------------------------------------
meta   <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
sp_raw <- read.csv(cpm_file,      row.names = 1, check.names = FALSE)

meta$Season         <- factor(meta$Season,         levels = c("Spring","Summer","Fall"))
meta$Salinity_class <- factor(meta$Salinity_class, levels = c("Low","Medium","High"))
meta$Bay            <- factor(meta$Bay,            levels = c("Chesapeake","Delaware"))

colnames(sp_raw) <- trimws(gsub("Bay_", "_", colnames(sp_raw)))
sp <- t(as.matrix(sp_raw))
common <- intersect(rownames(meta), rownames(sp))
meta <- meta[common, , drop = FALSE]
sp   <- sp[common, , drop = FALSE]

# standardized environmental predictors
env_vars <- c("Temp","Nitrate","Phosphate","Silicate",
              "Bacterial_production","Cells/mL","Salinity")
chemistry <- meta[, env_vars, drop = FALSE]
colnames(chemistry)[colnames(chemistry) == "Cells/mL"] <- "Cells_mL"
chem <- decostand(chemistry, "stand")
chem$Bay <- meta$Bay

# ---- 2. db-RDA model --------------------------------------------------------
dbRDA <- capscale(sp ~ Temp + Nitrate + Phosphate + Silicate +
                    Bacterial_production + Cells_mL + Salinity + Bay,
                  data = chem, distance = "robust.aitchison")

# ---- 3. Extract scores ------------------------------------------------------
perc  <- round(100 * summary(dbRDA)$cont$importance[2, 1:2], 2)
sc_si <- scores(dbRDA, display = "sites",   scaling = 1)
sc_sp <- scores(dbRDA, display = "species", scaling = 1)
sc_bp <- scores(dbRDA, display = "bp",      scaling = 1)

# top 11 high-loading MAGs to label
strength <- apply(abs(sc_sp[,1:2]), 1, max)
top_mags <- names(sort(strength, decreasing = TRUE))[1:11]

# ---- 4. Plot settings -------------------------------------------------------
season_cols  <- c(Spring = "#FFFACD", Summer = "#87CEEB", Fall = "#FFB6C1")
salinity_pch <- c(Low = 21, Medium = 22, High = 24)

all_x <- c(sc_sp[,1], sc_si[,1]); all_y <- c(sc_sp[,2], sc_si[,2])
xpad <- diff(range(all_x))*0.18;  ypad <- diff(range(all_y))*0.18
xlim <- c(min(all_x)-xpad, max(all_x)+xpad)
ylim <- c(min(all_y)-ypad, max(all_y)+ypad)

# Bay is categorical -> in the model, but not drawn as an arrow
sc_bp_cont <- sc_bp[!grepl("^Bay", rownames(sc_bp)), , drop = FALSE]
plot_r  <- min(diff(xlim), diff(ylim))/2
ascale  <- (plot_r*0.85) / max(sqrt(sc_bp_cont[,1]^2 + sc_bp_cont[,2]^2))
ax <- sc_bp_cont[,1]*ascale; ay <- sc_bp_cont[,2]*ascale

# ---- 5. Draw ----------------------------------------------------------------
pdf("Figure4_dbRDA.pdf", width = 10, height = 9)
par(mar = c(5,5,4,11))
plot(NA, xlim = xlim, ylim = ylim,
     xlab = paste0("CAP1 (", perc[1], "%)"),
     ylab = paste0("CAP2 (", perc[2], "%)"),
     main = "db-RDA: MAG community structure")

ordihull(sc_si, groups = meta$Season, draw = "polygon",
         col = season_cols[levels(meta$Season)],
         border = season_cols[levels(meta$Season)],
         alpha = 90, lty = "dotted", lwd = 2, label = FALSE)

points(sc_sp, pch = 24, bg = "firebrick", col = "white", cex = 1.1)

ang <- atan2(sc_sp[top_mags,2], sc_sp[top_mags,1]); nudge <- diff(xlim)*0.05
text(sc_sp[top_mags,1] + cos(ang)*nudge, sc_sp[top_mags,2] + sin(ang)*nudge,
     labels = top_mags, col = "grey15", cex = 0.62, font = 2)

for (s in levels(meta$Season)) for (sal in levels(meta$Salinity_class)) {
  idx <- which(meta$Season == s & meta$Salinity_class == sal)
  if (length(idx) > 0)
    points(sc_si[idx,,drop=FALSE], pch = salinity_pch[sal],
           bg = season_cols[s], col = "white", cex = 1.6, lwd = 0.9)
}

arrows(0, 0, ax, ay, length = 0.08, lwd = 2.2, col = "black")
text(ax*1.12, ay*1.12, labels = rownames(sc_bp_cont), col = "black", cex = 0.8, font = 2)
abline(h = 0, v = 0, col = "grey60", lty = 3)

legend(x = max(xlim)*1.05, y = max(ylim), xpd = TRUE, title = "Season",
       legend = levels(meta$Season), pt.bg = season_cols, pch = 21,
       pt.cex = 1.3, col = "black", bty = "n", cex = 0.8)
legend(x = max(xlim)*1.05, y = max(ylim) - diff(ylim)*0.3, xpd = TRUE, title = "Salinity",
       legend = levels(meta$Salinity_class), pch = salinity_pch, pt.bg = "grey70",
       pt.cex = 1.3, col = "black", bty = "n", cex = 0.8)
dev.off()

cat("Figure4_dbRDA.pdf saved.\n")
