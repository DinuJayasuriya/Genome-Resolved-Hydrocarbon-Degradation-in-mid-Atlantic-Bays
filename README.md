# Environmental Gradients Shape the Hydrocarbon-Degrading Microbiome in Two Mid-Atlantic Bays

Analysis code for the manuscript:

**Patabandige et al.** "Environmental Gradients Shape the Hydrocarbon-Degrading
Microbiome in Two Mid-Atlantic Bays." *mSystems* (in review).

Scripts to analyze hydrocarbon-degradation gene content and expression in
metagenome-assembled genomes (MAGs) from surface waters of the Chesapeake and
Delaware Bays. Each manuscript figure has its own script.

---

## Overview

The analyses combine two complementary approaches:

- **Community-level screening** of assembled metagenomic contigs and
  metatranscriptomic reads against a custom CANT-HYD DIAMOND database.
- **Genome-resolved (MAG-level) analysis** using DRAM/KEGG annotations linking
  MAGs to the hydrocarbon-degradation pathways they encode and express, and
  assessing functional redundancy and pathway co-expression.

---

## Scripts (by figure)

| Script | Produces | Description |
|--------|----------|-------------|
| `Figure1_MAG_relative_abundance.R` | Figure 1 | Relative abundance of top 20 MAGs (CPM), stacked bars by season |
| `Figure2_pathway_abundance_order.R` | Figure 2 | Mean HC-degradation pathway abundance by Order (log2 RPKG) |
| `Figure3_gene_abundance_heatmap.R` | Figure 3 | HC-degradation gene abundance within MAGs (log10 RPKG) |
| `Figure4_dbRDA.R` | Figure 4 | Distance-based redundancy analysis of HC-degrading MAGs |
| `Figure5_gene_environment_correlation.R` | Figure 5 | Spearman correlation matrix, genes vs environmental variables |
| `Figure6_expression_RPKGnormalized.R` | Figure 6 | Expression heatmap (normalized expression, row Z-score) |
| `Figure7_metabolic_network.R` | Figure 7 | Marker-gene presence/absence + MAG × category co-expression networks |

---

## Requirements

- R (>= 4.0)
- R packages: `dplyr`, `tidyr`, `stringr`, `pheatmap`, `RColorBrewer`,
  `paletteer`, `igraph`, `ggraph`, `readxl`, `vegan` (for dbRDA)

```r
install.packages(c("dplyr","tidyr","stringr","pheatmap","RColorBrewer",
                   "paletteer","igraph","ggraph","readxl","vegan"))
```

---

## Data availability

- **Sequence data:** NCBI Project PRJNA432171
- **Environmental / sampling metadata:** BCO-DMO,
  http://lod.bco-dmo.org/id/dataset/565451 (Campbell, 2016)
- **Custom DIAMOND database & archived code:** Jayasuriya D. (2026). Custom DIAMOND Protein Database for Hydrocarbon Degradation Genes in Estuarine Surface Water Metagenomes and Metatranscriptomes (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.20940772

---

## Citation

> Patabandige et al. Environmental Gradients Shape the Hydrocarbon-Degrading
> Microbiome in Two Mid-Atlantic Bays. mSystems (in review).

Archived code:

> [Author list]. [Year]. Analysis code for "Environmental Gradients Shape the
> Hydrocarbon-Degrading Microbiome in Two Mid-Atlantic Bays." Zenodo. DOI: [FILL IN]

---

## Contact

Dinuka Jayasuriya | Department of Biological Sciences, Clemson University
djayasu@g.clemson.edu
