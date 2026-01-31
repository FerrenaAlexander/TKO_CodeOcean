# =============================================================================
# Figure 6: Human Osteosarcoma Comparison
# =============================================================================
# Description: Generate Figure 6 panels showing human OS UMAPs, CNV violin
#              plot, marker dotplot, and module score heatmap comparing human
#              and murine datasets
# =============================================================================

# Load required libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(grid)

# Set seed for reproducibility
set.seed(2021)

# Source configuration file
source("config.R")

# -----------------------------------------------------------------------------
# Check Required Input Files
# -----------------------------------------------------------------------------

required_files <- c(
  HUMAN_OS_SEURAT,
  HUMAN_OS_MARKERS,
  SEURAT_INTEGRATED
)

if (!check_input_files(required_files)) {
  stop("Missing required input files. Please check INPUT_FILES_LIST.md")
}

# -----------------------------------------------------------------------------
# Read Input Data
# -----------------------------------------------------------------------------

message("Loading input data...")

# Human OS dataset with murine module scores
hum <- readRDS(HUMAN_OS_SEURAT)

# Human OS markers
pm <- read.csv(HUMAN_OS_MARKERS)

# Mouse integrated object (for levels)
sobjint <- readRDS(SEURAT_INTEGRATED)

message("Input data loaded successfully!")

# -----------------------------------------------------------------------------
# Define Custom Functions
# -----------------------------------------------------------------------------

theme_dimplot <- function(titlesize = 15) {
  base_size <- 11
  base_family <- ""
  base_line_size <- base_size / 22
  base_rect_size <- base_size / 22
  
  theme_grey(
    base_size = base_size, 
    base_family = base_family,
    base_line_size = base_line_size, 
    base_rect_size = base_rect_size
  ) %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = "black"), 
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey85", colour = "grey20"),
      legend.key = element_rect(fill = "white", colour = NA), 
      complete = TRUE,
      axis.ticks = element_blank(), 
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', 
                                size = titlesize)
    )
}

# -----------------------------------------------------------------------------
# Figure 6A: Human OS UMAPs
# -----------------------------------------------------------------------------

message("Generating Figure 6A: Human OS UMAPs...")

d1 <- DimPlot(hum, group.by = 'Condition') + 
  theme_dimplot() + 
  ggtitle('Human OS datasets')

d2 <- DimPlot(hum, group.by = 'CelltypeFinal', label.size = 2, label = TRUE, repel = TRUE) + 
  theme_dimplot() + 
  ggtitle('Re-annotated Celltypes')

pdfname <- file.path(BASE_OUTPUT_DIR, '6A_HumanOS_UMAPs.pdf')
pdf(pdfname, height = 3, width = 8)
print(d1 + d2)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 6B: CNV Violin and Marker Dotplot
# -----------------------------------------------------------------------------

message("Generating Figure 6B: CNV violin and marker dotplot...")

# Reorder cell types to match murine
hum$CelltypeFinal_reord <- factor(
  hum$CelltypeFinal,
  levels = c('Malignant', 'Fibroblast', 'Pericyte', 'Endothelial',
             'Macrophage', 'Osteoclast', 'Tcell', 'Bcell')
)

md <- hum@meta.data
axis.text.size <- 8
y.axis.text.size <- axis.text.size
plot.title.size <- 12

# CNV violin plot
vln_cnv_ext <- ggplot(md, aes(CelltypeFinal_reord, ExtremeCNV_genes, fill = CelltypeFinal_reord)) +
  geom_violin(scale = 'width') +
  theme_linedraw() +
  scale_y_continuous(breaks = c(0, 250, 500, 1000, 2000, 3000), 
                     minor_breaks = NULL, limits = c(0, 3500)) +
  xlab('') + 
  ylab('Extreme CNV Genes') +
  ggtitle('Number of Genes affected by Extreme CNVs') +
  theme(
    axis.text.x = element_text(size = axis.text.size, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = y.axis.text.size),
    plot.title = element_text(size = plot.title.size, face = 'bold')
  ) +
  NoLegend()

# Prepare markers
pm$score <- (pm$pct.1 - pm$pct.2) * pm$avg_log2FC
pm <- pm[pm$avg_log2FC > 0.5, ]
pm <- pm[pm$pct.1 > 0.5, ]
pm <- pm[pm$pct.2 < 0.5, ]
pm <- pm[pm$p_val_adj < 0.05, ]
pm$cluster <- factor(pm$cluster, levels = unique(pm$cluster))
pm <- pm[order(pm$score, decreasing = TRUE), ]
pm <- pm[!duplicated(pm$gene), ]
pm <- pm[order(pm$cluster), ]

# Manually selected genes
gene <- c(
  'SP7', 'SOX9', 'TPM2', 'NCAM1', 'MKI67', 'TOP2A',  # Malignant from mouse
  "CPE", "PTH1R", "MDFI",  # Malignant top markers
  "C1S", "EFEMP1", "MFAP5",  # Fibroblast top markers
  'RARRES2', 'ACTA2', 'FAP', 'COL1A1',  # Fibroblast mouse markers
  'TAGLN', 'THY1', 'PLAC9', 'RGS5', 'NDUFA4L2', 'NOTCH3',  # Pericyte
  "PLVAP", "VWF", "CLEC14A", 'PECAM1',  # Endothelial
  'PTPRC',  # Pan immune
  'MS4A6A', 'MS4A7', 'AIF1', 'CD68',  # Macrophage
  'NFATC1', 'OCSTAMP', 'ATP6V0D2', 'CTSK',  # Osteoclast
  'CD3D', 'TRBC2', 'CD2',  # T cell
  'IGHG1', 'JCHAIN', 'MZB1'  # B cell
)

y.axis.text.size <- 8

dp <- DotPlot(hum, features = rev(gene), group.by = "CelltypeFinal_reord") + 
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = axis.text.size),
    axis.text.y = element_text(size = y.axis.text.size),
    axis.title = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '6B_HumanOS_CNVviolin_and_MakrerDotPlot.pdf')
pdf(pdfname, height = 8, width = 8)
print(wrap_plots(list(vln_cnv_ext, dp), ncol = 1, heights = c(0.3, 0.7)))
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 6C: Module Score Heatmap
# -----------------------------------------------------------------------------

message("Generating Figure 6C: Module score heatmap...")

# Reorder to match murine
hum$CelltypeFinal_reord <- factor(
  hum$CelltypeFinal,
  levels = c('Malignant', 'Fibroblast', 'Pericyte', 'Endothelial',
             'Macrophage', 'Osteoclast', 'Tcell', 'Bcell')
)

cts <- c('Malignant', 'Fibroblast', 'Endothelial', 'Macrophage', 
         'Osteoclast', 'Neutrophil', 'Tcell', 'Bcell')

# Average heatmap
md <- hum@meta.data
gdf <- md[, c('CelltypeFinal_reord', cts)]
gdf <- aggregate(. ~ CelltypeFinal_reord, gdf, median)
rownames(gdf) <- gdf[, 1]
gdf <- gdf[, -1]
gdf <- t(gdf)

rownames(gdf) <- paste0('Murine ', rownames(gdf), ' Score')

hm <- ComplexHeatmap::Heatmap(
  gdf,
  column_names_rot = 45,
  name = 'Module Score',
  rect_gp = gpar(col = "white", lwd = 0.2),
  border_gp = gpar(col = "black", lwd = 1),
  column_title = 'Human OS Cell Types',
  column_title_side = 'bottom',
  row_title = 'Murine Signature Scores',
  row_title_side = 'right'
)

pdfname <- file.path(BASE_OUTPUT_DIR, '6C_HumanOS_ModuleScoreHeatmap.pdf')
pdf(pdfname, height = 6, width = 7)
draw(hm, heatmap_legend_side = 'left')
dev.off()
message("  Saved: ", pdfname)

# Heatmap without pericytes
md <- hum@meta.data
md <- md[md$CelltypeFinal != 'Pericyte', ]
hum_nop <- hum[, rownames(md)]

cts <- c('Malignant', 'Fibroblast', 'Endothelial', 'Macrophage', 'Osteoclast', 'Tcell', 'Bcell')

md <- hum_nop@meta.data
gdf <- md[, c('CelltypeFinal_reord', cts)]
gdf <- aggregate(. ~ CelltypeFinal_reord, gdf, median)
rownames(gdf) <- gdf[, 1]
gdf <- gdf[, -1]
gdf <- t(gdf)

rownames(gdf) <- paste0('Murine ', rownames(gdf), ' Score')

hm2 <- ComplexHeatmap::Heatmap(
  gdf,
  column_names_rot = 45,
  name = 'Module Score',
  rect_gp = gpar(col = "white", lwd = 0.2),
  border_gp = gpar(col = "black", lwd = 1),
  column_title = 'Human OS Cell Types',
  column_title_side = 'bottom',
  row_title = 'Murine Signature Scores',
  row_title_side = 'right'
)

pdfname <- file.path(BASE_OUTPUT_DIR, '6C_HumanOS_ModuleScoreHeatmap_IntersectingCelltypes.pdf')
pdf(pdfname, height = 6, width = 7)
draw(hm, heatmap_legend_side = 'right')
draw(hm2, heatmap_legend_side = 'right')
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 6 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
