# =============================================================================
# Figure 1: Integrated UMAP Analysis
# =============================================================================
# Description: Generate Figure 1 panels showing integrated UMAP colored by 
#              celltype and genotype, dotplot of marker genes, violin plot of 
#              CNVs, and PCA of CNVs
# =============================================================================

# Load required libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)

# Set seed for reproducibility
set.seed(2021)

# Source configuration file
source("config.R")

# -----------------------------------------------------------------------------
# Check Required Input Files
# -----------------------------------------------------------------------------

required_files <- c(
  SEURAT_INTEGRATED,
  CNV_PCA_ALL
)

if (!check_input_files(required_files)) {
  stop("Missing required input files. Please check INPUT_FILES_LIST.md")
}

# -----------------------------------------------------------------------------
# Read Input Data
# -----------------------------------------------------------------------------

message("Loading input data...")

# Global Seurat object integrated with RISC
sobjint <- readRDS(SEURAT_INTEGRATED)

# PCA computed from InferCNV HMM states
pca_all <- readRDS(CNV_PCA_ALL)

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
# Define Color Scheme for Cell Types
# -----------------------------------------------------------------------------

colordf <- data.frame(
  celltype = factor(c(
    "Bcell", "Endothelial", "Macrophage", "Fibroblast", "Neutrophil", 
    "Osteoclast", "Plasma", "RBC", "Tcell", "Malignant", "Malignant-MSC", 
    "Malignant?", "Unclassified"
  )),
  color = c(
    'cyan1', 'darkolivegreen1', 'dodgerblue', 'forestgreen', 'darkviolet', 
    'lightslateblue', 'darkturquoise', 'yellow4', 'royalblue4', 'firebrick', 
    'red', 'pink', 'grey'
  )
)

sampcolordf <- colordf[match(levels(sobjint$IntCelltype), colordf$celltype), ]

# -----------------------------------------------------------------------------
# Figure 1B: Integrated UMAP Colored by Celltype
# -----------------------------------------------------------------------------

message("Generating Figure 1B: UMAP by celltype...")

sobjint$Celltype <- sobjint$IntCelltype

plot_1b <- DimPlot(sobjint, group.by = 'Celltype', label = TRUE, repel = TRUE) + 
  scale_color_manual(values = sampcolordf$color) + 
  theme_dimplot()

pdfname <- file.path(BASE_OUTPUT_DIR, '1B_umap_celltype.pdf')
pdf(pdfname, height = 4, width = 5)
print(plot_1b)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 1C: DotPlot of Marker Genes
# -----------------------------------------------------------------------------

message("Generating Figure 1C: Marker gene dotplot...")

# Marker genes selected based on Seurat FindMarkers + canonical markers
genes <- c(
  'Sp7', 'Sox9', 'Tpm2', 'Ncam1', 'Mki67', 'Top2a',
  'Ccl11', 'Mfap5', 'Rarres2', 'Acta2', 'Fap', 'Col1a1',
  'Pecam1', 'Cdh5', 'Gpihbp1',
  'Ptprc', 'Aif1', 'Ms4a6c', 'Cd68',
  'Nfatc1', 'Ocstamp', 'Ctsk',
  'Csf3r', 'Lcn2', 'S100a9',
  'Cd3e', 'Cd3d', 'Trbc2',
  'Cd19', 'Ms4a1', 'Igkc'
)

sobjint$Celltype <- factor(
  sobjint$IntCelltype, 
  levels = c(
    'Malignant', 'Fibroblast', 'Endothelial',
    'Macrophage', 'Osteoclast', 'Neutrophil',
    'Tcell', 'Bcell'
  )
)

plot_1c <- DotPlot(sobjint, features = rev(genes), group.by = 'Celltype') + 
  coord_flip() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 15),
    legend.title = element_text(size = 10)
  ) +
  ggtitle('Marker Genes') + 
  xlab('') + 
  ylab('')

pdfname <- file.path(BASE_OUTPUT_DIR, '1C_dotplot_markers.pdf')
pdf(pdfname, height = 6, width = 7)
print(plot_1c)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 1D: Integrated UMAP Colored by Genotype
# -----------------------------------------------------------------------------

message("Generating Figure 1D: UMAP by genotype...")

plot_1d <- DimPlot(
  sobjint, 
  group.by = 'Genotype', 
  label = FALSE, 
  repel = FALSE, 
  shuffle = TRUE
) + 
  theme_dimplot() + 
  scale_color_manual(values = c('steelblue', 'orange', 'firebrick'))

pdfname <- file.path(BASE_OUTPUT_DIR, '1D_umap_genotype.pdf')
pdf(pdfname, height = 4, width = 5)
print(plot_1d)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 1E: Violin Plot of Extreme CNVs
# -----------------------------------------------------------------------------

message("Generating Figure 1E: Violin plot of extreme CNVs...")

# Subset data (exclude Macrophage as negative control)
md <- sobjint@meta.data
md <- md[md$IntCelltype != 'Macrophage', ]
sobjx <- sobjint[, rownames(md)]

levs <- levels(md$IntCelltype)
levs <- levs[levs %in% md$IntCelltype]
sampcolordf_sub <- colordf[match(levs, colordf$celltype), ]

plot_1e <- VlnPlot(sobjx, 'ExtremeCNV_genes', pt.size = 0) + 
  NoLegend() + 
  scale_fill_manual(values = sampcolordf_sub$color) + 
  scale_y_continuous(
    breaks = c(0, 100, 250, 500, 1000, 2000, 3000), 
    minor_breaks = NULL
  ) +
  theme_linedraw() +
  xlab('') + 
  ylab('Number of genes\naffected by CNVs') +
  ggtitle('Number of Genes affected by\n>=2X deletion or amplification') +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.title = element_text(size = 14, face = 'bold'),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '1E_violplot_extremeCNVs.pdf')
pdf(pdfname, height = 4, width = 5)
print(plot_1e)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 1F: PCA of CNVs Colored by Celltypes
# -----------------------------------------------------------------------------

message("Generating Figure 1F: PCA of CNVs by celltype...")

# Exclude macrophages
md <- sobjint@meta.data
md <- md[md$IntCelltype != 'Macrophage', ]
sobjx <- sobjint[, rownames(md)]

levs <- levels(md$IntCelltype)
levs <- levs[levs %in% md$IntCelltype]
sampcolordf_sub <- colordf[match(levs, colordf$celltype), ]

# Create PCA dataframe
pcadf <- as.data.frame(pca_all$x[, 1:2])
pcadf <- pcadf[rownames(pcadf) %in% rownames(md), ]
pcadf$Celltype <- md$IntCelltype

plot_1f <- ggplot(pcadf, aes(PC1, PC2, col = Celltype)) +
  geom_point(size = 0.5) +
  theme_dimplot() +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = sampcolordf_sub$color) +
  ggtitle('PCA computed from inferred CNVs') +
  theme(
    axis.text = element_text(),
    panel.grid = element_line(color = 'grey85')
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '1F_CNV_pca_allcells_colby_celltype.pdf')
pdf(pdfname, height = 4, width = 5)
print(plot_1f)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 1 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
