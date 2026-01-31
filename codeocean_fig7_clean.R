# =============================================================================
# Figure 7: Pathology Subtype Analysis
# =============================================================================
# Description: Generate Figure 7 panels showing pathology subtype label
#              transfer UMAPs, module scores, marker dotplots, and
#              compositional analysis by pathology type
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
  PATHOLOGY_SEURAT,
  PATHOLOGY_ATLAS_MARKERS,
  PATHOLOGY_MALIGNANT_MARKERS,
  SEURAT_INTEGRATED,
  METADATA
)

if (!check_input_files(required_files)) {
  stop("Missing required input files. Please check INPUT_FILES_LIST.md")
}

# -----------------------------------------------------------------------------
# Read Input Data
# -----------------------------------------------------------------------------

message("Loading input data...")

# Pathology-annotated malignant object
sobjx <- readRDS(PATHOLOGY_SEURAT)

# Global integrated object
sobjint <- readRDS(SEURAT_INTEGRATED)

# Atlas and malignant markers
m <- read.csv(PATHOLOGY_ATLAS_MARKERS)
mm <- read.csv(PATHOLOGY_MALIGNANT_MARKERS)

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
# Figure 7A-C: Label Transfer DimPlots
# -----------------------------------------------------------------------------

message("Generating Figure 7A-C: Label transfer DimPlots...")

d1 <- DimPlot(sobjx, group.by = 'Code', label = TRUE, repel = TRUE, label.size = 2) + 
  theme_dimplot() +
  xlab('UMAP1') + 
  ylab('UMAP2') +
  ggtitle('Sample')

d2 <- DimPlot(sobjx, group.by = 'Bone_LabelTransfer_simple2') + 
  theme_dimplot() +
  scale_color_brewer(palette = 'Set2') +
  xlab('UMAP1') + 
  ylab('UMAP2') +
  ggtitle('Inferred Pathologic Subtype')

d3 <- FeaturePlot(sobjx, 'Bone_LabelTransfer_simple_score') + 
  theme_dimplot() +
  xlab('UMAP1') + 
  ylab('UMAP2') +
  ggtitle('Label Transfer Score')

dp <- (d1 + d2 + d3)

pdfname <- file.path(BASE_OUTPUT_DIR, '7A_7B_7C_LabelTransfer_DimPlots.pdf')
pdf(pdfname, height = 3, width = 12)
print(dp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 7D: Label Transfer Module Violin Plots
# -----------------------------------------------------------------------------

message("Generating Figure 7D: Module score violins...")

cts <- levels(sobjx$Bone_LabelTransfer_simple2)
vpl <- lapply(cts, function(ct) {
  md <- sobjx@meta.data
  indf <- data.frame(
    CelltypeCall = md$Bone_LabelTransfer_simple2,
    ModuleScore = md[, ct],
    Code = md$Code,
    Genotype = md$Genotype
  )
  
  ctlab <- gsub('-like', '', ct)
  
  d1 <- ggplot(indf, aes(CelltypeCall, ModuleScore, fill = CelltypeCall)) +
    geom_violin() +
    geom_jitter(size = 0.00001, alpha = 0.01) +
    scale_fill_brewer(palette = 'Set2') +
    ggtitle(paste0(ctlab, ' module')) +
    xlab('Genotype') + 
    ylab('Module Score') +
    theme_light() +
    scale_y_continuous(limits = c(-4, 4)) +
    NoLegend() +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y = element_text(size = 6),
      axis.title.x = element_blank(),
      legend.title = element_text(size = 10),
      plot.margin = unit(c(0.2, 0.2, 0, 0), "cm")
    )
  
  return(d1)
})

vpp <- wrap_plots(vpl)

pdfname <- file.path(BASE_OUTPUT_DIR, '7D_LabelTransfer_VlnModules.pdf')
pdf(pdfname, height = 4, width = 6)
print(vpp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 7E: Label Transfer Marker Dotplot
# -----------------------------------------------------------------------------

message("Generating Figure 7E: Marker dotplot...")

# Top markers
n <- 10
top <- mm %>% group_by(cluster) %>% top_n(n = n, wt = score) %>% as.data.frame()
topsmall <- top[, c('cluster', 'gene')]

# Add canonical markers
addgenes <- data.frame(
  cluster = c('Osteo-like', 'Osteo-like', 'Fibro-like', 'Chondro-like', 
              'Chondro-like', 'Pericyte-like', 'MSC-like', 'Endothelial-like'),
  gene = c('Col1a1', 'Runx2', 'S100a4', 'Col2a1', 'Sox9', 'Rgs5', 'Lepr', 'Pecam1')
)

topsmall[topsmall$gene %in% addgenes$gene, ]
topsmall <- topsmall[!(topsmall$gene %in% addgenes$gene), ]
topsmall <- rbind(addgenes, topsmall)

cts <- levels(sobjx$Bone_LabelTransfer_simple2)
topsmall$cluster <- factor(topsmall$cluster, levels = cts)
topsmall <- topsmall[order(topsmall$cluster), ]

# Add lines between celltypes
xintlines <- cumsum(rev(table(topsmall$cluster))) + 0.5
xintlines <- xintlines[1:length(xintlines) - 1]

dp <- DotPlot(sobjx, features = rev(topsmall$gene), group.by = 'Bone_LabelTransfer_simple2') + 
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  geom_vline(xintercept = xintlines, linewidth = 0.1)

pdfname <- file.path(BASE_OUTPUT_DIR, '7E_LabelTransfer_DotPlot_malignantonly.pdf')
pdf(pdfname, height = 6, width = 5)
print(dp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 7F: Code Composition Barplot
# -----------------------------------------------------------------------------

message("Generating Figure 7F: Code barplot...")

md <- sobjx@meta.data
twt <- as.data.frame.matrix(table(md$Bone_LabelTransfer_simple2, md$Code))
twt <- cbind(rownames(twt), twt)
colnames(twt)[1] <- 'Bone_LabelTransfer_simple2'

twtl <- as.data.frame(tidyr::pivot_longer(twt, cols = colnames(twt)[-1]))
colnames(twtl)[2] <- 'Code'
twtl$Genotype <- str_split_fixed(twtl$Code, '_', n = 2)[, 1]
twtl$Code <- factor(twtl$Code, levels = levels(sobjx$Code))
twtl$Genotype <- factor(twtl$Genotype, levels = levels(sobjx$Genotype))
twtl$Bone_LabelTransfer_simple2 <- factor(twtl$Bone_LabelTransfer_simple2, 
                                          levels = levels(sobjx$Bone_LabelTransfer_simple2))

gbar <- ggplot(twtl, aes(Code, y = value, fill = Bone_LabelTransfer_simple2)) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_linedraw() +
  labs(fill = 'Inferred Pathology') +
  scale_fill_brewer(palette = 'Set2') +
  ylab('Proportion of Cell Types') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '7F_LabelTransfer_CodeBarPlot.pdf')
pdf(pdfname, height = 3, width = 5)
print(gbar)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 7G: Whole UMAP Split by Pathology
# -----------------------------------------------------------------------------

message("Generating Figure 7G: Whole UMAP split by pathology...")

# Add pathology to integrated object
pbmd <- as.data.frame(readxl::read_excel(METADATA))
pbmd <- pbmd[match(levels(sobjint$Code), pbmd$Code), ]
pbmd$LabelTransferPathology <- paste0(pbmd$LabelTransferPathology, '-like')

sobjint$InferredPathology <- sobjint$Code
sobjint$InferredPathology <- plyr::mapvalues(
  sobjint$InferredPathology,
  from = levels(sobjint$Code),
  to = pbmd$LabelTransferPathology
)

# Color scheme
colordf <- data.frame(
  celltype = factor(c("Bcell", "Endothelial", "Macrophage", "Fibroblast", 
                     "Neutrophil", "Osteoclast", "Plasma", "RBC", "Tcell", 
                     "Malignant", 'Malignant-MSC', 'Malignant?', "Unclassified")),
  color = c('cyan1', 'darkolivegreen1', 'dodgerblue', 'forestgreen', 
            'darkviolet', 'lightslateblue', 'darkturquoise', 'yellow4', 
            'royalblue4', 'firebrick', 'red', 'pink', 'grey')
)
sampcolordf <- colordf[match(levels(sobjint$IntCelltype), colordf$celltype), ]

sobjint$Celltype <- sobjint$IntCelltype
d <- DimPlot(sobjint, group.by = 'Celltype', split.by = 'InferredPathology', ncol = 3) + 
  scale_color_manual(values = sampcolordf$color) + 
  theme_dimplot() +
  xlab('UMAP 1') + 
  ylab('UMAP 2') +
  theme(
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, 'cm'),
    axis.title = element_text(size = 6),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '7G_LabelTransfer_WHoleUMAP_Split.pdf')
pdf(pdfname, height = 2, width = 6)
print(d)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 7H: Whole Data Composition Barplot
# -----------------------------------------------------------------------------

message("Generating Figure 7H: Whole data barplot...")

twt <- as.data.frame.matrix(table(sobjint$InferredPathology, sobjint$Celltype))
twt <- cbind(rownames(twt), twt)
colnames(twt)[1] <- 'InferredPathology'

twtl <- as.data.frame(tidyr::pivot_longer(twt, cols = colnames(twt)[-1]))
colnames(twtl)[2] <- 'Celltype'

twtl$InferredPathology <- factor(twtl$InferredPathology, 
                                 levels = levels(sobjint$InferredPathology))
twtl$Celltype <- factor(twtl$Celltype, levels = levels(sobjint$Celltype))

gbar <- ggplot(twtl, aes(InferredPathology, value, fill = Celltype)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual(values = sampcolordf$color) +
  theme_linedraw() +
  labs(fill = 'Celltype') +
  ylab('Proportion of Cell Types') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '7H_LabelTransfer_WholeData_Bar.pdf')
pdf(pdfname, height = 4, width = 3)
print(gbar)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 7I: Cytokine Expression by Pathology
# -----------------------------------------------------------------------------

message("Generating Figure 7I: Cytokine dotplot...")

pathos <- c("Osteo-like", "Fibro-like", "Chondro-like")
md <- sobjx@meta.data
md <- md[md$Bone_LabelTransfer_simple2 %in% c(pathos), ]
sobjxxx <- sobjx[, rownames(md)]
sobjxxx <- SetIdent(sobjxxx, value = sobjxxx$Bone_LabelTransfer_simple2)

genes <- c('Ccl2', 'Ccl7', 'Csf1', 'Tnfsf11')

dp <- DotPlot(sobjxxx, rev(genes)) + 
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '7I_LabelTransfer_cytokines.pdf')
pdf(pdfname, height = 4, width = 3)
print(dp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 7 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
