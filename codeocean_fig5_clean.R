# =============================================================================
# Figure 5: CNV Analysis
# =============================================================================
# Description: Generate Figure 5 panels showing CNV burden violin plots, CNV
#              PCA colored by genotype and code, CNV-DEG Venn diagrams, and
#              differential CNV heatmaps
# =============================================================================

# Load required libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(ggVennDiagram)
library(progress)

# Set seed for reproducibility
set.seed(2021)

# Source configuration file
source("config.R")

# -----------------------------------------------------------------------------
# Check Required Input Files
# -----------------------------------------------------------------------------

required_files <- c(
  SEURAT_INTEGRATED,
  CNV_MATRIX,
  CNV_PCA_ALL,
  CNV_PCA_MALIGNANT,
  CNV_CROSSCOND_MEAN,
  CNV_CROSSCOND_AMPS,
  DE_TKO_DKO_MALIGNANT,
  DE_DKOAA_DKO_MALIGNANT,
  METADATA
)

if (!check_input_files(required_files)) {
  stop("Missing required input files. Please check INPUT_FILES_LIST.md")
}

# -----------------------------------------------------------------------------
# Read Input Data
# -----------------------------------------------------------------------------

message("Loading input data...")

# Global Seurat object
sobjint <- readRDS(SEURAT_INTEGRATED)

# CNV data
hmm <- readRDS(CNV_MATRIX)
pca_all <- readRDS(CNV_PCA_ALL)
pca_m <- readRDS(CNV_PCA_MALIGNANT)

# Cross-condition CNV comparisons
mcnv <- readRDS(CNV_CROSSCOND_MEAN)
mcountcnv <- readRDS(CNV_CROSSCOND_AMPS)

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
# Calculate CNV Metrics
# -----------------------------------------------------------------------------

message("Calculating CNV metrics...")

# Calculate extreme deletions
pb <- progress::progress_bar$new(total = ncol(hmm), width = 50)
numCNV_extremedels <- sapply(1:ncol(hmm), FUN = function(i) {
  pb$tick()
  x <- hmm[, i]
  x <- factor(x, levels = c(0, 0.5, 1, 1.5, 2, 3))
  tab <- table(x)
  tab <- tab[names(tab) == 0]
  sum(tab)
})
names(numCNV_extremedels) <- colnames(hmm)

md <- sobjint@meta.data
md$ExtremeDEL_genes <- NA
md[match(names(numCNV_extremedels), rownames(md)), "ExtremeDEL_genes"] <- numCNV_extremedels
sobjint$ExtremeDEL_genes <- md$ExtremeDEL_genes

# Calculate extreme amplifications
pb <- progress::progress_bar$new(total = ncol(hmm), width = 50)
numCNV_extremeamps <- sapply(1:ncol(hmm), FUN = function(i) {
  pb$tick()
  x <- hmm[, i]
  x <- factor(x, levels = c(0, 0.5, 1, 1.5, 2, 3))
  tab <- table(x)
  tab <- tab[names(tab) %in% c(2, 3)]
  sum(tab)
})
names(numCNV_extremeamps) <- colnames(hmm)

md <- sobjint@meta.data
md$ExtremeAMP_genes <- NA
md[match(names(numCNV_extremeamps), rownames(md)), "ExtremeAMP_genes"] <- numCNV_extremeamps
sobjint$ExtremeAMP_genes <- md$ExtremeAMP_genes

# -----------------------------------------------------------------------------
# Figure 5B: Violin Plots of CNVs in Malignant Cells
# -----------------------------------------------------------------------------

message("Generating Figure 5B: CNV burden violin plots...")

# Subset to malignant only
md <- sobjint@meta.data
md <- md[md$IntCelltype == 'Malignant', ]
sobjx <- sobjint[, rownames(md)]
md <- sobjx@meta.data

plot.title.size <- 18
axis.text.size <- 15
y.axis.text.size <- 10

vln_cnv_ext <- ggplot(md, aes(Genotype, ExtremeCNV_genes, fill = Genotype)) +
  geom_violin() +
  theme_linedraw() +
  scale_fill_manual(values = c('steelblue', 'orange', 'firebrick')) +
  scale_y_continuous(breaks = c(0, 250, 500, 1000, 2000, 3000), 
                     minor_breaks = NULL, limits = c(0, 3500)) +
  xlab('') + 
  ylab('Extreme CNV Genes') +
  ggtitle('Extreme CNVs: Amps and Dels') +
  theme(
    axis.text.x = element_text(size = axis.text.size),
    axis.text.y = element_text(size = y.axis.text.size),
    plot.title = element_text(size = plot.title.size, face = 'bold')
  ) +
  NoLegend()

vln_cnv_amp <- ggplot(md, aes(Genotype, ExtremeAMP_genes, fill = Genotype)) +
  geom_violin() +
  theme_linedraw() +
  scale_fill_manual(values = c('steelblue', 'orange', 'firebrick')) +
  scale_y_continuous(breaks = c(0, 250, 500, 1000, 2000, 3000), 
                     minor_breaks = NULL, limits = c(0, 3500)) +
  xlab('') + 
  ylab('Extreme CNV Genes') +
  ggtitle('Extreme CNVs: Amplifications (>2X) only') +
  theme(
    axis.text.x = element_text(size = axis.text.size),
    axis.text.y = element_text(size = y.axis.text.size),
    plot.title = element_text(size = plot.title.size, face = 'bold')
  ) +
  NoLegend()

vln_cnv_del <- ggplot(md, aes(Genotype, ExtremeDEL_genes, fill = Genotype)) +
  geom_violin() +
  theme_linedraw() +
  scale_fill_manual(values = c('steelblue', 'orange', 'firebrick')) +
  scale_y_continuous(breaks = c(0, 250, 500, 1000, 2000, 3000), 
                     minor_breaks = NULL, limits = c(0, 3500)) +
  xlab('') + 
  ylab('Extreme CNV Genes') +
  ggtitle('Extreme CNVs: Deletions (both copies) only') +
  theme(
    axis.text.x = element_text(size = axis.text.size),
    axis.text.y = element_text(size = y.axis.text.size),
    plot.title = element_text(size = plot.title.size, face = 'bold')
  ) +
  NoLegend()

cnv_vl <- (vln_cnv_ext / vln_cnv_amp / vln_cnv_del)

pdfname <- file.path(BASE_OUTPUT_DIR, '5B_VLN_compare_CNV_counts.pdf')
pdf(pdfname, height = 8, width = 6)
print(cnv_vl)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 5C: PCA of Malignant Cells Colored by Genotype
# -----------------------------------------------------------------------------

message("Generating Figure 5C: CNV PCA by genotype...")

md <- sobjint@meta.data
md <- md[md$IntCelltype == 'Malignant', ]

pcadf <- as.data.frame(pca_m$x[, 1:2])
pcadf <- pcadf[rownames(pcadf) %in% rownames(md), ]
pcadf$Genotype <- md$Genotype
pcadf$Code <- md$Code

d <- ggplot(pcadf, aes(PC1, PC2, col = Genotype)) +
  geom_point(size = 0.5) +
  theme_dimplot() +
  scale_color_manual(values = c('steelblue', 'orange', 'firebrick')) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  ggtitle('PCA computed from inferred CNVs') +
  theme(
    axis.text = element_text(),
    panel.grid = element_line(color = 'grey85')
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '5C_CNV_pca_malignantcells_colby_genotype.pdf')
pdf(pdfname, height = 4, width = 5)
print(d + NoLegend())
print(d)
dev.off()
message("  Saved: ", pdfname)

# PCA split by genotype, colored by code
colpal <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
            "#E5C494", "#FB8072", "#8DA0CB", "#BC80BD")

d2 <- ggplot(pcadf, aes(PC1, PC2, col = Code)) +
  facet_wrap(~Genotype) +
  geom_point(size = 0.5) +
  theme_dimplot() +
  scale_color_manual(values = colpal) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(
    axis.text = element_text(),
    panel.grid = element_line(color = 'grey85')
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '5C_CNV_pca_malignantcells_splitbygenotype_colby_code.pdf')
pdf(pdfname, height = 2.5, width = 5)
print(d2 + NoLegend())
print(d2)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 5E: Venn Diagram of Differential CNVs and DEGs
# -----------------------------------------------------------------------------

message("Generating Figure 5E: CNV-DEG Venn diagrams...")

# TKO vs DKO
diffcnv <- mcnv$TKO_vs_DKO
diffcnv <- diffcnv[diffcnv$p < 0.1 & diffcnv$meandiff > 0, ]

diffcountcnv <- mcountcnv$TKO_vs_DKO
diffcountcnv <- diffcountcnv[diffcountcnv$p < 0.1 & diffcountcnv$meandiff > 0, ]

deres <- read.csv(DE_TKO_DKO_MALIGNANT)
degs <- deres[deres$FDR < 0.1 & deres$logFC > 0 & deres$pct.1 > 0.1, ]
colnames(degs)[1] <- 'gene'

inlist_rows <- list(
  DiffAmps_mean = diffcnv,
  DEGs = degs,
  DiffAmps_count = diffcountcnv
)
inlist <- lapply(inlist_rows, function(x) { x[, 'gene'] })

ggvenn_TKOvDKO <- ggVennDiagram::ggVennDiagram(inlist) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  NoLegend()

int_tko_dko <- Reduce(intersect, inlist)

# Get union of intersects
int1 <- intersect(inlist$DiffAmps_mean, inlist$DEGs)
int2 <- intersect(inlist$DiffAmps_count, inlist$DEGs)
tkodko_lenient_int <- union(int1, int2)

# DKOAA vs DKO
diffcnv <- mcnv$DKOAA_vs_DKO
diffcnv <- diffcnv[diffcnv$p < 0.1 & diffcnv$meandiff > 0, ]

diffcountcnv <- mcountcnv$DKOAA_vs_DKO
diffcountcnv <- diffcountcnv[diffcountcnv$p < 0.1 & diffcountcnv$meandiff > 0, ]

deres <- read.csv(DE_DKOAA_DKO_MALIGNANT)
degs <- deres[deres$FDR < 0.1 & deres$logFC > 0 & deres$pct.1 > 0.1, ]
colnames(degs)[1] <- 'gene'

inlist_rows <- list(
  DiffAmps_mean = diffcnv,
  DEGs = degs,
  DiffAmps_count = diffcountcnv
)
inlist <- lapply(inlist_rows, function(x) { x[, 'gene'] })

ggvenn_DKOAAvDKO <- ggVennDiagram::ggVennDiagram(inlist) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  NoLegend()

int_dkoaa_dko <- intersect(inlist$DEGs, inlist$DiffAmps_count)

pdfname <- file.path(BASE_OUTPUT_DIR, '5E_CNV-DEG_Venn.pdf')
pdf(pdfname, height = 4, width = 5)
print(ggvenn_TKOvDKO)
print(ggvenn_DKOAAvDKO)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 5G: Heatmap of Differential CNVs
# -----------------------------------------------------------------------------

message("Generating Figure 5G: Differential CNV heatmap...")

# Union of overlapping genes
intgenes <- c(tkodko_lenient_int, 'Pa2g4')

# Mean heatmap
cnvmat <- mcnv$TKO_vs_DKO
cnvmat <- cnvmat[cnvmat$gene %in% intgenes, ]

cnvmatdkoaa <- mcnv$DKOAA_vs_DKO
cnvmatdkoaa <- cnvmatdkoaa[cnvmatdkoaa$gene %in% intgenes, ]

cnvmat2 <- cnvmat[, grepl('DKO', colnames(cnvmat))]
cnvmat2 <- cbind(cnvmat2, cnvmatdkoaa[, grepl('DKOAA', colnames(cnvmatdkoaa))])
cnvmat2 <- cbind(cnvmat2, cnvmat[, grepl('TKO', colnames(cnvmat))])
cnvmat <- cnvmat2

col_fun <- circlize::colorRamp2(c(0, 1, 2), c("blue", "white", "red"))
column_split <- factor(c(rep('DKO', 3), rep('DKOAA', 3), rep('TKO', 4)), 
                      levels = c('DKO', 'DKOAA', 'TKO'))

cnvhm <- ComplexHeatmap::Heatmap(
  cnvmat, 
  show_row_names = TRUE, 
  name = "Average CNV state\nper sample\n(0 = deleted,\n1 = no CNV,\n2 = 2x amplified)", 
  rect_gp = grid::gpar(col = "black", lwd = 0.1),
  border_gp = grid::gpar(col = "black", lwd = 1),
  column_split = column_split,
  cluster_columns = FALSE,
  col = col_fun
)

pdfname <- file.path(BASE_OUTPUT_DIR, '4G_DiffCNVHeatmap_uniongenes.pdf')
pdf(pdfname, height = 5, width = 6)
print(cnvhm)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 5 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
