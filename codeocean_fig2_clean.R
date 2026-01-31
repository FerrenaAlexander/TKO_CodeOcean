# =============================================================================
# Figure 2: GSEA and EMT Analysis
# =============================================================================
# Description: Generate Figure 2 panels showing GSEA dotplots for TKO vs DKO
#              and DKOAA vs DKO comparisons, EMT/motility gene dotplot, and
#              metastasis-related gene signatures
# =============================================================================

# Load required libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(fgsea)

# Set seed for reproducibility
set.seed(2021)

# Source configuration file
source("config.R")

# -----------------------------------------------------------------------------
# Check Required Input Files
# -----------------------------------------------------------------------------

required_files <- c(
  SEURAT_CELLSUBTYPES,
  PATHWAY_TKO_DKO,
  PATHWAY_DKOAA_DKO
)

if (!check_input_files(required_files)) {
  stop("Missing required input files. Please check INPUT_FILES_LIST.md")
}

# -----------------------------------------------------------------------------
# Read Input Data
# -----------------------------------------------------------------------------

message("Loading input data...")

# Malignant subsetted object
sobjx <- readRDS(SEURAT_CELLSUBTYPES)$Malignant

# GSEA results TKO/DKO
res <- read.csv(PATHWAY_TKO_DKO)

# GSEA results DKOAA/DKO
dkoaa_dko_res <- read.csv(PATHWAY_DKOAA_DKO)

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
# Figure 2A and 2B: TKO vs DKO and DKOAA vs DKO Malignant Cell GSEA Dotplot
# -----------------------------------------------------------------------------

message("Generating Figure 2A and 2B: GSEA dotplots...")

# Remove extra column
res <- res[, 1:(ncol(res) - 1)]
dkoaa_dko_res <- dkoaa_dko_res[, 1:(ncol(dkoaa_dko_res) - 1)]

# Threshold by padj
res <- res[res$padj < 0.1, ]
dkoaa_dko_res <- dkoaa_dko_res[dkoaa_dko_res$padj < 0.1, ]

# Order by NES * pval
res$ord <- sign(res$NES) * -log10(res$pval)
res <- res[order(res$ord, decreasing = TRUE), ]

dkoaa_dko_res$ord <- sign(dkoaa_dko_res$NES) * -log10(dkoaa_dko_res$pval)
dkoaa_dko_res <- dkoaa_dko_res[order(dkoaa_dko_res$ord, decreasing = TRUE), ]

# Get top significant pathways
res <- rbind(head(res, 20), tail(res, 20))
dkoaa_dko_res <- rbind(head(dkoaa_dko_res, 20), tail(dkoaa_dko_res, 20))

# Clean pathway names
res$pathway <- gsub('_', ' ', res$pathway)
res$pathway <- gsub('HALLMARK ', '', res$pathway)
res$pathway <- factor(res$pathway, levels = unique(rev(res$pathway)))

# Create TKO vs DKO plot
d <- ggplot(res, aes(x = -log10(pval), y = pathway, size = size, col = NES)) +
  geom_point() +
  scale_color_gradient2(
    low = 'steelblue', 
    high = 'red', 
    mid = 'white', 
    midpoint = 0, 
    name = 'Normalized\nEnrichment\nScore'
  ) +
  theme_linedraw() +
  scale_size(range = c(2, 6)) +
  ylab('') + 
  labs(title = 'Hallmarks GSEA TKO vs DKO') +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 15),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

# Clean DKOAA pathway names
dkoaa_dko_res$pathway <- gsub('_', ' ', dkoaa_dko_res$pathway)
dkoaa_dko_res$pathway <- gsub('HALLMARK ', '', dkoaa_dko_res$pathway)
dkoaa_dko_res$pathway <- factor(
  dkoaa_dko_res$pathway, 
  levels = unique(rev(dkoaa_dko_res$pathway))
)

# Create DKOAA vs DKO plot
d2 <- ggplot(dkoaa_dko_res, aes(x = -log10(pval), y = pathway, size = size, col = NES)) +
  geom_point() +
  scale_color_gradient2(
    low = 'steelblue', 
    high = 'orange', 
    mid = 'white', 
    midpoint = 0, 
    name = 'Normalized\nEnrichment\nScore'
  ) +
  theme_linedraw() +
  scale_size(range = c(2, 6)) +
  ylab('') + 
  labs(title = 'Hallmarks GSEA DKOAA vs DKO') +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 15),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

plot_2a_2b <- d / d2

pdfname <- file.path(BASE_OUTPUT_DIR, '2A_2B_DE_malignant_GSEA.pdf')
pdf(pdfname, height = 12, width = 4.5)
print(plot_2a_2b)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 2C: Dotplot of EMT/Motility Related Genes
# -----------------------------------------------------------------------------

message("Generating Figure 2C: EMT/motility gene dotplot...")

genes <- c(
  'Mmp2', 'Mmp9', 'Mmp14', 'Mmp23', 'Pcolce', 
  'Timp1', 'Loxl2', 'Postn', 'Fap', 'Sparc',
  'Col1a1', 'Col2a1', 'Col3a1', 'Col4a1',
  'Rac1', 'Rac3', 'Rhoa', 'Rhoq', 'Cdc42bpa', 'Enah',
  'Cdh1', 'Epcam', 'Cdh2', 'Vim',
  'Zeb1', 'Zeb2', 'Twist1', 'Twist2', 'Snai1', 'Snai2'
)

plot_2c <- DotPlot(sobjx, features = rev(genes), group.by = 'Genotype') + 
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '2C_Malignant_MotilityEMT_genes_DOPTPLOT.pdf')
pdf(pdfname, height = 5, width = 4)
print(plot_2c)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 2D/2E: Metastasis Signature Module Scores
# -----------------------------------------------------------------------------

message("Generating Figure 2D/2E: Metastasis gene module scores...")

# Define metastasis signatures from literature
# Qin et al., Aging 2023: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10415573/
metsigs <- list(
  Qin_AgingAlbany2023_ProMetastatic = c(
    'Sema5a', 'Hilpda', 'Wif1', 'Alpl', 'Lox', 'Sec61g'
  ),
  Qin_AgingAlbany2023_AntiMetastatic = c(
    'Aspa', 'Rbm24', 'Gpr174', 'Nlgn1', 'Igf2', 'Pdlim3'
  )
)

# Calculate module scores
mdcols <- ncol(sobjx@meta.data)
sobjx <- AddModuleScore(sobjx, features = metsigs, name = names(metsigs))

mods <- sobjx@meta.data[, (mdcols + 1):ncol(sobjx@meta.data)]
sobjx@meta.data <- sobjx@meta.data[, 1:mdcols]
colnames(mods) <- names(metsigs)
sobjx@meta.data <- cbind(sobjx@meta.data, mods)

# Create violin plots for each signature
vlns <- lapply(1:length(metsigs), function(i) {
  md <- sobjx@meta.data
  md <- md[, c('Genotype', names(metsigs)[i])]
  colnames(md)[2] <- 'sig'
  
  vln <- ggplot(md, aes(Genotype, sig, fill = Genotype)) +
    geom_violin(alpha = 1) +
    geom_jitter(width = 0.1, size = 0.1, alpha = 0.2) +
    theme_linedraw() +
    scale_fill_manual(values = c('steelblue', 'orange', 'firebrick')) +
    ylab('') + 
    xlab('') +
    labs(title = names(metsigs)[i]) +  # FIXED: Added + operator here (line 276 bug)
    theme(axis.title.y = element_text(size = 8)) +
    NoLegend()
  
  return(vln)
})

# Add axis label to first plot
vlns[[1]] <- vlns[[1]] + 
  ylab('Module Score') + 
  xlab('') +
  theme(axis.title.y = element_text(size = 8))

plot_2d_2e <- patchwork::wrap_plots(vlns, ncol = 2)

pdfname <- file.path(BASE_OUTPUT_DIR, '2B_Malignant_Met_genes_module.pdf')
pdf(pdfname, height = 3, width = 6)
print(plot_2d_2e)
dev.off()

message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 2 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
