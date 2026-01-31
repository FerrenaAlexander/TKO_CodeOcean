# =============================================================================
# Figure 8: Myogenesis Analysis
# =============================================================================
# Description: Generate Figure 8 panels showing myogenesis GSEA, leading edge
#              gene heatmap, transcription factor dotplots, and SCENIC regulon
#              activity for myogenic factors
# =============================================================================

# Load required libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(fgsea)
library(SCENIC)
library(SCopeLoomR)
library(edgeR)

# Set seed for reproducibility
set.seed(2021)

# Source configuration file
source("config.R")

# -----------------------------------------------------------------------------
# Check Required Input Files
# -----------------------------------------------------------------------------

required_files <- c(
  SEURAT_CELLSUBTYPES,
  MSIGDB_PATHWAYS,
  PATHWAY_OBJECTS,
  DE_TKO_DKO_MALIGNANT,
  DE_DKOAA_DKO_MALIGNANT,
  DE_TKO_DKOAA_MALIGNANT,
  SCENIC_LOOM,
  METADATA
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

# Pathways
pathways <- readRDS(MSIGDB_PATHWAYS)
pathways <- as.data.frame(pathways)

# Pathway objects
pwayresfile <- PATHWAY_OBJECTS
de_l <- readRDS(pwayresfile)

# DE filepaths
filepaths <- c(DE_TKO_DKO_MALIGNANT, DE_DKOAA_DKO_MALIGNANT, DE_TKO_DKOAA_MALIGNANT)
names(filepaths) <- c('TKO_vs_DKO', 'DKOAA_vs_DKO', 'TKO_vs_DKOAA')

# Get myogenesis leading edge genes
# pway_enrichment <- readRDS(PATHWAY_OBJECTS)
# if ("pathway" %in% colnames(pway_enrichment)) {
#   leading_edge <- pway_enrichment[pway_enrichment$pathway == "HALLMARK_MYOGENESIS", "leadingEdge"]
#   if (length(leading_edge) > 0) {
#     leading_edge <- str_split_fixed(leading_edge, '/', Inf)[1, ]
#   }
# } else {
#   # Try reading from the TKO vs DKO pathway table
#   pway_enrichment <- read.csv(file.path(
#     dirname(PATHWAY_OBJECTS), 
#     "TKO_vs_DKO/HALLMARK/Malignant/pathwaytable.csv"
#   ))
#   leading_edge <- pway_enrichment[pway_enrichment$pathway == "HALLMARK_MYOGENESIS", "leadingEdge"]
#   leading_edge <- str_split_fixed(leading_edge, '/', Inf)[1, ]
# }
leading_edge <- de_l$pathway_analysis_mainlist_comps$TKO_vs_DKO$HALLMARK$Malignant$gseares
leading_edge <- leading_edge[leading_edge$pathway == 'HALLMARK_MYOGENESIS', "leadingEdge"]
leading_edge <- leading_edge[[1]]

# Also get DKOAA leading edge (from TKO vs DKOAA comparison)
dkoaa_le <- de_l$pathway_analysis_mainlist_comps$TKO_vs_DKOAA$HALLMARK$Malignant$gseares
dkoaa_le <- dkoaa_le[dkoaa_le$pathway == "HALLMARK_MYOGENESIS", "leadingEdge"]
dkoaa_le <- dkoaa_le[[1]]

# Get intersect of leading edges
if (length(leading_edge) > 0 & length(dkoaa_le) > 0) {
  leading_edge <- intersect(leading_edge, dkoaa_le)
}

# SCENIC loom
loom <- open_loom(SCENIC_LOOM)

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
# Figure 8A-B: Myogenesis GSEA Plots
# -----------------------------------------------------------------------------

message("Generating Figure 8A-B: Myogenesis GSEA plots...")

pwayname <- "HALLMARK_MYOGENESIS"
pway_to_plot <- list(pathways[pathways$gs_name == pwayname, "gene_symbol"])
names(pway_to_plot) <- pwayname
pway_to_plot_l <- pway_to_plot
pway_to_plot <- pway_to_plot[[1]]

# Generate GSEA plots for each comparison
gseares <- lapply(1:length(filepaths), function(i) {
  fp <- filepaths[i]
  comp <- names(filepaths)[i]
  nicecomp <- gsub('_', ' ', comp)
  deres <- read.csv(fp)
  
  scores <- setNames(-log10(deres$PValue) * sign(deres$logFC), nm = deres$gene_symbol)
  scores <- sort(scores, decreasing = TRUE)
  
  d <- fgsea::plotEnrichment(pathway = pway_to_plot, stats = scores) +
    ggtitle(paste0(nicecomp)) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 10),
      axis.title.y = element_text(size = 7),
      legend.title = element_text(size = 7),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    ylab('Enrichment Score') + 
    xlab('')
  
  gr <- fgsea(pathways = pway_to_plot_l, stats = scores)
  
  # Get original NES/FDR if available
  if (!is.null(de_l$pathway_analysis_mainlist_comps[[comp]])) {
    gseares <- de_l$pathway_analysis_mainlist_comps[[comp]]$HALLMARK$Malignant$gseares
    if (pwayname %in% gseares$pathway) {
      gr <- gseares[gseares$pathway == pwayname, ]
    }
  }
  
  list(d = d, gr = gr)
})

names(gseares) <- names(filepaths)

pdfname <- file.path(BASE_OUTPUT_DIR, '8A_8B_myogenesis_malignant_GSEA.pdf')
pdf(pdfname, height = 5, width = 3)
print(wrap_plots(list(
  gseares$TKO_vs_DKO$d,
  gseares$DKOAA_vs_DKO$d,
  gseares$TKO_vs_DKOAA$d
), ncol = 1))
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 8C: Myogenesis Leading Edge Heatmap
# -----------------------------------------------------------------------------

message("Generating Figure 8C: Myogenesis leading edge heatmap...")

# Get pseudobulk matrix
mat <- SDAP::pseudobulk(sobjx, grouping_colname_in_md = "Code", 
                        assay = 'RISC', slot = 'counts')

# EdgeR normalize
group <- colnames(mat)
group <- str_split_fixed(group, '_', 2)[, 1]

eobj <- DGEList(counts = mat, group = group)
eobj <- calcNormFactors(eobj)
mat <- cpm(eobj)
mat <- log1p(mat)
rm(eobj, group)

# Get leading edge genes
matgenes <- rownames(mat)
le_genes <- matgenes[matgenes %in% leading_edge]
mat <- mat[match(le_genes, rownames(mat)), ]
mat <- as.matrix(mat)
mat <- t(scale(t(mat)))

# Labels
rowlabs <- data.frame(
  genes = rownames(mat),
  idx = 1:nrow(mat),
  newlabs = rownames(mat)
)

hm <- ComplexHeatmap::Heatmap(
  mat, 
  row_names_gp = grid::gpar(fontsize = 8),
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  column_split = c(rep(c("DKO", 'DKOAA', 'TKO'), each = 3), 'TKO'),
  row_gap = unit(1.5, "mm"),
  column_gap = unit(1.5, "mm"),
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 7),
  column_title = 'Hallmark Myogenesis Leading Edge Genes',
  name = 'Scaled\nLog1p\nEdgeR-norm\nPseudobulk\nCounts'
)

pdfname <- file.path(BASE_OUTPUT_DIR, '8C_Myogenesis_leadingEdge_Heatmap.pdf')
pdf(pdfname, height = 4.2, width = 5.2)
print(hm)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 8D: Myogenesis TF Dotplot
# -----------------------------------------------------------------------------

message("Generating Figure 8D: Myogenesis TF dotplot...")

genes <- c('Pax7', 'Myod1', 'Myf5', 'Myf6', 'Myog')
genes <- genes[genes %in% rownames(sobjx)]

dp <- DotPlot(sobjx, rev(genes), group.by = 'Genotype') + 
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  ylab('Genotype')

pdfname <- file.path(BASE_OUTPUT_DIR, '8D_Myogenesis_TF_dotplot.pdf')
pdf(pdfname, height = 4, width = 4)
print(dp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 8E-F: Myogenesis TF SCENIC Regulons and Gene Expression
# -----------------------------------------------------------------------------

message("Generating Figure 8E-F: SCENIC regulon and gene expression...")

# Extract SCENIC data
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
close_loom(loom)

# Get regulon AUC score matrix
regmatseurat <- regulonAUC@assays@data$AUC
regmatseurat <- regmatseurat[, match(colnames(sobjx), colnames(regmatseurat))]
regmatseurat <- Seurat::CreateAssayObject(regmatseurat)
sobjx[["SCENIC"]] <- regmatseurat

# Order codes by pathology
pbmd <- as.data.frame(readxl::read_excel(METADATA))
sobjx$Code_Ord <- factor(
  sobjx$Code, 
  levels = c('DKO_1', 'DKO_3', 'DKOAA_1', 'TKO_1', 'TKO_3',  # Osteo
             'DKO_2', 'DKOAA_3', 'TKO_2', 'TKO_4',
             'DKOAA_2')
)

genes <- c("Pax7(+)", 'Myod1(+)', 'Myf5(+)', 'Myf6(+)', 'Myog(+)')
xintlines <- c(5, 9) + 0.5

dp_regulons <- DotPlot(sobjx, features = rev(genes), assay = 'SCENIC', group.by = 'Code_Ord') +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0.1), "cm")
  ) +
  ylab('Genotype') +
  geom_hline(yintercept = xintlines, linewidth = 0.2)

actualgenes <- str_split_fixed(genes, pattern = '\\(', n = 2)[, 1]

dp_gene <- DotPlot(sobjx, features = rev(c(actualgenes)), assay = 'RISC', group.by = 'Code_Ord') +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0.1, 0, 0), "cm")
  ) +
  ylab('Genotype') +
  geom_hline(yintercept = xintlines, linewidth = 0.2)

dp_patch <- dp_gene + dp_regulons

pdfname <- file.path(BASE_OUTPUT_DIR, '8E_8F_Myogenesis_TFs_DOTPLOT.pdf')
pdf(pdfname, height = 4, width = 9)
print(dp_patch)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 8H: SKP2-related Genes Dotplot
# -----------------------------------------------------------------------------

message("Generating Figure 8H: SKP2-related genes dotplot...")

genes <- c('Skp2', 'Myod1', 'Cdkn1c', 'Fbxl13', 'Asb5')

dp_skp2genes <- DotPlot(sobjx, features = rev(c(genes)), assay = 'RISC', group.by = 'Code_Ord') +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0.1, 0, 0), "cm")
  ) +
  ylab('Genotype') +
  geom_hline(yintercept = xintlines, linewidth = 0.2)

pdfname <- file.path(BASE_OUTPUT_DIR, '8H_Myogenesis_skp2genes_DOTPLOT.pdf')
pdf(pdfname, height = 4, width = 4.5)
print(dp_skp2genes)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 8 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
