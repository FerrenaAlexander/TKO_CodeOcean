# =============================================================================
# Figure 4: E2F/Myc/SCENIC Analysis
# =============================================================================
# Description: Generate Figure 4 panels showing E2F regulon heatmap, apoptosis
#              gene dotplot, Myc GSEA, and Myc leading edge gene analysis
# =============================================================================

# Load required libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(edgeR)
library(SCENIC)
library(SCopeLoomR)

# Set seed for reproducibility
set.seed(2021)

# Source configuration file
source("config.R")

# -----------------------------------------------------------------------------
# Check Required Input Files
# -----------------------------------------------------------------------------

required_files <- c(
  SCENIC_LOOM,
  SEURAT_CELLSUBTYPES,
  PATHWAY_OBJECTS,
  MSIGDB_PATHWAYS,
  PATHWAY_TKO_DKO,
  PATHWAY_DKOAA_DKO,
  MYCHIGH_TKO_DKO_NCELLS,
  MYCHIGH_TKO_DKO_PROP,
  MYCHIGH_TKO_DKO_COMP,
  MYCHIGH_DKOAA_DKO_NCELLS,
  MYCHIGH_DKOAA_DKO_PROP,
  MYCHIGH_DKOAA_DKO_COMP,
  MYCHIGH_TKO_DKOAA_COMP
)

if (!check_input_files(required_files)) {
  stop("Missing required input files. Please check INPUT_FILES_LIST.md")
}

# -----------------------------------------------------------------------------
# Read Input Data
# -----------------------------------------------------------------------------

message("Loading input data...")

# SCENIC loom file
loom <- open_loom(SCENIC_LOOM)

# Malignant subsetted object
sobjx <- readRDS(SEURAT_CELLSUBTYPES)$Malignant

# DE pathway list
de_res_pway_list <- readRDS(PATHWAY_OBJECTS)
m_bycluster_crosscondition_de_comps <- de_res_pway_list$m_bycluster_crosscondition_de_comps

# Pathways
pathways <- readRDS(MSIGDB_PATHWAYS)
pathways <- as.data.frame(pathways)

# Myc leading edge genes
pway_enrichment <- read.csv(PATHWAY_TKO_DKO)
leading_edge <- pway_enrichment[pway_enrichment$pathway == "HALLMARK_MYC_TARGETS_V1", "leadingEdge"]
leading_edge <- str_split_fixed(leading_edge, '/', Inf)[1, ]

pway_enrichment_dkoaa <- read.csv(PATHWAY_DKOAA_DKO)
dkoaa_le <- pway_enrichment_dkoaa[pway_enrichment_dkoaa$pathway == "HALLMARK_MYC_TARGETS_V1", "leadingEdge"]
dkoaa_le <- str_split_fixed(dkoaa_le, '/', Inf)[1, ]

leading_edge <- intersect(leading_edge, dkoaa_le)

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
# Figure 4A: E2F Regulon SCENIC Heatmap
# -----------------------------------------------------------------------------

message("Generating Figure 4A: E2F regulon SCENIC heatmap...")

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

# Get E2F and Myc regulons
regmat <- regmatseurat@data
regscores <- rownames(regmat)[grepl('E2f', rownames(regmat))]
regscores2 <- rownames(regmat)[grepl('Myc', rownames(regmat))]
regscores <- c(regscores, regscores2)
genes <- str_sort(regscores, numeric = TRUE)

# Average expression heatmap
avgmat <- AverageExpression(sobjx, assays = 'SCENIC', features = regscores, 
                           layer = 'data', group.by = 'Code')[[1]]
avgmat2 <- as.matrix(avgmat)
avgmat2 <- t(scale(t(avgmat2)))

colsplit <- colnames(avgmat2)
colsplit <- str_split_fixed(colsplit, pattern = '-', 2)[, 1]

scenic_HM <- ComplexHeatmap::Heatmap(
  avgmat2, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  column_split = colsplit,
  row_gap = unit(1.5, "mm"),
  column_gap = unit(1.5, "mm"),
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 5),
  name = 'Scaled\nSCENIC\nAUCell\nRegulon\nScore',
  column_title = 'SCENIC E2F regulon scores',
  column_title_gp = grid::gpar(fontsize = 10)
)

pdfname <- file.path(BASE_OUTPUT_DIR, '4A_E2f_heatmap_SCENIC_e2fregulons_withmyc.pdf')
pdf(pdfname, height = 2.5, width = 3.5)
print(scenic_HM)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 4B: Myc/E2F/Apoptosis Dotplot
# -----------------------------------------------------------------------------

message("Generating Figure 4B: Apoptosis gene dotplot...")

genes <- c('Bbc3', 'Bid', 'Casp3', 'Dapk2', 'Pmaip1')

dp <- DotPlot(sobjx, features = rev(genes), group.by = 'Code') +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  xlab('') + 
  ylab('')

pdfname <- file.path(BASE_OUTPUT_DIR, '4B_Myc_E2f_apoptosis_dotplot.pdf')
pdf(pdfname, height = 4, width = 5)
print(dp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 4C and 4D: Myc GSEA
# -----------------------------------------------------------------------------

message("Generating Figure 4C and 4D: Myc GSEA plots...")

# Get TKO vs DKO DE results
deres <- m_bycluster_crosscondition_de_comps$TKO_vs_DKO$Malignant

# Get Myc pathway genes
mycgenes <- pathways[pathways$gs_name == 'HALLMARK_MYC_TARGETS_V1', 'gene_symbol']

# Create ranked gene list
scores <- setNames(-log10(deres$PValue) * sign(deres$logFC), nm = deres$gene_symbol)
scores <- sort(scores, decreasing = TRUE)

# GSEA enrichment plot
myc_pway_list <- list(HALLMARK_MYC_TARGETS_V1 = mycgenes)
gsea_plot <- fgsea::plotEnrichment(pathway = mycgenes, stats = scores) +
  ggtitle('Hallmark Myc Targets V1\nTKO vs DKO') +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 12),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 10)
  ) +
  ylab('Enrichment Score') + 
  xlab('')

pdfname <- file.path(BASE_OUTPUT_DIR, '4C_4D_malignant_Myc_gsea.pdf')
pdf(pdfname, height = 4, width = 5)
print(gsea_plot)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 4E: Myc Leading Edge Heatmap
# -----------------------------------------------------------------------------

message("Generating Figure 4E: Myc leading edge heatmap...")

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

# Create heatmap
hm_avg <- ComplexHeatmap::Heatmap(
  mat, 
  row_names_gp = grid::gpar(fontsize = 7),
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  column_split = c(rep("DKO", 3), rep('DKOAA', 3), rep('TKO', 4)),
  row_gap = unit(1.5, "mm"),
  column_gap = unit(1.5, "mm"),
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 7),
  column_title = 'Hallmark Myc Targets V1 Leading Edge Genes',
  name = 'Scaled\nPseudobulk\nCounts'
)

pdfname <- file.path(BASE_OUTPUT_DIR, '4E_Myc_leadingedge_heatmap.pdf')
pdf(pdfname, height = 4.2, width = 5.2)
print(hm_avg)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 2F: Myc-Hi Compositional Analysis
# -----------------------------------------------------------------------------

message("Generating Figure 2F: Myc-Hi compositional analysis...")

# Read compositional data
cellstab <- read.csv(MYCHIGH_TKO_DKO_NCELLS)
proptab <- read.csv(MYCHIGH_TKO_DKO_PROP)
pres <- read.csv(MYCHIGH_TKO_DKO_COMP)

d_cellstab <- read.csv(MYCHIGH_DKOAA_DKO_NCELLS)
d_proptab <- read.csv(MYCHIGH_DKOAA_DKO_PROP)
d_pres <- read.csv(MYCHIGH_DKOAA_DKO_COMP)

# Prepare proportion table
comp_proptab <- proptab
rownames(comp_proptab) <- pres$BaselineProp.clusters
comp_proptab <- comp_proptab[1, ]
pres <- pres[1, ]

# Make new proportion table with DKO, DKOAA, TKO
new_comp_proptab <- proptab[1, 2:4]
new_comp_proptab <- c(new_comp_proptab, d_proptab[1, 5:7])
new_comp_proptab <- unlist(new_comp_proptab)
new_comp_proptab <- unlist(c(new_comp_proptab, proptab[1, 5:8]))
comp_proptab <- new_comp_proptab

col_fun <- circlize::colorRamp2(c(0, max(comp_proptab)), c("white", "red"))
comp_proptab <- t(as.matrix(comp_proptab))

hmprop_comp <- Heatmap(
  comp_proptab, 
  name = "Proportion", 
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  border_gp = gpar(col = "black", lwd = 1),
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  column_split = c(1, 1, 1, 2, 2, 2, 4, 4, 4, 4),
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 6),
  column_title = 'Myc-Response Hi Proportions',
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", comp_proptab[i, j]), x, y, gp = gpar(fontsize = 8))
  }
)

pdfname <- file.path(BASE_OUTPUT_DIR, '2F_Myc_Hi_CompositionalAnalysis.pdf')
pdf(pdfname, height = 1.3, width = 3.5)
print(hmprop_comp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 2G and H: Myc Module Score Violin Plots
# -----------------------------------------------------------------------------

message("Generating Figure 2G and H: Myc module score violin plots...")

# Calculate Myc module scores
pathways_df <- readRDS(MSIGDB_PATHWAYS)
pathways_df <- as.data.frame(pathways_df)
pathways_df <- pathways_df[pathways_df$gs_cat == 'H', ]
mycgenes <- pathways_df[pathways_df$gs_name == 'HALLMARK_MYC_TARGETS_V1', 'gene_symbol']
mycgenes2 <- pathways_df[pathways_df$gs_name == 'HALLMARK_MYC_TARGETS_V2', 'gene_symbol']

mycgenes_list <- list(
  HALLMARK_MYC_TARGETS_V1 = mycgenes,
  HALLMARK_MYC_TARGETS_V2 = mycgenes2
)

sobjx <- AddModuleScore(sobjx, features = mycgenes_list, name = 'HALLMARK_MYC_TARGETS_V')

# Module score violin plot
sobjx$Genotype <- factor(sobjx$Genotype, levels = c('DKO', 'DKOAA', 'TKO'))
vp <- ggplot(sobjx@meta.data, aes(Genotype, HALLMARK_MYC_TARGETS_V1, fill = Genotype)) +
  geom_violin() +
  geom_jitter(size = 0.00001, alpha = 0.3) +
  scale_fill_manual(values = c('steelblue', 'orange', 'firebrick')) +
  ggtitle('Hallmark Myc Targets V1') +
  xlab('') + 
  ylab('Module Score') +
  theme_light() +
  NoLegend() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 15),
    legend.title = element_text(size = 10)
  )

# Gene expression violin plot
md <- sobjx@meta.data
gene <- 'Ccna2'
md[, 'Expression'] <- sobjx@assays$RISC@data[gene, ]

vp2 <- ggplot(md, aes(Genotype, Expression, fill = Genotype)) +
  geom_violin() +
  geom_jitter(size = 0.00001, alpha = 0.3) +
  scale_fill_manual(values = c('steelblue', 'orange', 'firebrick')) +
  ggtitle(paste0(gene, ' expression')) +
  xlab('Genotype') + 
  ylab('Expression') +
  theme_light() +
  NoLegend() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 15),
    legend.title = element_text(size = 10)
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '2G_H_MycModuleScore_ViolinPlot.pdf')
pdf(pdfname, height = 4, width = 4)
print(vp)
print(vp2)
print(vp / vp2)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 4 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
