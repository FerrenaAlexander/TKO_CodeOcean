# =============================================================================
# Figure 3: T Cell Analysis
# =============================================================================
# Description: Generate Figure 3 panels showing GSEA across celltypes, T cell
#              UMAP and markers, TIL signatures, and compositional analysis
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
  PATHWAY_OBJECTS,
  MSIGDB_PATHWAYS,
  SEURAT_CELLSUBTYPES,
  TCELL_ANNOTATION,
  CELLTYPE_MARKERS,
  TIL_MARKERS,
  TCELL_COMP_NCELLS,
  TCELL_COMP_PROP,
  TCELL_COMP_TKO_DKO,
  TCELL_COMP_DKOAA_DKO,
  METADATA
)

if (!check_input_files(required_files)) {
  stop("Missing required input files. Please check INPUT_FILES_LIST.md")
}

# -----------------------------------------------------------------------------
# Read Input Data
# -----------------------------------------------------------------------------

message("Loading input data...")

# Read in TKO v DKO malignant DEGs
de_res_pway_list <- readRDS(PATHWAY_OBJECTS)

# MSigDB pathways
pathways <- readRDS(MSIGDB_PATHWAYS)

# T cell object
sobjx <- readRDS(SEURAT_CELLSUBTYPES)$Tcell

# T cell annotation
annot <- as.data.frame(readxl::read_excel(TCELL_ANNOTATION))

# Markers of subclusterings
mlist <- readRDS(CELLTYPE_MARKERS)

# Tumor infiltrating lymphocyte markers
tilmarkers <- as.data.frame(readxl::read_excel(TIL_MARKERS, sheet = 'Mouse_signatures', skip = 2))

# Compositional analysis data
cellstab <- read.csv(TCELL_COMP_NCELLS)
proptab <- read.csv(TCELL_COMP_PROP)
pres <- read.csv(TCELL_COMP_TKO_DKO)
pres_dkoaa <- read.csv(TCELL_COMP_DKOAA_DKO)

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

fix_underflow <- function(scores, logFC_vec) {
  if (missing(logFC_vec)) {
    logFC_vec <- NULL
    warning('Optimal to provide logFC values to logFC_vec for sorting')
  }
  
  bool <- (scores == Inf | scores == -Inf)
  tbl <- table(factor(bool, levels = c(FALSE, TRUE)))
  message(tbl['TRUE'], ' underflow genes detected')
  
  # Underflow, positive
  if (any(scores == Inf)) {
    scores_uf <- scores[scores == Inf]
    
    # Sort INF values by lfc
    if (!is.null(logFC_vec)) {
      logFC_vec_uf <- logFC_vec[names(logFC_vec) %in% names(scores_uf)]
      logFC_vec_uf <- sort(abs(logFC_vec_uf), decreasing = TRUE)
      scores_uf <- scores_uf[match(names(logFC_vec_uf), names(scores_uf))]
    }
    
    last <- scores[length(scores_uf) + 1]
    replacement <- c()
    for (i in 1:length(scores_uf)) {
      to_replace <- ifelse(i == 1, yes = last, no = replacement[i - 1])
      to_replace <- to_replace + 1  # For pos + 1
      replacement <- c(replacement, to_replace)
    }
    
    # Reverse, lo to hi
    replacement <- rev(replacement)
    names(replacement) <- names(scores_uf)
    scores[names(replacement)] <- replacement
  }
  
  # Underflow, negative
  if (any(scores == -Inf)) {
    scores <- rev(scores)
    scores_uf <- scores[scores == -Inf]
    
    # Sort INF values by lfc
    if (!is.null(logFC_vec)) {
      logFC_vec_uf <- logFC_vec[names(logFC_vec) %in% names(scores_uf)]
      logFC_vec_uf <- sort(abs(logFC_vec_uf), decreasing = TRUE)
      scores_uf <- scores_uf[names(logFC_vec_uf)]
    }
    
    last <- scores[length(scores_uf) + 1]
    replacement <- c()
    for (i in 1:length(scores_uf)) {
      to_replace <- ifelse(i == 1, yes = last, no = replacement[i - 1])
      to_replace <- to_replace - 1  # For negative only, minus 1
      replacement <- c(replacement, to_replace)
    }
    
    # For negative only
    replacement <- rev(replacement)
    names(replacement) <- names(scores_uf)
    scores[names(replacement)] <- replacement
    scores <- rev(scores)
  }
  
  scores <- sort(scores, decreasing = TRUE)
  return(scores)
}

# -----------------------------------------------------------------------------
# Figure 3A: GSEA Dotplot of TKO Up Signatures in All Celltypes
# -----------------------------------------------------------------------------

message("Generating Figure 3A: GSEA across celltypes...")

pathway_analysis_mainlist_comps <- de_res_pway_list$pathway_analysis_mainlist_comps
comps <- de_res_pway_list$comps

cp.font.size <- 8
compidx <- 1

# Get pathway analysis
pathway_analysis_mainlist <- pathway_analysis_mainlist_comps[[compidx]]

# Get comparison condition levels
c1 <- comps[compidx, 1]
c2 <- comps[compidx, 2]
lab <- comps[compidx, 3]

# Extract tables from all categories
cat_cpres_list <- lapply(pathway_analysis_mainlist, function(pwaycatlist) {
  clust_cpres <- lapply(1:length(pwaycatlist), function(clustidx) {
    clustname <- names(pwaycatlist)[clustidx]
    pwayres_DE_across_conditions_per_cluster <- pwaycatlist[[clustidx]]
    gseares_plot <- pwayres_DE_across_conditions_per_cluster$dp$data
    gseares_plot$cluster <- clustname
    gseares_plot$condition <- c1
    gseares_plot[sign(gseares_plot$NES) == -1, "condition"] <- c2
    return(gseares_plot)
  })
  dplyr::bind_rows(clust_cpres)
})

catdex <- 1
cpres_cat <- cat_cpres_list[[catdex]]
catname <- names(cat_cpres_list)[catdex]
cond <- c1

# Get result table for this condition
cpres_cat_cond <- cpres_cat[cpres_cat$condition == cond, , drop = FALSE]

# If no conditions, skip
if (nrow(cpres_cat_cond) > 0) {
  # Subselect categories if more than 30 total
  if (nrow(cpres_cat_cond) > 30) {
    cpres_cat_cond_sub <- cpres_cat_cond %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = -log10(padj)) %>%
      top_n(n = 5, wt = abs(NES)) %>%
      as.data.frame()
    
    cpres_cat_cond <- cpres_cat_cond[cpres_cat_cond$pathway %in% cpres_cat_cond_sub$pathway, ]
  }
  
  # Make sure orders are proper
  cpres_cat_cond$cluster <- factor(cpres_cat_cond$cluster, levels = unique(cpres_cat_cond$cluster))
  cpres_cat_cond$pathway <- factor(cpres_cat_cond$pathway, levels = rev(unique(cpres_cat_cond$pathway)))
  cpres_cat_cond$pathway2 <- plyr::mapvalues(
    cpres_cat_cond$pathway, 
    from = levels(cpres_cat_cond$pathway), 
    to = gsub(x = levels(cpres_cat_cond$pathway), pattern = 'HALLMARK ', '')
  )
  
  plot_3a <- ggplot(cpres_cat_cond, aes(x = cluster, y = pathway2, size = -log10(pval), col = NES)) +
    geom_point() +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = cp.font.size),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    scale_color_gradient2(
      low = 'steelblue', 
      high = 'red', 
      mid = 'white', 
      midpoint = 0, 
      name = 'Normalized\nEnrichment\nScore'
    ) +
    scale_size(range = c(2, 6), name = '-log10(pval)') +
    xlab('') + 
    ylab('') +
    ggtitle('TKO upregulated Hallmarks\ngenesets across celltypes')
  
  pdfname <- file.path(BASE_OUTPUT_DIR, '3A_GSEA_Celltypes.pdf')
  pdf(pdfname, height = 4, width = 6)
  print(plot_3a)
  dev.off()
  
  message("  Saved: ", pdfname)
}

# -----------------------------------------------------------------------------
# Figure 3B-3D: T Cell UMAP, Markers, and Dotplots
# -----------------------------------------------------------------------------

message("Generating Figure 3B-3D: T cell characterization...")

# Add annotations
annot$clust_assign <- paste0(annot$Cluster, ' - ', annot$Assignment)
sobjx$CellSubType <- sobjx$seurat_clusters
sobjx$CellSubType <- plyr::mapvalues(sobjx$CellSubType, from = annot$Cluster, to = annot$clust_assign)

# UMAP
plot_3b <- DimPlot(sobjx, group.by = 'CellSubType', label = TRUE, repel = TRUE) +
  theme_dimplot() +
  ggtitle('T Cell Subtypes')

pdfname <- file.path(BASE_OUTPUT_DIR, '3B_Tcell_UMAP.pdf')
pdf(pdfname, height = 4, width = 5)
print(plot_3b)
dev.off()
message("  Saved: ", pdfname)

# Markers
tm <- mlist$Tcell
tm$score <- (tm$pct.1 - tm$pct.2) * tm$avg_log2FC
tm <- tm[tm$avg_log2FC > 0.5, ]
tm <- tm[tm$pct.1 > 0.5, ]
tm <- tm[tm$pct.2 < 0.5, ]
tm <- tm[tm$p_val_adj < 0.05, ]
tm <- tm[order(tm$score, decreasing = TRUE), ]
tm <- tm[!duplicated(tm$gene), ]
tm <- tm[order(tm$cluster), ]

n <- 3
top <- tm %>% group_by(cluster) %>% top_n(n = n, wt = score) %>% as.data.frame()

plot_3c <- DotPlot(sobjx, features = rev(unique(top$gene)), group.by = 'CellSubType') +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    axis.title = element_blank(),
    legend.title = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

pdfname <- file.path(BASE_OUTPUT_DIR, '3C_Tcell_markers.pdf')
pdf(pdfname, height = 6, width = 5)
print(plot_3c)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 3E: TIL Module Scores
# -----------------------------------------------------------------------------

message("Generating Figure 3E: TIL signatures...")

tilmarkers <- as.list(tilmarkers)
tilmarkers <- lapply(tilmarkers, function(x) { x[!is.na(x)] })

sobjx_withmod <- AddModuleScore(sobjx, tilmarkers, name = names(tilmarkers))

rel_modcols <- colnames(sobjx_withmod@meta.data)[(length(colnames(sobjx@meta.data)) + 1):length(colnames(sobjx_withmod@meta.data))]
moddf <- sobjx_withmod@meta.data[, rel_modcols]
colnames(moddf) <- names(tilmarkers)
sobjx@meta.data <- cbind(sobjx@meta.data, moddf)

vln1 <- VlnPlot(sobjx, c('CD8_Tex'), ncol = 1) + 
  NoLegend() + 
  ggtitle('Andreatta et al CD8 Exhausted Signature') + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 9),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  xlab('') + 
  ylab('Module Score')

vln2 <- VlnPlot(sobjx, c('CD8_Tpex'), ncol = 1) + 
  NoLegend() +
  ggtitle('Precursor-Exhausted Signature') + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 9),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  xlab('') + 
  ylab('Module Score')

vln3 <- VlnPlot(sobjx, c('CD8_EffectorMemory'), ncol = 1) + 
  NoLegend() +
  ggtitle('CD8 Effector-Memory Signature') + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 9),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  xlab('') + 
  ylab('Module Score')

vln4 <- VlnPlot(sobjx, c('Treg'), ncol = 1) + 
  NoLegend() +
  ggtitle('Regulatory T Cell Signature') + 
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1, face = 'bold', size = 9),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  xlab('') + 
  ylab('Module Score')

pdfname <- file.path(BASE_OUTPUT_DIR, '3E_Tcell_TIL_modules.pdf')
pdf(pdfname, height = 5, width = 5)
(vln1 | vln2) / (vln3 | vln4)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Figure 3F: Compositional Analysis of TKO vs DKO T cells
# -----------------------------------------------------------------------------

message("Generating Figure 3F: T cell compositional analysis...")

# Prepare data
rownames(proptab) <- proptab$X
proptab <- proptab[, -1]

# Read metadata for ordering
pbmd <- as.data.frame(readxl::read_excel(METADATA))
pbmd$Condition <- factor(pbmd$Condition, c('DKO', 'DKOAA', 'TKO'))
pbmd <- pbmd[order(pbmd$Condition, decreasing = TRUE), ]

comp_proptab <- proptab[, pbmd$Code]
condition_vector_ordering <- factor(pbmd$Condition, levels = rev(levels(pbmd$Condition)))

# Prepare row annotations
jointpres <- pres[, c('X', 'PropMean.TKO', 'PropMean.DKO', 'PropRatio', 'P.Value', 'FDR')]
jointpres2 <- pres_dkoaa[, c('X', 'PropMean.DKOAA', 'PropMean.DKO', 'PropRatio', 'P.Value', 'FDR')]

colnames(jointpres) <- paste0('TKO_', colnames(jointpres))
colnames(jointpres2) <- paste0('DKOAA_', colnames(jointpres2))

jointpres <- jointpres[match(rownames(comp_proptab), jointpres$TKO_X), ]
jointpres2 <- jointpres2[match(rownames(comp_proptab), jointpres2$DKOAA_X), ]

jointpres <- cbind(jointpres, jointpres2)
rownames(jointpres) <- jointpres$TKO_X

# Meta-analysis
jointpres$avgPropTab_TKODKOAA_vs_DKO <- rowMeans(jointpres[, c('TKO_PropRatio', 'DKOAA_PropRatio')])
jointpres$avgPval <- rowMeans(jointpres[, c('TKO_P.Value', 'DKOAA_P.Value')])
jointpres$JointX <- jointpres$TKO_X
jointpres[jointpres$avgPval < 0.05, 'JointX'] <- paste0('* ', jointpres[jointpres$avgPval < 0.05, 'JointX'], ' *')

jointpres$score <- log(jointpres$avgPropTab_TKODKOAA_vs_DKO) * -log(jointpres$avgPval)
jointpres <- jointpres[order(jointpres$score, decreasing = TRUE), ]

propmat_annot <- jointpres[, c('TKO_PropMean.TKO', 'DKOAA_PropMean.DKOAA', 'TKO_PropMean.DKO')]

row_ha <- rowAnnotation(
  "Proportions" = anno_barplot(
    propmat_annot,
    gp = gpar(fill = c('firebrick', 'orange', 'steelblue'), 
              col = c('firebrick', 'orange', 'steelblue'))
  )
)

comp_proptab <- comp_proptab[match(rownames(jointpres), rownames(comp_proptab)), ]
rownames(comp_proptab) <- jointpres$JointX

col_fun <- circlize::colorRamp2(c(0, max(proptab)), c("white", "red"))

hmprop_comp <- Heatmap(
  comp_proptab, 
  name = "Proportion", 
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  border_gp = gpar(col = "black", lwd = 1),
  column_split = condition_vector_ordering,
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 7),
  right_annotation = row_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", comp_proptab[i, j]), x, y, gp = gpar(fontsize = 6))
  }
)

pdfname <- file.path(BASE_OUTPUT_DIR, '3F_Tcell_comp.pdf')
pdf(pdfname, height = 4, width = 5)
print(hmprop_comp)
dev.off()
message("  Saved: ", pdfname)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message("\n========================================")
message("Figure 3 generation complete!")
message("All outputs saved to: ", BASE_OUTPUT_DIR)
message("========================================\n")
