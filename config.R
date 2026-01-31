# =============================================================================
# Configuration File for CodeOcean Capsule - SIMPLIFIED VERSION
# =============================================================================
# This file centralizes all input/output paths for the analysis scripts
# ALL INPUT FILES GO IN A SINGLE "DATA" FOLDER
# Created: 2026-01-27
# Updated: 2026-01-28 - Simplified to single DATA folder
# =============================================================================

# -----------------------------------------------------------------------------
# Directory Paths
# -----------------------------------------------------------------------------

# Base input directory - UPDATE THIS PATH
# For local use: set to your DATA folder location
# For CodeOcean: set to "../data" 
BASE_INPUT_DIR <- "DATA"  # <-- CHANGE THIS TO YOUR DATA FOLDER PATH

# Base output directory
BASE_OUTPUT_DIR <- "results/figures"

# Create output directory if it doesn't exist
if (!dir.exists(BASE_OUTPUT_DIR)) {
  dir.create(BASE_OUTPUT_DIR, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Main Seurat Objects (2 files)
# -----------------------------------------------------------------------------

# Integrated Seurat object with CNVs
SEURAT_INTEGRATED <- file.path(BASE_INPUT_DIR, "sobjint_WITH_CNVs.rds")

# Cell subtype objects list
SEURAT_CELLSUBTYPES <- file.path(BASE_INPUT_DIR, "ctsobjlist.rds")

# -----------------------------------------------------------------------------
# InferCNV Data (5 files)
# -----------------------------------------------------------------------------

# PCA of CNVs - all cells
CNV_PCA_ALL <- file.path(BASE_INPUT_DIR, "harmonized_pca_allcells.rds")

# PCA of CNVs - malignant only
CNV_PCA_MALIGNANT <- file.path(BASE_INPUT_DIR, "harmonized_pca_maligonly.rds")

# CNV matrix
CNV_MATRIX <- file.path(BASE_INPUT_DIR, "harmonized_CNV_matrix.rds")

# Cross-condition CNV comparisons
CNV_CROSSCOND_MEAN <- file.path(BASE_INPUT_DIR, "CROSSCOND_mean.rds")
CNV_CROSSCOND_AMPS <- file.path(BASE_INPUT_DIR, "CROSSCOND_counts_amps.rds")

# -----------------------------------------------------------------------------
# Pathway Analysis Results (5 files)
# -----------------------------------------------------------------------------

# GSEA pathway tables
PATHWAY_TKO_DKO <- file.path(BASE_INPUT_DIR, "pathwaytable_TKODKO_HALLMARK.csv")
PATHWAY_DKOAA_DKO <- file.path(BASE_INPUT_DIR, "pathwaytable_DKOAADKO_HALLMARK.csv")
PATHWAY_TKO_DKOAA <- file.path(BASE_INPUT_DIR, "pathwaytable_TKODKOAA_HALLMARK.csv")

# Pathway analysis objects
PATHWAY_OBJECTS <- file.path(BASE_INPUT_DIR, "DE_pathways_plot_objects_list.rds")

# MSigDB pathways
MSIGDB_PATHWAYS <- file.path(BASE_INPUT_DIR, "msigdb_pathways.rds")

# -----------------------------------------------------------------------------
# Differential Expression Results (3 files)
# -----------------------------------------------------------------------------

DE_TKO_DKO_MALIGNANT <- file.path(BASE_INPUT_DIR, "DE_TKO_vs_DKO_Malignant.csv")
DE_DKOAA_DKO_MALIGNANT <- file.path(BASE_INPUT_DIR, "DE_DKOAA_vs_DKO_Malignant.csv")
DE_TKO_DKOAA_MALIGNANT <- file.path(BASE_INPUT_DIR, "DE_TKO_vs_DKOAA_Malignant.csv")

# -----------------------------------------------------------------------------
# T Cell Analysis (6 files)
# -----------------------------------------------------------------------------

TCELL_ANNOTATION <- file.path(BASE_INPUT_DIR, "Tcell_annotation.xlsx")
CELLTYPE_MARKERS <- file.path(BASE_INPUT_DIR, "ctmarkerlist.rds")
TIL_MARKERS <- file.path(BASE_INPUT_DIR, "41467_2021_23324_MOESM4_ESM.xlsx")

# Compositional analysis
TCELL_COMP_NCELLS <- file.path(BASE_INPUT_DIR, "Tcell_NumberCells.csv")
TCELL_COMP_PROP <- file.path(BASE_INPUT_DIR, "Tcell_ProportionCells.csv")
TCELL_COMP_TKO_DKO <- file.path(BASE_INPUT_DIR, "Tcell_Comp_TKO_vs_DKO.csv")
TCELL_COMP_DKOAA_DKO <- file.path(BASE_INPUT_DIR, "Tcell_Comp_DKOAA_vs_DKO.csv")

# -----------------------------------------------------------------------------
# SCENIC Results (1 file)
# -----------------------------------------------------------------------------

SCENIC_LOOM <- file.path(BASE_INPUT_DIR, "outscenic.loom")

# -----------------------------------------------------------------------------
# Myc-High Cell Analysis (7 files)
# -----------------------------------------------------------------------------

# TKO vs DKO
MYCHIGH_TKO_DKO_NCELLS <- file.path(BASE_INPUT_DIR, "MycHi_TKO_DKO_NumberCells.csv")
MYCHIGH_TKO_DKO_PROP <- file.path(BASE_INPUT_DIR, "MycHi_TKO_DKO_ProportionCells.csv")
MYCHIGH_TKO_DKO_COMP <- file.path(BASE_INPUT_DIR, "MycHi_TKO_DKO_Comp.csv")

# DKOAA vs DKO
MYCHIGH_DKOAA_DKO_NCELLS <- file.path(BASE_INPUT_DIR, "MycHi_DKOAA_DKO_NumberCells.csv")
MYCHIGH_DKOAA_DKO_PROP <- file.path(BASE_INPUT_DIR, "MycHi_DKOAA_DKO_ProportionCells.csv")
MYCHIGH_DKOAA_DKO_COMP <- file.path(BASE_INPUT_DIR, "MycHi_DKOAA_DKO_Comp.csv")

# TKO vs DKOAA
MYCHIGH_TKO_DKOAA_COMP <- file.path(BASE_INPUT_DIR, "MycHi_TKO_DKOAA_Comp.csv")

# -----------------------------------------------------------------------------
# Human Osteosarcoma Data (2 files)
# -----------------------------------------------------------------------------

HUMAN_OS_SEURAT <- file.path(BASE_INPUT_DIR, "human_OS_sobjint_withMurineModuleScores.rds")
HUMAN_OS_MARKERS <- file.path(BASE_INPUT_DIR, "human_OS_Markers_Celltypes.csv")

# -----------------------------------------------------------------------------
# Pathology Subtype Analysis (3 files)
# -----------------------------------------------------------------------------

PATHOLOGY_SEURAT <- file.path(BASE_INPUT_DIR, "pathology_sobj_malignantonly_WithLabelTransResOnly.rds")
PATHOLOGY_ATLAS_MARKERS <- file.path(BASE_INPUT_DIR, "pathology_atlasmarkers.csv")
PATHOLOGY_MALIGNANT_MARKERS <- file.path(BASE_INPUT_DIR, "pathology_malignant_pathmarkers.csv")

# -----------------------------------------------------------------------------
# Metadata (1 file)
# -----------------------------------------------------------------------------

METADATA <- file.path(BASE_INPUT_DIR, "metadata.xlsx")

# -----------------------------------------------------------------------------
# Helper Function to Check File Existence
# -----------------------------------------------------------------------------

check_input_files <- function(file_paths) {
  missing_files <- c()
  for (file_path in file_paths) {
    if (!file.exists(file_path)) {
      missing_files <- c(missing_files, file_path)
    }
  }
  
  if (length(missing_files) > 0) {
    message("\n=================================================================")
    message("WARNING: The following input files are missing:")
    message("=================================================================")
    for (f in missing_files) {
      message("  ✗ ", basename(f))
    }
    message("=================================================================")
    message("Expected location: ", BASE_INPUT_DIR)
    message("Please check that all files are copied to your DATA folder.")
    message("See INPUT_FILES_LIST_SIMPLIFIED.md for the complete list.")
    message("=================================================================\n")
    return(FALSE)
  } else {
    message("\n=================================================================")
    message("✓ All required input files found!")
    message("=================================================================\n")
    return(TRUE)
  }
}

# -----------------------------------------------------------------------------
# Print Configuration Summary
# -----------------------------------------------------------------------------

print_config_summary <- function() {
  cat("\n")
  cat("=================================================================\n")
  cat("Configuration Summary\n")
  cat("=================================================================\n")
  cat("Base Input Directory:  ", BASE_INPUT_DIR, "\n")
  cat("Base Output Directory: ", BASE_OUTPUT_DIR, "\n")
  cat("=================================================================\n")
  cat("\n")
}

# Print summary when config is loaded
print_config_summary()

# -----------------------------------------------------------------------------
# Quick Start Instructions
# -----------------------------------------------------------------------------

message("=================================================================")
message("QUICK START:")
message("=================================================================")
message("1. Copy all 36 input files to your DATA folder")
message("2. Update BASE_INPUT_DIR in this file if needed:")
message("   - For local use: BASE_INPUT_DIR <- 'path/to/your/DATA'")
message("   - For CodeOcean: BASE_INPUT_DIR <- '../data'")
message("3. Run a figure script: source('codeocean_fig1_clean.R')")
message("=================================================================\n")
