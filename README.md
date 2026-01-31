# CodeOcean Capsule - Figure Generation Scripts

## Overview

This repository contains cleaned and standardized R scripts for generating all main scRNA-seq figures for the manuscript. All scripts have been organized with a centralized configuration system for easy deployment to CodeOcean.

## Quick Start

1. **Update Configuration**
   ```r
   # Edit config.R and change this line:
   BASE_INPUT_DIR <- "path/to/codeoceaninput"
   # to your actual CodeOcean data directory path
   ```

2. **Verify Input Files**
   ```r
   source("config.R")
   # Will automatically check if all paths are correctly configured
   ```

3. **Run Individual Figure**
   ```r
   source("codeocean_fig1_clean.R")
   ```

4. **Run All Figures**
   ```r
   source("run_all_figures.R")
   ```

## File Structure

```
.
├── README.md                     # This file
├── config.R                      # Configuration with all paths
├── INPUT_FILES_LIST.md          # Complete list of required input files (36 files)
├── CLEANUP_SUMMARY.md           # Details on changes, bugs fixed, improvements
├── run_all_figures.R            # Master script to run all figures
├── codeocean_fig1_clean.R       # Figure 1: Integrated UMAP Analysis
├── codeocean_fig2_clean.R       # Figure 2: GSEA and EMT Analysis
├── codeocean_fig3_clean.R       # Figure 3: T Cell Analysis
├── codeocean_fig4_clean.R       # Figure 4: E2F/Myc/SCENIC Analysis
├── codeocean_fig5_clean.R       # Figure 5: CNV Analysis
├── codeocean_fig6_clean.R       # Figure 6: Human OS Comparison
├── codeocean_fig7_clean.R       # Figure 7: Pathology Subtypes
└── codeocean_fig8_clean.R       # Figure 8: Myogenesis Analysis
```

## Required Inputs

All scripts require 36 input files organized in a specific directory structure. See `INPUT_FILES_LIST.md` for the complete list and organization.

### Key Input Files:
- **Main Seurat objects**: Integrated and cell subtype-specific
- **InferCNV results**: CNV matrices and PCA results
- **Pathway analysis**: GSEA results for multiple comparisons
- **Differential expression**: DE results for TKO/DKO/DKOAA comparisons
- **SCENIC results**: Regulon activity loom file
- **Metadata**: Sample annotation spreadsheet

## Required R Packages

```r
# CRAN packages
install.packages(c("tidyverse", "patchwork", "readxl", "ggVennDiagram", 
                   "circlize", "grid", "RColorBrewer", "plyr", "progress"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "ComplexHeatmap", "fgsea", "edgeR", 
                       "SCENIC", "SCopeLoomR"))
```

## Key Improvements

### 1. Centralized Configuration
- All file paths defined in `config.R`
- Single variable to update for CodeOcean deployment
- Automatic input file validation

### 2. Bug Fixes
- **Figure 2**: Fixed missing `+` operator (line 276)
- **All scripts**: Removed hardcoded absolute paths
- **Consistent output structure**: All outputs to `results/figures/`

### 3. Code Quality
- Consistent style and formatting
- Extensive comments and section headers
- Progress messages for debugging
- Error handling and validation

### 4. Documentation
- Each script has clear header documentation
- Function definitions with descriptions
- Input/output specifications

## Output Structure

All figures are saved to `results/figures/` with descriptive filenames:

```
results/figures/
├── 1B_umap_celltype.pdf
├── 1C_dotplot_markers.pdf
├── 1D_umap_genotype.pdf
├── 1E_violplot_extremeCNVs.pdf
├── 1F_CNV_pca_allcells_colby_celltype.pdf
├── 2A_2B_DE_malignant_GSEA.pdf
├── ... (and many more)
└── execution_summary.csv  (if using run_all_figures.R)
```

## Execution Time

Approximate runtime per script (depends on hardware):
- Figure 1: ~2-5 minutes
- Figure 2: ~3-7 minutes
- Figure 3: ~5-15 minutes (most complex)
- Figure 4: ~5-10 minutes (SCENIC processing)
- Figure 5: ~10-20 minutes (CNV calculations)
- Figure 6: ~3-7 minutes
- Figure 7: ~5-10 minutes
- Figure 8: ~5-10 minutes

**Total: ~40-90 minutes** for all figures

## System Requirements

- **R version**: >= 4.0
- **RAM**: >= 32GB recommended
- **Storage**: ~50GB for input data + outputs
- **CPU**: Multi-core beneficial but not required

## Troubleshooting

### Common Issues:

1. **Missing input files**
   ```r
   source("config.R")
   # Check error messages for which files are missing
   ```

2. **Package installation errors**
   - SCENIC requires HDF5 libraries: `sudo apt-get install libhdf5-dev`
   - Check Bioconductor version compatibility

3. **Memory errors**
   - Close other applications
   - Consider running figures individually rather than all at once
   - Increase swap space if on Linux

4. **Path issues**
   - Ensure `BASE_INPUT_DIR` in `config.R` is correct
   - Use forward slashes (/) even on Windows
   - Avoid spaces in directory names

### Getting Help:

1. Review `CLEANUP_SUMMARY.md` for detailed information
2. Check `INPUT_FILES_LIST.md` for file requirements
3. Examine script-specific comments for details
4. Verify all required packages are installed

## Script Descriptions

### Figure 1: Integrated UMAP Analysis
- **Panels**: 1B-1F
- **Content**: Cell type UMAP, marker dotplot, genotype UMAP, CNV violin, CNV PCA
- **Key inputs**: Integrated Seurat object, CNV PCA results

### Figure 2: GSEA and EMT Analysis
- **Panels**: 2A-2E
- **Content**: GSEA dotplots (TKO/DKOAA vs DKO), EMT gene dotplot, metastasis signatures
- **Key inputs**: GSEA results, malignant subset object

### Figure 3: T Cell Analysis
- **Panels**: 3A-3F
- **Content**: GSEA across celltypes, T cell markers, TIL signatures, compositional analysis
- **Key inputs**: T cell object, annotations, compositional data

### Figure 4: E2F/Myc/SCENIC Analysis
- **Panels**: 4A-4D + supplemental
- **Content**: E2F regulon heatmap, apoptosis dotplot, Myc GSEA, Myc leading edge
- **Key inputs**: SCENIC loom, pathway data, malignant object

### Figure 5: CNV Analysis
- **Panels**: 4B-4G (CNV-focused)
- **Content**: CNV burden violins, CNV PCA, CNV-DEG Venn diagrams, differential CNV heatmap
- **Key inputs**: CNV matrices, inferCNV results, DE results

### Figure 6: Human OS Comparison
- **Panels**: 6A-6C
- **Content**: Human OS UMAP, marker dotplot, module score heatmap
- **Key inputs**: Human OS Seurat object with module scores

### Figure 7: Pathology Subtypes
- **Panels**: 7A-7I
- **Content**: Label transfer UMAPs, module violins, marker dotplots, compositional barplots
- **Key inputs**: Pathology-annotated object, atlas markers

### Figure 8: Myogenesis Analysis
- **Panels**: 8A-8H
- **Content**: Myogenesis GSEA, leading edge heatmap, TF dotplots, SCENIC regulons
- **Key inputs**: Malignant object, SCENIC loom, pathway data

## Citation

If you use these scripts, please cite the original manuscript (details to be added).

## License

[Add license information]

## Contact

For questions or issues:
- Review documentation files first
- Check `CLEANUP_SUMMARY.md` for known issues
- Contact: [Add contact information]

## Version History

- **v1.0** (2026-01-27): Initial cleaned version
  - Centralized configuration system
  - Fixed critical bugs
  - Standardized all scripts
  - Comprehensive documentation

---

**Note**: These scripts are optimized for CodeOcean deployment but can be run in any R environment with the required packages and data.
