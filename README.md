# Brain Meta-Analysis Pipeline

A reproducible, scalable pipeline for performing both disease-specific and cross-phenotype meta-analyses of brain-related traits using bulk RNA-seq, GWAS summary statistics, and single-cell transcriptomic data.

---

## Table of Contents
1. [Motivation](#motivation)
2. [Features](#features)
3. [Repository Structure](#repository-structure)
4. [Installation](#installation)
5. [Data Preparation](#data-preparation)
6. [Workflow Overview](#workflow-overview)
   - [1. Disease-Specific Analysis](#1-disease-specific-analysis)
   - [2. Mega Analysis Across Phenotypes](#2-mega-analysis-across-phenotypes)
7. [Outputs & Interpretation](#outputs--interpretation)
8. [Customization & Extensions](#customization--extensions)
9. [Best Practices](#best-practices)
10. [Contributing](#contributing)
11. [License](#license)

---

## Motivation
Brain disorders and traits share complex molecular underpinnings. While individual studies yield important insights, integrating multiple high-throughput datasets (RNA-seq, GWAS, and single-cell data) across several phenotypes uncovers shared pathways, cell-type–specific drivers, and potential therapeutic targets.

This pipeline:
- Automates uniform preprocessing and statistical modeling for each phenotype.  
- Aggregates results for comparative meta-analysis to highlight convergent and divergent biology.  
- Supports reproducible research by adopting standardized directory layouts and container-ready scripts.

---

## Features
- **Phenotype-Specific Workflows**: Separate modules for Alzheimer’s, Parkinson’s, Schizophrenia, and other brain-related conditions.
- **Cross-Phenotype Integration**: Combine differential expression outputs, GWAS hits, and cell-type enrichment to map shared mechanisms.
- **Modular Code**: R scripts organized by function (preprocessing, DE analysis, integration) for easy customization.
- **Data-Driven Defaults**: Assumes common file formats but allows user-defined parameter overrides.
- **Visualization Templates**: Volcano plots, heatmaps, network diagrams, and dot plots for quick publication‑ready figures.

---

## Repository Structure
```
Brain_meta_analysis/
├── data/
│   ├── bulk_rnaseq/        # Raw counts + metadata subfolders per phenotype
│   ├── gwas/               # GWAS summary stats per phenotype
│   └── sc_rnaseq/          # Single-cell expression & metadata
│
├── disease_specific/       # Per-phenotype analysis pipelines
│   ├── Alzheimer/
│   │   ├── config.yml      # Paths & parameters
│   │   ├── preprocessing.R
│   │   ├── differential_expression.R
│   │   └── results/        # DE tables, plots
│   ├── Parkinson/          # ... similar structure
│   └── ...                 # Additional phenotypes
│
├── mega_analysis/
│   ├── codes/
│   │   ├── integrate_rnaseq_gwas.R    # Correlate gene expression with GWAS
│   │   ├── sc_integration.R           # Map bulk results onto cell types
│   │   └── summarize_results.R        # Aggregate cross-phenotype summaries
│   └── results/                       # Combined tables & figures
│
├── scripts/              # Utility scripts (e.g., download data, format conversion)
├── notebooks/            # Jupyter notebooks for exploratory analyses
├── renv.lock             # R dependency snapshot (optional)
└── README.md             # This documentation
```

---

## Installation
Ensure you have:
- **R** (>= 4.1.0) and **Command-line Git**
- **R packages**: `tidyverse`, `data.table`, `limma`, `edgeR`, `DESeq2`, `SingleCellExperiment`, `Matrix`, `readr`, `ggplot2`, etc.
- (Optional) **renv** for isolated environments

```bash
git clone https://github.com/isadeghi87/Brain_meta_analysis.git
cd Brain_meta_analysis
# (Optional) Initialize project library
Rscript -e "install.packages('renv'); renv::init(); renv::restore()"
```

---

## Data Preparation
1. **Bulk RNA-seq**: Place raw count matrices (CSV/TSV) and sample metadata in `data/bulk_rnaseq/PHENOTYPE/`.  
2. **GWAS Summary Statistics**: Store harmonized `.tsv` or `.txt` summary files in `data/gwas/PHENOTYPE/`.  
3. **Single-Cell Data**: Provide processed expression objects (`.rds` or `.h5ad`) and cell metadata in `data/sc_rnaseq/PHENOTYPE/`.  

Use consistent phenotype folder names (e.g., `Alzheimer`, `Parkinson`, `Schizophrenia`). Customize paths in each `config.yml` if needed.

---

## Workflow Overview
The pipeline consists of two main phases:

### 1. Disease-Specific Analysis
For each phenotype:
```bash
Rscript disease_specific/PHENOTYPE/preprocessing.R      # Data filtering & normalization
Rscript disease_specific/PHENOTYPE/differential_expression.R  # Identify DE genes
```
- **Output**: normalized matrices, DE gene tables, volcano plots in `disease_specific/PHENOTYPE/results/`.

### 2. Mega Analysis Across Phenotypes
Once all phenotypes finish:
```bash
Rscript mega_analysis/codes/integrate_rnaseq_gwas.R
Rscript mega_analysis/codes/sc_integration.R
Rscript mega_analysis/codes/summarize_results.R
```
- **Output**: cross-phenotype gene overlap tables, cell-type enrichment heatmaps, summary networks in `mega_analysis/results/`.

---

## Outputs & Interpretation
- **DE Gene Tables**: Fold changes, p-values, and adjusted p-values per phenotype
- **Volcano & MA Plots**: Quick visual check of expression shifts
- **Cross-Phenotype Overlap**: Venn diagrams or upset plots for shared genes
- **Cell-Type Mapping**: Dot plots showing phenotype-specific expression enriched in cell clusters
- **Integrated GWAS Links**: Correlation matrices linking expression signatures to GWAS hits

Refer to the generated HTML reports or RMarkdown outputs for interactive exploration.

---

## Customization & Extensions
- Edit individual `config.yml` files to tweak filtering thresholds, fold-change cutoffs, or annotation sources.
- Add new phenotype modules by copying an existing directory under `disease_specific/` and updating file paths.
- Extend `notebooks/` for bespoke visualizations or downstream analyses.

---

## Best Practices
- Keep raw data **immutable**; work only on copies in `results/`.
- Use **version control** branches for major modifications.
- Document any parameter changes in `config.yml` or as comments in the R scripts.

---

## Contributing
Contributions are welcome! Please:
1. Fork the repository  
2. Create a feature branch (`git checkout -b feature/your-feature`)  
3. Commit your changes (`git commit -m "Add new analysis module"`)  
4. Push and open a Pull Request  

See [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) for guidelines.

---

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
