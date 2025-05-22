# Brain Meta-Analysis

## Overview

This repository provides a comprehensive meta-analysis of eight brain-related phenotypes using bulk RNA-seq, GWAS summary statistics, and single-cell transcriptomic data. The analysis is structured into two key phases:

1. **Disease-Specific Analysis**: Individual pipelines for each brain phenotype in `disease_specific/`.
2. **Cross-Phenotype Mega-Analysis**: Aggregation and comparative analysis across all phenotypes in `mega_analysis/codes/`.

## Table of Contents

- [Installation](#installation)
- [Data Requirements](#data-requirements)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
  - [1. Disease-Specific Analysis](#1-disease-specific-analysis)
  - [2. Cross-Phenotype Analysis](#2-cross-phenotype-analysis)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

## Installation

Ensure you have the following installed:

- R (>= 4.1.0)
- R packages: `tidyverse`, `data.table`, `limma`, `edgeR`, `DESeq2`, `singleCellExperiment`, etc.
- Optional: `renv` for project-local dependency management

```bash
git clone https://github.com/isadeghi87/Brain_meta_analysis.git
cd Brain_meta_analysis
# (Optional) Initialize renv
Rscript -e "install.packages('renv'); renv::init()"
Rscript -e "renv::restore()"
```

## Data Requirements

Prior to running the analysis, download and organize your data as follows:

- **Bulk RNA-seq**: Place raw counts and metadata in `data/bulk_rnaseq/PHENOTYPE/`
- **GWAS Summary Statistics**: Place summary `.txt` or `.tsv` files in `data/gwas/PHENOTYPE/`
- **Single-Cell RNA-seq**: Place expression matrices (e.g., `.h5ad`, `.rds`) and cell metadata in `data/sc_rnaseq/PHENOTYPE/`

Replace `PHENOTYPE` with the actual phenotype name (e.g., `Alzheimer`, `Parkinson`, etc.).

## Directory Structure

```text
├── disease_specific/          # Phenotype-specific analysis pipelines
│   ├── Alzheimer/
│   │   ├── preprocessing.R
│   │   ├── differential_expression.R
│   │   └── results/
│   ├── Parkinson/
│   └── ...                   # Other phenotypes
├── mega_analysis/             # Cross-phenotype analysis codes
│   └── codes/
│       ├── integrate_rnaseq_gwas.R
│       ├── sc_integration.R
│       └── summarize_results.R
├── data/                      # Symlink or placeholders for input data directories
│   ├── bulk_rnaseq/
│   ├── gwas/
│   └── sc_rnaseq/
└── README.md                  # This file
```

## Usage

### 1. Disease-Specific Analysis

For each phenotype, run the provided R scripts:

```bash
# Example for Alzheimer's analysis
Rscript disease_specific/Alzheimer/preprocessing.R
Rscript disease_specific/Alzheimer/differential_expression.R
```

Outputs are saved in `disease_specific/Alzheimer/results/`.

### 2. Cross-Phenotype Analysis

Once all disease-specific pipelines have run, perform the mega-analysis:

```bash
Rscript mega_analysis/codes/integrate_rnaseq_gwas.R
Rscript mega_analysis/codes/sc_integration.R
Rscript mega_analysis/codes/summarize_results.R
```

Results will be stored in `mega_analysis/results/`.

## Results

- **Volcano plots**, **heatmaps**, and **network analyses** can be found in the respective `results/` directories.
- Summary tables of overlapping differentially expressed genes and GWAS hits are provided as `.csv`.

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/YourFeature`
3. Commit your changes: `git commit -m 'Add feature'`
4. Push to the branch: `git push origin feature/YourFeature`
5. Open a Pull Request.

Please follow the [Contributor Covenant](CODE_OF_CONDUCT.md).

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
