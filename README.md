# Microbiome Trm Crosstalk scATAC-seq

## Project Overview
This repository contains a complete, end-to-end single-cell ATAC-seq (Assay for Transposase-Accessible Chromatin) pipeline. The primary objective is to process raw epigenetic sequencing data to identify and characterize mucosal Tissue-Resident Memory T (Trm) cells within the small intestine microenvironment.

The workflow translates raw DNA fragment coordinates into biologically meaningful gene activity, culminating in the physical visualization of chromatin accessibility at key tissue-residency loci. The environment is properly configured for reproducible RStudio execution, utilizing standard UTF-8 encoding and two-space tabulation.

## Pipeline Architecture

* **Script 01: Epigenetic Quality Control**
    * Calculates Nucleosome Signal to assess the physical integrity of DNA wrapping.
    * Calculates Transcriptional Start Site (TSS) Enrichment to isolate living cells with high signal-to-noise ratios.
    * Filters out technical artifacts and dying cells.

* **Script 02: Dimensionality Reduction and Clustering**
    * Applies TF-IDF normalization to address the zero-inflated nature of scATAC-seq data.
    * Performs Singular Value Decomposition (SVD) and Latent Semantic Indexing (LSI).
    * Generates a 2D UMAP projection to identify distinct cellular neighborhoods.

* **Script 03: Gene Activity Matrix Generation**
    * Maps raw open chromatin peaks against the *Mus musculus* reference genome.
    * Estimates transcriptional activity to translate epigenetic coordinates into readable gene expression profiles.
    * Visualizes critical mucosal immunity markers (*Cd8a* and *Itgae*).

* **Script 04: Locus-Specific Coverage Visualization**
    * Maps physical transposase cut sites directly onto the genome.
    * Generates publication-ready coverage tracks proving chromatin unspooling at the *Itgae* (CD103) locus across distinct cellular clusters.

## Visual Results
All generated plots proving data integrity and cellular identity are securely saved in the `results/figures/` directory. The coverage plot specifically demonstrates robust accessibility at the *Itgae* promoter, confirming the tissue-resident phenotype of the target clusters.

## Dependencies & Reproducibility
This project utilizes the `renv` package manager to guarantee strict computational reproducibility. The primary R packages required for this epigenetic analysis include:

* **Signac & Seurat:** The core engines for single-cell chromatin and multimodal analysis.
* **EnsDb.Mmusculus.v79 & ensembldb:** Bioconductor databases providing the *Mus musculus* genomic architecture and mapping coordinates.
* **hdf5r & Rsamtools:** Essential dependencies for reading massive hierarchical binary matrices and indexing fragment files (`.tbi`).
* **ggplot2:** Utilized for generating all high-resolution, publication-ready data visualizations.

## How to Run the Pipeline
To replicate this analysis:
1. Clone this repository to your local machine.
2. Open the `Microbiome-Trm-Crosstalk-scRNAseq.Rproj` file in RStudio.
3. Run `renv::restore()` in the R console. This will automatically read the `renv.lock` file and install the exact package versions required.
4. Execute the code in the `scripts/` directory in sequential order (01 through 04).