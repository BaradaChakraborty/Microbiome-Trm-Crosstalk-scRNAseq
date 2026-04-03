# Microbiome Trm Crosstalk scATAC-seq

## Project Overview
This repository contains a complete, end-to-end single-cell ATAC-seq (Assay for Transposase-Accessible Chromatin) pipeline. The primary objective is to process raw epigenetic sequencing data to identify and characterize mucosal Tissue-Resident Memory T (Trm) cells and Secretory IgA Plasma cells within the small intestine microenvironment.

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

* **Script 05: Host-Microbe Trans-Kingdom Correlation**
    * Implements a statistical framework simulating microbial abundance (16S OTU reads) against host single-cell epigenetics.
    * Calculates Spearman's Rho to paramatize the crosstalk between bacterial load and physical *Itgae* chromatin accessibility.

* **Script 06: Secretory IgA Plasma Cell Epigenetic Mapping**
    * Computationally isolates the antibody-secreting B-cell compartment within the mucosal microenvironment.
    * Generates targeted coverage tracks mapping the euchromatin state of the *Igj* (J-chain) locus, confirming the presence of functional, gut-adapted Secretory IgA plasma cells.

## Visual Results
All generated plots proving data integrity and cellular identity are securely saved in the `results/figures/` directory. The coverage plots specifically demonstrate robust accessibility at the *Itgae* promoter for Trm cells, and the *Igj* promoter for Plasma cells, confirming the tissue-resident phenotype of the target clusters.

## Future Directions: In Vivo Validation Proposal
To physically validate the computational host-microbe crosstalk and IgA induction observed in this dataset, a controlled *in vivo* gnotobiotic experiment is proposed. Germ-free mice will be colonized with a defined microbial consortium (or challenged with a pathogen like *Salmonella*) to establish differential microbial loads, providing a living biological model of the simulated parametric correlation.

Following colonization kinetics, intestinal tissue microdissection and lamina propria cell isolation will be performed. High-parameter flow cytometry (FACS) will be utilized to physically sort the target immune populations identified in the scATAC-seq pipeline: specifically, the CD8+ CD103+ Trm compartment and the IgJ+ secretory plasma cell population.

Finally, the computational chromatin accessibility profiles will be validated using orthogonal wet-lab techniques. Utilizing confocal microscopy for spatial cluster mapping within the gut mucosa, alongside localized qRT-PCR, this approach will definitively link microbial metabolic pressure to the physical remodeling of the *Itgae* and *Igj* loci, bridging high-dimensional epigenetics with functional mucosal defense.

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
4. Execute the code in the `scripts/` directory in sequential order (01 through 06).