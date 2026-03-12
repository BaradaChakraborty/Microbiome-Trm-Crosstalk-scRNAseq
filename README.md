Microbiome-Trm-Crosstalk-scRNAseq
Author: Barada Chakraborty

Affiliation: Independent Researcher
Project Overview
This repository contains the computational pipeline for mapping microbiota-dependent cytotoxic Tissue-Resident Memory (Trm) cell programming using single-cell RNA sequencing and bulk transcriptomics.

Biological Rationale
The intestinal microbiota actively dictates the transcriptional programming required for circulating effector T cells to differentiate into protective mucosal Trm cells. This project investigates how transient depletion or alteration of the microbiota reshapes the transcriptomic trajectory of these T cells, either impairing their tissue residency markers (such as CD103 and CD69) or skewing their functional cytokine profiles. This aligns directly with contemporary research on mucosal immunity and host-microbe bidirectional crosstalk.

Computational Workflow
Data Acquisition: Mining the Gene Expression Omnibus (GEO) for scRNA-seq datasets of intestinal lamina propria lymphocytes from wild-type versus germ-free or antibiotic-treated murine models.

Single-Cell Processing: Utilizing Seurat for quality control, filtering, and clustering to isolate cytotoxic T cell populations.

Marker Identification: Identifying Trm-specific signatures, specifically the upregulation of residency genes (Itgae, Cd69, Prdm1) and the downregulation of tissue-egress receptors (S1pr1, Sell).

Trajectory Inference: Applying Monocle3 to map the pseudotime developmental trajectory from circulating effector cells into fully differentiated mucosal Trm cells.

Integration: Running DESeq2 on mucosal bulk RNA-seq data to calculate differential gene expression of broad immune signaling pathways (e.g., TGF-beta signaling).

Verified Literature Context
This computational approach is grounded in leading mucosal immunology literature, specifically focusing on the mechanisms by which intestinal microbes undergo rapid transcriptional and metabolic adaptation to host immune activation, and how transient microbiota depletion enhances mucosal immunity.
