library(Signac)
library(Seurat)
library(hdf5r)

message("Step A: Loading the DNA peak matrix for Small Intestine Trm cells")
counts <- Read10X_h5("data/raw/GSM6214535_4_SI_IEL_filtered_peak_bc_matrix.h5")

message("Step B: Creating the Epigenetic Chromatin Assay")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "data/raw/GSM6214535_4_SI_IEL_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

message("Step C: Initializing the scATAC-seq Seurat Object")
si_trm_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = "SI_IEL"
)

message("Object created successfully! Ready for epigenetic quality control.")
si_trm_atac

library(EnsDb.Mmusculus.v79)
library(ensembldb)

message("Step D: Extracting the official Mus musculus genome annotations")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

message("Step E: Formatting chromosome names to match our mucosal dataset (e.g., changing '1' to 'chr1')")
seqlevelsStyle(annotations) <- "UCSC"

message("Step F: Attaching the genomic street map directly to our T cell object")
Annotation(si_trm_atac) <- annotations

message("Genome successfully attached! The object now knows exactly where the immunity genes are located.")


message("Step G: Calculating Nucleosome Signal (Checking for healthy DNA wrapping)")
si_trm_atac <- NucleosomeSignal(object = si_trm_atac)

message("Step H: Calculating TSS Enrichment (Using the FAST approximation method)")
si_trm_atac <- TSSEnrichment(object = si_trm_atac, fast = TRUE)

message("Step I: Filtering out dead cells and technical noise")
si_trm_atac <- subset(
  x = si_trm_atac,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

message("Epigenetic Quality Control Complete! Ready for visualization.")
si_trm_atac

library(ggplot2)

message("Step J: Creating the figures directory inside the results folder")
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

message("Step K: Generating the Quality Control Violin Plot")
qc_plot <- VlnPlot(
  object = si_trm_atac,
  features = c("nCount_peaks", "nucleosome_signal", "TSS.enrichment"),
  pt.size = 0.05,
  ncol = 3
)

message("Step L: Saving the high-resolution plot")
ggsave(
  filename = "results/figures/01_scATAC_QC_Metrics.png",
  plot = qc_plot,
  width = 12, 
  height = 5, 
  dpi = 300
)

message("QC Plot perfectly saved! Check your results/figures folder.")