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