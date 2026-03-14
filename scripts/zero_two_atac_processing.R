library(Signac)
library(Seurat)
library(ggplot2)

message("Step A: Performing TF-IDF Normalization to balance sparse DNA data")
si_trm_atac <- RunTFIDF(si_trm_atac)

message("Step B: Finding the most informative DNA peaks")
si_trm_atac <- FindTopFeatures(si_trm_atac, min.cutoff = 'q0')

message("Step C: Running Singular Value Decomposition (SVD)")
si_trm_atac <- RunSVD(si_trm_atac)

message("Step D: Generating UMAP and Identifying Cell Clusters")
si_trm_atac <- RunUMAP(object = si_trm_atac, reduction = 'lsi', dims = 2:30)
si_trm_atac <- FindNeighbors(object = si_trm_atac, reduction = 'lsi', dims = 2:30)
si_trm_atac <- FindClusters(object = si_trm_atac, verbose = FALSE, algorithm = 3)

message("Data processing complete! The mucosal T cells are now mathematically clustered.")

message("Step E: Visualizing the T Cell Clusters on a 2D Map")
umap_plot <- DimPlot(object = si_trm_atac, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("scATAC-seq UMAP: Small Intestine Trm Cells")

message("Step F: Saving the high-resolution UMAP plot")
ggsave(
  filename = "results/figures/02_scATAC_UMAP_Clusters.png",
  plot = umap_plot,
  width = 8, 
  height = 6, 
  dpi = 300
)

message("UMAP Plot perfectly saved! Check your results/figures folder.")