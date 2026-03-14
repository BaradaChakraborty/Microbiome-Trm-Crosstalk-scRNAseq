message("Step A: Calculating the Gene Activity Matrix")
gene.activities <- GeneActivity(si_trm_atac)

message("Step B: Adding the Gene Activity data to our T cell object")
si_trm_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)

message("Step C: Normalizing the new Gene Activity data")
DefaultAssay(si_trm_atac) <- 'RNA'
si_trm_atac <- NormalizeData(
  object = si_trm_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = 10000
)

message("Gene Activity Matrix successfully built! Ready to identify the clusters.")

library(ggplot2)

message("Step D: Visualizing specific Gut-Resident Memory Marker Genes")
marker_plot <- FeaturePlot(
  object = si_trm_atac,
  features = c("Cd8a", "Itgae"),
  max.cutoff = "q95",
  pt.size = 0.5,
  ncol = 2
)

message("Step E: Saving the Gene Activity plots")
ggsave(
  filename = "results/figures/03_scATAC_Marker_Genes.png",
  plot = marker_plot,
  width = 10, 
  height = 5, 
  dpi = 300
)

message("Marker plots perfectly saved! Check your results/figures folder.")