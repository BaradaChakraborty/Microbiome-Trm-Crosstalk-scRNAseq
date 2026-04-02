library(Signac)
library(Seurat)
library(ggplot2)

message("Step A: Switching to the Gene Activity (RNA) tab")
DefaultAssay(si_trm_atac) <- "RNA"

message("Step B: Visualizing B-cell and IgA Plasma Cell Markers")
# Cd19 is the classic marker for B-cells. 
# Igha is the actual blueprint for the IgA antibody.
iga_feature_plot <- FeaturePlot(
  object = si_trm_atac,
  features = c("Cd19", "Igha"),
  max.cutoff = "q95",
  pt.size = 0.5,
  ncol = 2
)

message("Step C: Saving the IgA Marker Plot")
ggsave(
  filename = "results/figures/06_scATAC_IgA_Markers.png", 
  plot = iga_feature_plot, 
  width = 10, 
  height = 5, 
  dpi = 300
)

message("Step D: Switching back to the physical DNA tab")
# We must switch back to the raw DNA to draw the coverage mountains
DefaultAssay(si_trm_atac) <- "peaks"

message("Step E: Drawing the epigenetic mountains for the Igj gene (IgA Secreting Cells)")
iga_cov_plot <- CoveragePlot(
  object = si_trm_atac,
  region = "Igj",
  extend.upstream = 5000,
  extend.downstream = 5000
)

message("Step F: Saving the Plasma Cell Coverage Plot")
ggsave(
  filename = "results/figures/07_scATAC_Igj_Coverage.png", 
  plot = iga_cov_plot, 
  width = 8, 
  height = 10, 
  dpi = 300
)

message("IgA Pipeline Complete! Check your results/figures folder for the two new plots.")