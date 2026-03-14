library(Signac)
library(Seurat)
library(ggplot2)

message("Step A: Switching back to the Epigenetic DNA tab")
DefaultAssay(si_trm_atac) <- "peaks"

message("Step B: Generating the epigenetic mountain ranges for the Itgae (CD103) gene")
cov_plot <- CoveragePlot(
  object = si_trm_atac,
  region = "Itgae",
  extend.upstream = 5000,
  extend.downstream = 5000
)

message("Step C: Saving the high-resolution Coverage Plot")
ggsave(
  filename = "results/figures/04_scATAC_Itgae_Coverage.png",
  plot = cov_plot,
  width = 8, 
  height = 10, 
  dpi = 300
)

message("Coverage Plot perfectly saved! Check your results/figures folder.")