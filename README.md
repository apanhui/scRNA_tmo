## Pre-install software and R packages

### Software

+ CellRanger 8.0.0
+ R-4.3.2

### R packages

+ Seurat 5.1.0
+ DoubletFinder-2.0.3
+ harmony 1.2.0
+ ggplot2_3.3.5
+ patchwork_1.1.1
+ dplyr_1.0.9

## Scripts

### 1. CellRanger & Seurat

1.cellRanger_seurat/work.sh -- the main shell for run Cellranger and Seurat.

Main content

+ CellRanger analysis
+ DoubletFinder analysis
+ Filter Cells
+ Normalization
+ PCA
+ Reduce dimension (Umap and tSNE)
+ RunHarmony for Batch correction
+ DoFindClusters
+ Draw figures
+ Find marker genes
+ Plot marker genes

### 2. plot

2.plot/work.sh -- the main shell for plot

Main content

+ dotplot plot
+ violin plot