
args <- commandArgs(T)
if(length(args) < 3) {
  stop("Rscript-3.5.1_conda xxx.R <file> <glist_file> <outdir> [dotplot,violin] [seurat_clusters|orig.ident|Groups] [group_order]\n")
}

file   <- args[1]
glist  <- args[2]
outdir <- args[3]
type   <- ifelse(length(args) > 3, args[4], 'feature')
group  <- ifelse(length(args) > 4, args[5], 'seurat_clusters')
group_ord <- ifelse(length(args) > 5, args[6], NA)

run <- list(
  dotplot = grepl('dotplot', type),
  violin  = grepl('violin',  type),
  group.by = group
)
print(run)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

source("Seurat_lib.R")

## read gene list
if (file.exists(glist)) {
  genes <- unique(readLines(glist))
  genes <- gsub('_', '-', genes)
  print(genes)
} else {
  stop(paste(glist, 'not found!'))
}

### Creat Seurat Object
message( "==>Reading 10x data<==" )
load(file)
if (!'obj' %in% ls()) {
  obj <- object
  rm('object')
}
if (!is.na(group_ord)) {
  group_ord <- readLines(group_ord)
  obj <- obj[, as.character(obj@meta.data[[group]]) %in% group_ord]
  obj@meta.data[[group]] <- factor(obj@meta.data[[group]], levels = group_ord)
  print(str(obj@meta.data[[group]]))
}

### set work dir
system(paste("mkdir -p", outdir))
if ( ! is.na(outdir) ) setwd(outdir)

### input glist
if (glist != 'none') {
  message( "==>display markers<==" )
  print(str(genes))
  outpre <- 'Marker'
} else {
  q()
}

### dotplot
if (run$dotplot) {
  message( "==> dotplot <==" )
  pdf_out <- paste0(outpre, ".DotPlot.pdf")
  png_out <- paste0(outpre, ".DotPlot.png")
  p <- PlotDotPlot(obj, features = unique(genes), group.by = run$group.by, outfile = NULL)
  p <- p + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')))
  w <- max(6, length(unique(obj@meta.data[[run$group.by]])) * 0.4 + 2)
  h <- max(6, ceiling(length(genes)) * 0.35 + 2)
  #p <- p + coord_flip()
  ggsave(p, file = pdf_out, width = h, height = w)
  system(paste("convert -density 300", pdf_out, png_out))
}

### violin plot
if (run$violin) {
  message( "==> violin plot<==" )
  dir.create("ViolinPlot/", showWarnings = F, recursive = T)
  unlink("ViolinPlot/*", recursive = T)
  PlotVlnPlot(obj, unique(genes), outpref = "ViolinPlot/ViolinPlot", group.by = run$group.by)
  system("for i in ViolinPlot/*.pdf; do convert -density 200 $i ${i/pdf/png}; done")
}

### Hasta la vista, baby
system("rm -rf Rplot.pdf")
message( "==>All Done!<==" )

