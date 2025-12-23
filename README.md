# Pre-install software and R packages

## Software

``` 
CellRanger	8.0.0
R	4.3.2
python	3.9.12
```

## Rpackages

	clusterProfiler	3.14.3
	cowplot	1.1.1
	data.table	1.14.2
	DOSE	3.12.0
	DoubletFinder	2.0.4
	dplyr	1.0.9
	enrichplot	1.6.1
	fgsea	1.12.0
	ggfun	0.0.9
	ggplot2	3.3.5
	ggpubr	0.4.0
	ggrepel	0.8.2
	harmony	1.2.0
	igraph	1.3.0
	magrittr	2.0.1
	monocle	2.22.0
	patchwork	1.1.1
	purrr	0.3.4
	RColorBrewer	1.1.2
	reshape2	1.4.3
	Seurat	5.1.0
	stringr	1.4.0
	tibble	3.2.1
	tidyr	1.1.3
# Scripts

## 1. scRNAseq script

- Package require:

  ```
  Seurat, DoubletFinder, dplyr, ggplot2, patchwork, future
  ```

- Step1 Cellranger

  ``` 
  ref_indir=''
  sample=''
  fq_indir=''
  cellranger count --id flowers-a --transcriptome ${ref_indir}/ --disable-ui --fastqs ${fq_indir}/${sample} --create-bam=true --sample ${sample} --expect-cells 3000  --chemistry threeprime --include-introns true --localcores 4 --localmem 100
  ```

- Step2 DoubleFinder

  ```
  indir='‘		# Triplet File
  sample=''
  outdir=''
  Rscript Doublet.R ${indir}/ ${sample} Seurat_lib.R ${outdir}/
  ```

- Step3 Seurat

  ```
  outdir=''
  Rscript Seurat.R parameter.yaml ${outdir}/ Seurat_lib.R
  ```

- Step4 plot

  ``` 
  input='‘		# seurat_object
  glist=''		# target gene list
  outdir=''
  Rscript seurat_plot.R ${input} ${glist} ${outdir} [dotplot,violin] [seurat_clusters|orig.ident|Groups] [group_order]
  ```

## 2. scRNAseq Diff script

- Package require: 

  ```
  Seurat, dplyr, ggfun, ggplot2, ggpubr, ggrepel, magrittr, purrr, RColorBrewer, reshape2, tibble, tidyr
  ```

- Step1 Group Diff 

  ``` 
  input = '‘		# seurat_object
  outdir = ''
  Rscript scDiff.R ${input} ${outdir}/
  ```

- Step2 Plot

  ```
  # Violin
  Rscript-3.6.3 plot.R plot.violin.draw.xls plot.violin.cell.list outdir/
  # Volcano
  Rscript bar.volcano.R bar.volcano.draw.xls bar.volcano.show.names.xls outprefx
  ```

## 3. Enrichment

- Package require: 

  ```
  clusterProfiler, data.table, DOSE, dplyr, enrichplot, fgsea, ggplot2, ggrepel,  patchwork, RColorBrewer, stringr, tidyr
  ```

- Step1 GO/KEGG enrichment

  ```
  outdir = '' 
  Rscript enrichment_analysis.R example.degene.glist example.go.annot example.kegg.annot ${outdir}/
  ```

- Step2 GSEA

  ```
  $input = ''		# DifferMarker with log2FC and pvalue
  outdir = ''
  Rscript GSEA_analysis.R ${input} go.gmt kegg.gmt ${outdir}/
  ```

- Step3 Plot

  ``` 
  # dotplot
  Rscript plot_dot.R example.enrichment.draw.xls outprefx 20 "Gene" "GOterm" "Top 20 of GO Enrichment" 0 all
  # barplot
  Rscript plot_bar.R example.enrichment.draw.xls outprefx 20 "Gene" "Pathway" "Top 20 of KEGG Enrichment" all "none"
  # GSEA plot
  Rscript plot_GSEA.r example.GSEA.draw.xls example.GSEA.ES.xls example.GSEA.ranked_gene.xls KO00790 "KO00790	Folate biosynthesis" output.pdf no
  ```

## 4. monocle2 scricpt

- Package require: 

  ```
  dplyr, igraph, monocle, Seurat, cowplot
  ```

- Step1 monole2

  ```
  input = '‘		# seurat_object
  outdir = ''
  Rscript monocle2.R  ${input} ${outdir}/
  ```

  

