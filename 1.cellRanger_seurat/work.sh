### Step1 cellranger
ref_indir=''
fq_indir=''
# cellranger-8.0.0
cellranger count --id flowers-a --transcriptome ${ref_indir}/ --disable-ui --fastqs ${fq_indir}/flowers-a --create-bam=true --sample flowers-a --expect-cells 3000  --chemistry threeprime --include-introns true --localcores 4 --localmem 100
cellranger count --id flowers-b --transcriptome ${ref_indir}/ --disable-ui --fastqs ${fq_indir}/flowers-b --create-bam=true --sample flowers-b --expect-cells 3000  --chemistry threeprime --include-introns true --localcores 4 --localmem 100
cellranger count --id leaves-a --transcriptome ${ref_indir}/ --disable-ui --fastqs ${fq_indir}/leaves-a --create-bam=true --sample leaves-a --expect-cells 3000  --chemistry threeprime --include-introns true --localcores 4 --localmem 100
cellranger count --id leaves-b --transcriptome ${ref_indir}/ --disable-ui --fastqs ${fq_indir}/leaves-b --create-bam=true --sample leaves-b --expect-cells 3000  --chemistry threeprime --include-introns true --localcores 4 --localmem 100

### Step2 DoubleFinder
indir=''
outdir=''
Rscript Doublet.R ${indir}/flowers-a flowers-a Seurat_lib.R ${outdir}/flowers-a
Rscript Doublet.R ${indir}/flowers-b flowers-b Seurat_lib.R ${outdir}/flowers-b
Rscript Doublet.R ${indir}/leaves-a leaves-a Seurat_lib.R ${outdir}/leaves-a
Rscript Doublet.R ${indir}/leaves-b leaves-b Seurat_lib.R ${outdir}/leaves-b

### Step3 Seurat
outdir=''
Rscript Seurat.R parameter.yaml ${outdir} Seurat_lib.R

