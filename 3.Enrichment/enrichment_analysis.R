library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(RColorBrewer)
library(ggrepel)


# get diff glist
read_gene_lists <- function(directory) {
  gene_lists <- list()
  gene_files <- list.files(directory, pattern = "\\.glist$", full.names = TRUE)
  
  for (file in gene_files) {
    list_name <- tools::file_path_sans_ext(basename(file))
    genes <- readLines(file, header = FALSE, col.names = "gene")
    
    gene_vector <- unique(na.omit(genes$gene))
    gene_vector <- gene_vector[gene_vector != ""]
    
    if (length(gene_vector) > 0) {
      gene_lists[[list_name]] <- gene_vector
    }
  }
  
  return(gene_lists)
}

# read annot
read_custom_annotations <- function(go_file, kegg_file) {
  annotations <- list()
  
  # GO annot
  if (file.exists(go_file)) {
    go_annot <- read.table(go_file, header = FALSE, col.names = c("gene", "go_id", "description"))
    annotations$go <- go_annot
  } else {
    warning(paste("Can't not find",go_file))
  }
  
  # KEGG annot
  if (file.exists(kegg_file)) {
    kegg_annot <- read.table(kegg_file, header = FALSE, col.names = c("gene", "pathway_id", "description"))
    annotations$kegg <- kegg_annot
  } else {
    warning(paste("Can't not find", kegg_file))
  }
  
  return(annotations)
}

# do GO enrichment
custom_go_enrichment <- function(gene_list, go_annotations, background_genes = NULL, 
                                p_cutoff = 0.05, q_cutoff = 0.2, min_gs_size = 10, max_gs_size = 500) {
  
  if (is.null(go_annotations)) {
    stop()
  }
  
  term2gene <- go_annotations[, c("go_id", "gene")]
  term2name <- unique(go_annotations[, c("go_id", "description")])
  
  ego <- enricher(gene = gene_list,
                  pvalueCutoff = p_cutoff,
                  pAdjustMethod = "BH",
                  minGSSize = min_gs_size,
                  maxGSSize = max_gs_size,
                  qvalueCutoff = q_cutoff,
                  TERM2GENE = term2gene,
                  TERM2NAME = term2name)
  
  return(ego)
}

# do KEGG enrichment
custom_kegg_enrichment <- function(gene_list, kegg_annotations, background_genes = NULL,
                                 p_cutoff = 0.05, q_cutoff = 0.2, min_gs_size = 10, max_gs_size = 500) {
  
  if (is.null(kegg_annotations)) {
    stop()
  }
  
  term2gene <- kegg_annotations[, c("pathway_id", "gene")]
  term2name <- unique(kegg_annotations[, c("pathway_id", "description")])
  
  ekg <- enricher(gene = gene_list,
                  pvalueCutoff = p_cutoff,
                  pAdjustMethod = "BH",
                  minGSSize = min_gs_size,
                  maxGSSize = max_gs_size,
                  qvalueCutoff = q_cutoff,
                  TERM2GENE = term2gene,
                  TERM2NAME = term2name)
  
  return(ekg)
}

# summarize result
summarize_enrichment <- function(enrich_result, analysis_name, db_type) {
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    return(data.frame())
  }
  
  summary_df <- enrich_result@result %>%
    mutate(
      analysis = analysis_name,
      db = db_type,
      gene_count = as.numeric(str_extract(GeneRatio, "\\d+")),
      bg_count = as.numeric(str_extract(BgRatio, "\\d+")),
      gene_ratio = gene_count / Count,
      enrichment_factor = (gene_count / Count) / (bg_count / as.numeric(str_extract(BgRatio, "(?<=/)\\d+")))
    ) %>%
    select(ID, Description, analysis, db, pvalue, p.adjust, qvalue, Count, gene_count, 
           bg_count, gene_ratio, enrichment_factor, geneID)
  
  return(summary_df)
}

### main
args <- commandArgs[T]

if (length(args) < 5) {
  cat("Rscript enrichment_analysis.R <gene_dir> <go_file> <kegg_file> <output_dir>\n")
}

gene_dir <- args[1]
go_file <- args[2]
kegg_file <- args[3]
output_dir <- args[4]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1.read glist
gene_lists <- read_gene_lists(gene_dir)

# 2.read annot
custom_annot <- read_custom_annotations(go_file, kegg_file)

# 3.run enrichment
all_results <- list()
all_summaries <- data.frame()

for (list_name in names(gene_lists)) {
  gene_list <- gene_lists[[list_name]]
  # GO
  if (!is.null(custom_annot$go)) {
    go_result <- custom_go_enrichment(gene_list, custom_annot$go)
    
    if (!is.null(go_result) && nrow(go_result@result) > 0) {
      all_results[[paste0(list_name, "_GO")]] <- go_result
      go_summary <- summarize_enrichment(go_result, list_name, "GO")
      all_summaries <- rbind(all_summaries, go_summary)
      write.csv(go_result@result, file.path(output_dir, "tables", paste0(list_name, "_GO_enrichment.csv")), row.names = FALSE)
    } else {
      cat("No Diff\n")
    }
  }
  
  # KEGG
  if (!is.null(custom_annot$kegg)) {
    kegg_result <- custom_kegg_enrichment(gene_list, custom_annot$kegg)
    
    if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
      all_results[[paste0(list_name, "_KEGG")]] <- kegg_result
      kegg_summary <- summarize_enrichment(kegg_result, list_name, "KEGG")
      all_summaries <- rbind(all_summaries, kegg_summary)
      write.csv(kegg_result@result, file.path(output_dir, "tables", paste0(list_name, "_KEGG_enrichment.csv")), row.names = FALSE)
    } else {
      cat("No Diff\n")
    }
  }
}

