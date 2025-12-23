library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(fgsea)
library(msigdbr)

# prepare gmt
read_gmt_file <- function(gmt_file) {
  gmt_data <- readLines(gmt_file)
  
  if (!all(c("term", "gene") %in% colnames(gmt_data))) {
    stop()
  }
  
  gene_sets <- split(gmt_data$gene, gmt_data$term)
  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
  return(gene_sets)
}

# order genes
prepare_gene_rank <- function(deg_table, rank_column = "log2FC") {
  deg_table <- deg_table[!is.na(deg_table[[rank_column]]), ]
  deg_table <- deg_table[!is.na(deg_table$gene), ]
  
  deg_table <- deg_table %>%
    group_by(gene) %>%
    arrange(desc(abs(.data[[rank_column]]))) %>%
    slice(1) %>%
    ungroup()
  
  gene_rank <- setNames(deg_table[[rank_column]], deg_table$gene)
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  return(gene_rank)
}

# GSEA
run_gsea_analysis <- function(gene_rank, gene_sets, gmt_name = "gene_set", 
                              min_size = 15, max_size = 500, n_perm = 10000) {
  
  # filter genes
  gs_sizes <- sapply(gene_sets, length)
  valid_sets <- names(gene_sets)[gs_sizes >= min_size & gs_sizes <= max_size]
  gene_sets_filtered <- gene_sets[valid_sets]
  
  gsea_result <- GSEA(gene_rank,
                      TERM2GENE = data.frame(
                        term = rep(names(gene_sets_filtered), 
                                   sapply(gene_sets_filtered, length)),
                        gene = unlist(gene_sets_filtered)
                      ),
                      nPerm = n_perm,
                      minGSSize = min_size,
                      maxGSSize = max_size,
                      pvalueCutoff = 0.05,
                      verbose = FALSE)
  
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    warning("No Result")
    return(NULL)
  }
  
  return(gsea_result)
}

# format GSEA results
format_gsea_results <- function(gsea_result, gmt_name = "gene_set") {
  if (is.null(gsea_result)) {
    return(data.frame())
  }
  
  result_df <- gsea_result@result %>%
    mutate(
      database = gmt_name,
      set_size = as.numeric(setSize),
      enrichment_score = as.numeric(NES),
      normalized_enrichment_score = as.numeric(NES),
      p_value = as.numeric(pvalue),
      p_adjust = as.numeric(p.adjust),
      q_value = as.numeric(qvalues),
      leading_edge_size = as.numeric(core_enrichment %>% 
                                       sapply(function(x) length(strsplit(x, "/")[[1]]))),
      leading_edge_percent = (leading_edge_size / set_size) * 100
    ) %>%
    select(
      ID, Description, database, set_size, enrichment_score, normalized_enrichment_score,
      p_value, p_adjust, q_value, 
      leading_edge_size, leading_edge_percent,
      core_enrichment, everything()
    ) %>%
    arrange(p_adjust, p_value)
  
  return(result_df)
}

# save result
save_gsea_results <- function(gsea_results, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_results <- data.frame()
  
  for (i in seq_along(gsea_results)) {
    gmt_name <- names(gsea_results)[i]
    result_df <- gsea_results[[i]]
    
    if (nrow(result_df) > 0) {
      output_file <- file.path(output_dir, paste0("GSEA_", gmt_name, "_detailed.csv"))
      write.csv(result_df, output_file, row.names = FALSE)
      
      sig_results <- result_df %>%
        filter(p_adjust < 0.25) %>%
        arrange(p_adjust, desc(abs(normalized_enrichment_score)))
      
      if (nrow(sig_results) > 0) {
        sig_file <- file.path(output_dir, paste0("GSEA_", gmt_name, "_significant_FDR25.csv"))
        write.csv(sig_results, sig_file, row.names = FALSE)
      } else {
        cat("No sig diff\n")
      }
      
      all_results <- bind_rows(all_results, result_df)
    }
  }
  
  if (nrow(all_results) > 0) {
    combined_file <- file.path(output_dir, "GSEA_all_results_combined.csv")
    write.csv(all_results, combined_file, row.names = FALSE)
  }
  
  return(all_results)
}


### main
args <- commandArgs[T]

if (length(args) < 4) {
  cat("Rscript gsea_analysis.R <deg_file> <go_gmt> <kegg_gmt> <output_dir> [rank_column]\n")
}

deg_file <- args[1]
go_gmt_file <- args[2]
kegg_gmt_file <- args[3]
output_dir <- args[4]
rank_column <- if (length(args) > 4) args[5] else "log2FC"

# 1.read deg_file
if (file.exists(deg_file)) {
  deg_data <- read.table(deg_file, sep = "\t", quote = "", header = T)
  # check input
  if (!"gene" %in% colnames(deg_data)) {
    possible_gene_cols <- c("Gene", "gene_name", "gene_id", "GeneID", "GeneSymbol")
    gene_col <- intersect(possible_gene_cols, colnames(deg_data))
    if (length(gene_col) > 0) {
      deg_data$gene <- deg_data[[gene_col[1]]]
    } else {
      stop("Not find gene col")
    }
  }
} else {
  stop("Not find input:", deg_file)
}

# 2.prepare rank
gene_rank <- prepare_gene_rank(deg_data, rank_column)

# 3.read gmt
gmt_files <- c(GO = go_gmt_file, KEGG = kegg_gmt_file)
gene_sets_list <- list()

for (gmt_name in names(gmt_files)) {
  gmt_file <- gmt_files[[gmt_name]]
  if (file.exists(gmt_file)) {
    gene_sets_list[[gmt_name]] <- read_gmt_file(gmt_file)
  } else {
    warning("Not find:", gmt_file)
  }
}

# 4.run GSEA
gsea_results_list <- list()

for (gmt_name in names(gene_sets_list)) {
  gsea_result <- run_gsea_analysis(gene_rank, 
                                   gene_sets_list[[gmt_name]], 
                                   gmt_name = gmt_name)
  
  if (!is.null(gsea_result)) {
    formatted_results <- format_gsea_results(gsea_result, gmt_name)
    gsea_results_list[[gmt_name]] <- formatted_results
  } else {
    cat(gmt_name, "No result!\n")
  }
}

all_results <- save_gsea_results(gsea_results_list, output_dir)

