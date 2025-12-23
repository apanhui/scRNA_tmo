
library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)

perform_pairwise_de_analysis <- function(obj, 
                                         cluster_col = "seurat_clusters",
                                         group_col = "group",
                                         ident.1 = NULL,
                                         ident.2 = NULL,
                                         min.pct = 0.1,
                                         logfc.threshold = 0.25,
                                         test.use = "wilcox",
                                         only.pos = FALSE,
                                         output_dir = "./DE_results") {
  
  # get info
  clusters <- unique(obj@meta.data[[cluster_col]])
  clusters <- sort(clusters)
  groups <- unique(obj@meta.data[[group_col]])
  groups <- sort(groups)
  group_pairs <- combn(groups, 2, simplify = FALSE)
  

  all_results <- list()
  # run Diff
  for (cluster in clusters) {
    cluster_cells <- Cells(obj)[obj@meta.data[[cluster_col]] == cluster]
    obj_cluster <- subset(obj, cells = cluster_cells)
    
    Idents(obj_cluster) <- obj_cluster@meta.data[[group_col]]
    for (pair in group_pairs) {
      group1 <- pair[1]
      group2 <- pair[2]
      
      cat(paste0(group1, " vs ", group2, "\n"))
      
      tryCatch({
        de_genes <- FindMarkers(object = obj_cluster,
          ident.1 = group1, ident.2 = group2,
          min.pct = min.pct,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          only.pos = only.pos,
          verbose = FALSE
        )
        
        if (nrow(de_genes) > 0) {
          de_genes$cluster <- cluster
          de_genes$group1 <- group1
          de_genes$group2 <- group2
          de_genes$comparison <- paste0(group1, "_vs_", group2)
          de_genes$gene <- rownames(de_genes)
          
          de_genes <- de_genes %>%
            select(cluster, comparison, group1, group2, gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, everything())
          
          # save result
          result_name <- paste0("cluster", cluster, "_", group1, "_vs_", group2)
          all_results[[result_name]] <- de_genes
          
          filename <- paste0("cluster", cluster, "_", group1, "_vs_", group2, ".tsv")
          write.table(de_genes, file = file.path(output_dir, filename), row.names = FALSE)
        } else {
          cat(paste0("No Differ\n"))
        }
        
      }, error = function(e) {
        cat(paste0("Error: ", conditionMessage(e), "\n"))
      })
    }
  }
  
  # combined all Differ
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    return(combined_results)
  } else {
    cat("\nNo Differ\n")
    return(NULL)
  }
}

### main
args <- commandArgs[T]
obj <- args[1]    # seurat object
outdir <- args[2]   # ouput directory

load(obj)
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

de_results <- perform_pairwise_de_analysis(
  obj = obj,
  cluster_col = "seurat_clusters",
  group_col = "group",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  only.pos = FALSE,
  output_dir = outdir
)
write.table(de_results, file = file.path(output_dir, "all_DE_genes_combined.tsv"), row.names = FALSE)
