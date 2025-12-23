library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

args <- commandArgs(T)
expr <- read.table(args[1], sep = "\t", quote = "", header = T, check.names = F, stringsAsFactors = F)
cluster_df <- read.table(args[2], sep = "\t", quote = "", header = F, check.names = F, stringsAsFactors = F)
outdir <- args[3]

colnames(cluster_df) <- c('Cells', 'Cluster')

data_long <- melt(expr, variable.name = "Gene") %>% left_join(cluster_df, by="Cells")
head(data_long)
data_long$Cluster <- factor(data_long$Cluster, levels = c("leaves_03","flowers_03"))

get_sig <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}


## each gene
results_df <- data.frame(gene = character(), p_value = numeric(), signif = character(), stringsAsFactors = FALSE)

for (current_gene in unique(data_long$Gene)) {
  gene_data <- data_long %>% filter(Gene == current_gene)
  test_result <- wilcox.test(value ~ Cluster, data = gene_data)
  p_value <- test_result$p.value
  sig <- get_sig(p_value)
  results_df <- rbind(results_df, data.frame(Gene = current_gene, p_value = p_value, signif = sig))

  p <- ggplot(gene_data, aes(x=Cluster, y=value, fill=Cluster)) +
    geom_violin(scale="width", trim=TRUE) +
    #ylim(0, max(gene_data$value)/3) +
    scale_fill_manual(values=c("leaves_03"="#50c94f", "flowers_03"="#d681ec")) +
    labs(x="Cluster", y="Expression", title=current_gene) +
    theme_bw() +
    stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("leaves_03", "flowers_03")),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***","**","*", "ns"))) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
  
  ggsave(p, filename=paste0(outdir,"/", current_gene, "_violin.pdf"), width=6, height=5)
  ggsave(p, filename=paste0(outdir,"/", current_gene, "_violin.png"), width=6, height=5, dpi=300)
}
write.table(results_df, file = 'diff.xls', sep = "\t", quote = F, col.names = T, row.names = F)

#results_df <- read.table('diff.xls', sep = "\t", quote = "", header = T, check.names = F, stringsAsFactors = F)
# all genes
data_long <- left_join(data_long, results_df, by='Gene')
plot <- ggplot(data_long, aes(x=Gene, y=value, fill=Cluster)) +
  geom_violin(scale="width", trim=TRUE, position="dodge") +
  scale_fill_manual(values=c("leaves_03"="#50c94f", "flowers_03"="#d681ec")) +
  #geom_text(aes(x = Gene, y = max(data_long$value), label = signif), family = 'sans') +
  geom_text(aes(x = Gene, y = 20, label = signif), family = 'sans') +
  ylim(0, 20) +
  #coord_cartesian(ylim = c(0, 20)) +
  labs(x="Gene", y="Expression", title="") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())

ggsave(paste0(outdir, "/all_violin.pdf"), plot, width=14, height=6)
ggsave(paste0(outdir, "/all_violin.png"), plot, width=14, height=6, dpi=300)
