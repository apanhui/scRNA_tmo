library(ggfun)
library(ggrepel)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tibble)

args <- commandArgs(T)

fin <- args[1]
flabel <- args[2]
outpre <- args[3]

dd1 <- read.table(fin, sep = '\t', header = T, check.names = F, quote = "", stringsAsFactors = F)
dd1$Group <- factor(dd1$Group, levels = c("EC","VC", "MC"))
dd1 <- dd1[order(factor(dd1$significance, levels = c("nosig", "up", "down"))), ]
str(dd1)

label1 <- read.table(flabel, sep = '\t', header = T, check.names = F)

# get jitter pos
set.seed(122)

mycolor <- c("#a495c2","#f8ee53","#79c360") 
p <- ggplot(dd1, aes(x = Group, y = logFC)) +
     geom_jitter(mapping = aes(color = significance), size = .5, width = 0.45) +
     geom_jitter(label1, mapping = aes(x = Group, y = logFC), size = 2, width = 0.45, color = 'black') +
     geom_text_repel(label1, mapping = aes(label = Gene), size = 4, box.padding = 0, point.padding = 0, min.segment.length = 0, max.overlaps = Inf, family = 'sans', position = position_jitter(seed = 123)) + 
     scale_color_manual(values = c(up = '#d42525', nosig = '#bebebe', down = "#3186bd")) +
     labs(y = 'log2FC', x = 'Cluster', color = 'Significant') +
     guides(color = guide_legend(override.aes = list(size = 3))) +
     geom_tile(dd1 %>% dplyr::select(Group) %>% dplyr::distinct(Group), mapping = aes(x = Group, y = 0), height = 0.6, color = mycolor, fill = mycolor, alpha = 1, show.legend = F) +
     theme_bw() +
     theme(axis.text = element_text(color = "#000000", size = 12), axis.title = element_text(color = "#000000", size = 15))

pdfout <- paste0(outpre,".pdf")
pngout <- paste0(outpre,".png")
ggsave(p, file = pdfout, width = 9, height = 8)
system(paste("convert -density 300", pdfout, pngout))
