library(ggplot2)

args = commandArgs(TRUE)

term_info = args[1]
all_info = args[2]
rnk_info = args[3]
name = args[4]
title = args[5]
outpdf = args[6]
slope = args[7]


replotGSEA <- function(term_info, all_info, rnk_info, gene.set, class.name, metric.range,
                       enrichment.score.range,png=png,pdf=pdf) {
  cex_value=2
  if (missing(gene.set)) {
    stop("Gene set argument is required")
  }

  ## Load .rnk data
  gsea.rnk <- read.delim(file = rnk_info, header = TRUE)
  gsea.rnk = gsea.rnk[,c("NAME","SCORE")]
  colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
  metric_list = gsea.rnk$metric[!is.na(gsea.rnk$metric) & gsea.rnk$metric!=-Inf & gsea.rnk$metric!=Inf]
  if (missing(metric.range)) {
     metric.range <- c(min(metric_list), max(metric_list))
  }  
  
  ## read es info 
  all_info_data = read.delim(file = all_info, header = TRUE, stringsAsFactors = FALSE)
  # Get enrichment score
  gsea.enrichment.score <- as.numeric(all_info_data[all_info_data[,1]==gene.set,"ES"])
  color = ifelse(gsea.enrichment.score>0,"red","green")
  
  # Get gene set name
  gsea.normalized.enrichment.score <- as.numeric(all_info_data[all_info_data[,1]==gene.set,"NES"])

  # Get nominal p-value
  gsea.p.value <- as.numeric(all_info_data[all_info_data[,1]==gene.set,"NOM.p.val"])
  
  # Get FDR
  gsea.fdr <- as.numeric(all_info_data[all_info_data[,1]==gene.set,"FDR.q.val"])

  # gsea.metric
  gsea.metric = "SCORE"

  ## term infos 
  term_data = read.delim(file = term_info, header = TRUE, stringsAsFactors = FALSE)
  
  # Get hit indices
  gsea.hit.indices = term_data$RANK.IN.GENE.LIST
  
  # Get ES profile
  gsea.es.profile = term_data$RUNNING.ES
  
  # Set enrichment score range
  if (missing(enrichment.score.range)) {
    enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
    enrichment.score.range[1] = enrichment.score.range[1]-(enrichment.score.range[2]-enrichment.score.range[1])*0.1
    enrichment.score.range[2] = enrichment.score.range[2]+(enrichment.score.range[2]-enrichment.score.range[1])*0.1
  }
  
  ## save figs
  if(!missing(pdf)) pdf(pdf)
  if(!missing(png)) png(png)
  
  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
#  dev.new(width = 3, height = 3)

  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(2.5, 0.7, 0.15, 2))
#  layout.show(gsea.layout)
  xlim = c(1, length(gsea.rnk$metric)*1.02)

#  x_pos = c(1,rep(gsea.hit.indices,each=2),length(gsea.rnk$metric))
  x_pos = c(1, gsea.hit.indices, length(gsea.rnk$metric))
  y_pos = c(0, gsea.es.profile,0)
  len_y = length(y_pos)
#  slope = -0.05
  if(slope == "yes"){
	  y_pos0 = rep(slope*(x_pos[2:len_y]-x_pos[2:len_y-1]),each=2) * c(0,0.001)
	  x_pos = c(1,rep(gsea.hit.indices,each=2),length(gsea.rnk$metric))
#  y_pos = c(0,gsea.es.profile[1],rep(y_pos[1:len_y],each=2)+y_pos0)
	  y_pos = c(rep(y_pos[2:len_y-1],each=2)+y_pos0)
  }
  ylim = c(min(y_pos),max(y_pos)) + c(-1,1) * (max(y_pos)-min(y_pos))*0.1
#  y_pos[length(y_pos)] = 0
  
  # Create plots
  par(mar = c(0, 5, 4, 2))
  plot(x_pos,y_pos,
        type = "l", col = color, lwd = 5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",cex.lab=cex_value,
       ylim = ylim,
	   xlim = xlim,
       main = list(title, font = 2, cex = 2),
	   cex.axis=1.5,
       panel.first = {
          abline(h = seq(round(enrichment.score.range[1], digits = 1),
                         enrichment.score.range[2], 0.1),
                 col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score), adj = c(0, 0),cex=1.5)
  } else {
    text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score, "\n"), adj = c(1, 1),cex=1.5)
  }
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = xlim)
  abline(v = gsea.hit.indices, lwd = 1.2,col="black")
  
  par(mar = c(0, 5, 0, 2))
  rank.colors = gsea.rnk$metric
  quar_a = c(0.9,0.6,0.4,0.2)
  quar_b = c(0.8,0.6,0.4,0.1)
  quar = c(as.vector(quantile(gsea.rnk$metric[gsea.rnk$metric>0],probs = quar_a)),0,
           as.vector(quantile(gsea.rnk$metric[gsea.rnk$metric<0],probs = quar_b)),min(gsea.rnk$metric))
  for(i in 1:length(quar)){
	rank.colors[rank.colors>=quar[i]] = i-100
  }
#  print(rank.colors)
#  print(gsea.rnk$metric)
#  print(quar)
  rank.colors = rank.colors + 100
  color_num = length(unique(rank.colors))
#  print(rank.colors)
#  print(colorRampPalette(c("blue", "white", "red"))(color_num))
  tryCatch({
    rank.colors <- colorRampPalette(c("red", "white", "blue"))(color_num)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = xlim)
  box()
#  par(new=TRUE)
  abline(v = gsea.hit.indices, lwd = 0.75,col="grey60")
#  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
#       ylab = "", xlim = xlim)
#  abline(v = gsea.hit.indices, lwd = 0.75)
#  text(length(gsea.rnk$metric) / 2, 0.7,
  text(length(gsea.rnk$metric) * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(length(gsea.rnk$metric) * 0.99, 0.7, "Negative", adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(round(gsea.rnk$metric, digits = 2))
  plot(gsea.rnk$metric, type = "n", xaxs = "i",
	     xlab = "Rank in ordered gene list", xlim = xlim,
	     ylim = metric.range, yaxs = "i",
	     ylab = if(gsea.metric == "None") {"Ranking metric"} else {gsea.metric},
		 cex.lab=cex_value,
		 cex.axis=1.5)

  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
       xlab = "", xlim = xlim,
       ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
       cex.axis=1.5,
       ylab = "", space = 0, add = TRUE)
  box()
  zero_pos = length(gsea.rnk$metric[gsea.rnk$metric>0])
  abline(v = zero_pos, lwd = 0.75,col="grey60",lty=2)
  text(zero_pos*0.7, metric.range[2]*0.05, paste("Zero cross at",zero_pos), adj = c(0,0),offset=0.5,cex=2)

  
  # Reset to default
  par(def.par)
  if(!missing(png) || !missing(pdf)){
	dev.off()
	}
}

### main
replotGSEA(term_info, all_info, rnk_info, name,"tt",pdf=outpdf)

