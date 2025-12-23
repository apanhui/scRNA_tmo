
args<-commandArgs(T) 

infile = args[1]
outpfx = args[2]
max_line = as.numeric(args[3])
name_key = args[4]
ylab = args[5]
title = args[6]
rev = as.numeric(args[7])
pfc = args[8]

library(ggplot2)
mat <- read.table(infile, sep = "\t", check.names = 0, header = T, quote="",stringsAsFactors=FALSE)
if(pfc != "all"){
	mat = mat[GetPFCData(pfc,mat$class),]
}
if(nrow(mat) < 1){
	q()
}
if(nrow(mat) < max_line){max_line = nrow(mat)}
mat = mat[1:max_line,,drop=F]
if("qvalue" %in% colnames(mat)){
	if(min(mat$qvalue) < 1e-10) mat$qvalue = signif(mat$qvalue,digits=2)
}else{
	if(min(mat$pvalue) < 1e-10) mat$pvalue = signif(mat$pvalue,digits=2)
}

id_len = apply(data.frame(mat$id),1,nchar)
big260 = id_len > 60
sepline_num = length(id_len[big260])

for(i in 1:nrow(mat)){
    mat$id[i] = sepline(mat$id[i],60)
}
max_len = GetSepMaxLen(mat$id)

pq_name = colnames(mat)[3]
colnames(mat)[3] = "pq"

matx = max(mat$pq)
mati = min(mat$pq)

# rev 
if(rev == 1){
	mat$id = factor(mat$id,levels=mat$id)
}else{
	mat$id = factor(mat$id,levels=rev(mat$id))
}

#head(mat)
p <- ggplot(mat, aes(ratio, id))
if (matx > mati){
	fig = p + geom_point(aes(size = num,colour = pq)) + 
	scale_colour_continuous(pq_name, low="red", high = "blue")+
	guides(colour = guide_colorbar(order=2),size = guide_legend(order=1))
}else{
	fig = p + geom_point(aes(size = num,colour = pq))+
	scale_colour_continuous(pq_name)+
	guides(colour = guide_legend(order=2),size = guide_legend(order=1))
}
fig = fig + xlim(min(mat$ratio)-0.05,max(mat$ratio)+0.05)
fig = fig + scale_size(paste(name_key,"Number",sep="")) + labs(title = title, x = "RichFactor", y = ylab) 
#    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face="bold"))
width = 5+max_len/10
height = 8+sepline_num*0.15

SaveFig(fig=fig,outpfx,width=width,height=height)
