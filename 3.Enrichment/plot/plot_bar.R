
args<-commandArgs(T) 

infile = args[1]
outfile = args[2]
max_line = as.numeric(args[3])
name_key = args[4]
xlab = args[5]
title = args[6]
pfc = args[7]
color = args[8]

ylab = paste(name_key,"Percent(%)",sep=" ")


#library(Rmisc)
library(ggplot2)

data = read.table(infile,sep="\t",header=T,stringsAsFactors=FALSE,quote="")
colnames(data)[1] = c("id")
if(pfc != "all"){
	data = data[GetPFCData(pfc,data$class),]
}
if(nrow(data) < 1){
	cat("No any data\n");
	q()
}

if(nrow(data) < max_line){max_line = nrow(data)}
data = data[1:max_line,,drop=F]
if("qvalue" %in% colnames(data)){
	data$qvalue = signif(data$qvalue,digits=2)
}else if("pvalue" %in% colnames(data)){
	data$pvalue = signif(data$pvalue,digits=2)
}else{
	data[,3] = signif(data[,3],digits=2)
}
if(!("per" %in% colnames(data))){
	colnames(data)[4] = "per"
}

## get max len 
id_len = apply(data.frame(data$id),1,nchar)
big260 = id_len > 60
sepline_num = length(id_len[big260])

#head(data)
for(i in 1:nrow(data)){
	data$id[i] = sepline(data$id[i],60)
}
max_len = GetSepMaxLen(data$id)

data$id = factor(x=data$id,levels = rev(data$id))
data$text = paste(data$num," (",round(data[,3],4),")",sep="")
data$text = paste(data$num," (",data[,3],")",sep="")
ymax = data$per*1.2

pq_name = colnames(data)[3]
colnames(data)[3] = "pq"
#color = "#6495ED"
#ymax = ifelse(max(data$per)>50,max(data$per)+50,100)
ymax = max(data$per)*1.4

fig = ggplot(data,aes(id,per,fill=pq))
fig = fig + geom_bar(stat = "identity",width=0.5)+ coord_flip()+ylim(0,ymax)+
      labs(y=ylab,x=xlab,color="",fill=pq_name,title = title)+
	  geom_text(aes(id,label=text,hjust=-.1),size=3,color="black")
if(color == "red_blue") fig = fig + scale_fill_continuous(pq_name, low="red", high = "blue")

width = 5+max_len/10
height = 8+sepline_num*0.15

ggsave(paste(outfile,".png",sep=""),width=8,height=6,fig,limitsize = FALSE)
ggsave(paste(outfile,".pdf",sep=""),width=width,height=height,fig,limitsize = FALSE)
