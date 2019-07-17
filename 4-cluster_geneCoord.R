## plot chromosomes
library(ggplot2)
library(ggrepel)
args = commandArgs(trailingOnly=TRUE)
chrlen<-read.delim("chr-length.txt", sep=" ", header=F, col.names=c("chr", "length"))
segmsall<-read.delim(args[1], sep=" ", header = F, col.names=c("chr", "block", "func","name", "genes","fdr", "start", "end", "chrlen"))
# longest chr
maxchr<-(max(chrlen$length)+2000000)/1000000
segms<-segmsall[which(segmsall$genes>=3),]
#barplot(chrlen$len/1000000, horiz=TRUE, names.arg=c(1:7, "ZW"), col="white", xlab="Length (Mb)", xlim=c(0, 100))
pdf(file=paste0(args[1], "_coord.pdf"),width = 11.7, height = 6.3)
p <- ggplot(data=chrlen, aes(chrlen$chr, chrlen$length/1000000)) + geom_bar(stat="identity", fill="lightgrey", width = 0.4, position = position_stack(reverse = TRUE)) +coord_flip()+labs(title="Clusters",x="Chromosome", y="Coordinate (Mb)")+ylim(0,maxchr)
p<-p + geom_segment(data=segms, aes(x=chr, xend=chr, y=start/1000000, yend=end/1000000, colour=func), size=9) + geom_text_repel(data=segms, mapping=aes(x=chr, y=(start/1000000+end/1000000)/2, label=paste0(name, " (", genes, ")")), size=1.8, angle=30, vjust=-1.5, hjust=-0.1, color="black") + theme_bw() 
p+theme(legend.position = "none")
dev.off()
#geom_segment(aes(colour=pfam)) + theme(legend.position = "top") #different colors based on pfam groups
