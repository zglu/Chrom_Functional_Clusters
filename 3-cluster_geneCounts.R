library(ggplot2)
library(ggrepel)
args = commandArgs(trailingOnly=TRUE)
segmsall<-read.delim(args[1], sep=" ", header = F, col.names=c("chr", "block", "func", "name", "genes", "fdr", "start", "end", "chrlen"))
#eg SM_V7_1 0 PF00169 PH 3 0.0450167 2212203 4703829 88881357
segms<-segmsall[which(segmsall$genes>=3),]

maxgenes<-max(segms$genes)
pdf(file=paste0(args[1], "_geneCounts.pdf"), width=11.3, height=7.7)
# facet
pp<-ggplot(data=segms, aes(x=start/1000000, y=genes)) + geom_bar(stat="identity", position="dodge", fill="blue", width=0.6) + geom_text_repel(data=segms[which(segms$genes>=5),], mapping=aes(x=start/1000000, y=genes,label=paste0(name, " (", genes, ")")), size=3, color="black")  + geom_rect(mapping=aes(xmin=0, xmax=chrlen/1000000, ymin=0, ymax=maxgenes), fill="white", color="red", alpha=0.01, size=0.5) + facet_grid(rows=vars(chr)) 
pp<-pp+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + labs(x="Coordinate (Mb)", y="Genes") + ylim(0, maxgenes) + theme(legend.position="none") # + theme_classic() / theme_bw()

pp

dev.off()
