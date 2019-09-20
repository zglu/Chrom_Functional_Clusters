# Fisher's Exact Test

args = commandArgs(trailingOnly=TRUE)
fishertable <- read.delim(args[1], sep=" ", header=F)
no.tests<-nrow(fishertable)


expotable<-data.frame()

for (i in 1:no.tests) {
  fisherp<-fisher.test(matrix(c(fishertable[i,'V3'],fishertable[i,'V4'],fishertable[i,'V5'],fishertable[i,'V6']),nrow=2), alternative="greater")$p.value  # alternative ="greater" (doesn's make a difference)
  newline<-cbind(fishertable[i,], round(fisherp, digits = 4)) # 
  expotable<-rbind(expotable, newline)
}

colnames(expotable)<-c("block","func", "anno_block", "anno_all", "nonanno_block", "nonanno_all", "p_value")

## p-value adjustment for multiple comparisons ?p.adjust
### control of the family-wise error rate: bonferroni, holm, hochberg, hommel
expotable$Bonferroni <- p.adjust(expotable$p_value, method = "bonferroni") # p-values are multiplied by the number of comparisons (no.test)
expotable$Holm <- p.adjust(expotable$p_value, method = "holm")
### control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses. Less stringent 
expotable$FDR <- p.adjust(expotable$p_value, method = "BH") # Benjamini & Hochberg
#expotable$FDR<-NULL # drop a column
#expotable<-expotable[order(expotable$FDR),]
#expotable<-expotable[ which(expotable$FDR<0.05 & expotable$anno_block>2),]
expotable<-subset(expotable, (as.numeric(expotable[,"FDR"])<0.05 | grepl("<", expotable[,"FDR"])) & (expotable[,"anno_block"]>2))
write.table(expotable, file=paste0("FisherResults-", args[1]), row.names = F, sep="\t", quote = F)
