library(ggplot2)
## needs to be sorted by chr and coord
ggo<-read.delim("gene-chr-start.txt", header=F, sep=" ", col.names=c("id", "chr", "start"))

## GO on each chromosome
args = commandArgs(trailingOnly=TRUE)
chr=paste("SM_V7_", args[1], sep="") 

blockSize <- 50 # no. genes in each block

ggosub<-ggo[grep(chr, ggo$chr),]
rownames(ggosub)<-1:nrow(ggosub) # re-assign rownames as line numbers
ggosub.split<-split(ggosub, (as.numeric((rownames(ggosub)))-1) %/% blockSize) # make 100-gene blocks
spleng<-length(ggosub.split)
#lapply(ggosub.split, nrow) # show no. of genes in each split

## split into blocks 
i<-0;
for (i in 0:(spleng-1)) {
  blockstart<-paste("!BlockStart:", chr, ggosub[i*blockSize+1,"start"], sep=" ") # start of each 100-gene block
  testGenes<-ggosub[(i*blockSize+1):((i+1)*blockSize),"id"] # gene ids in that block
  outfile<-paste("splitChr/", chr,"-", i+1, ".txt", sep="")
  ## some codes for fisher's text
  
  write.table(blockstart, file=outfile, quote=F, row.names = F, col.names = F)
  write.table(testGenes, file=outfile, quote=F, row.names = F, col.names = F, append = T)
}  
