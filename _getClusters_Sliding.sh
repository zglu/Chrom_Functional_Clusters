#!/bin/bash

## find enriched gene clusters on chromosomes
## making 50-gene sliding windows and test in each block functional enrichment
## usage: ./getClusters_Sliding.sh [WINDOWSIZE] [SLIDINGSIZE]
## eg. ./getClusters_Sliding.sh 100 10


##--*--## INPUT FILES ####
#1 gene-chr-start.txt	gene chromosome start_coord (only for chromosomes)
# !!!sorted by chr and coord
# Smp_329140 SM_V7_1 88327

#2 gene-func.txt	gene and function annotations
# Smp_000040	PF13374,PF13424

#3 func-names.txt	function id and name SHOULD NOT CONTAIN : ' or "
# PF00001	7tm_1

#4 chr-length.txt	chromosome length
# SM_V7_1 88881357


# CHECK ARGUMENTS
if [[ $# -lt 2 ]] ; then
    echo 'Usage: ./getClusters_Sliding.sh [WINDOWSIZE] [SLIDINGSIZE]'
    exit 0
fi

# CHECK SPECIAL CHARACTERS IN FUNC-NAMES.TXT
if grep -q "'" func-names.txt; then
	echo "func-names.txt contains single quotes"
elif grep -q '"' func-names.txt; then
	echo "func-names.txt contains double quotes"
elif grep -q ':' func-names.txt; then
	echo "func-names.txt contains :"
	exit 0
fi


echo "Generating files to use..."
##--*--## PRODUCE FILES TO USE ####
# total number of genes for each function (as a reference for the test)
cat gene-func.txt | sed 's/,/ /g'|awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}'|sort -u| awk '{print $2}'| sort | uniq -c | awk '{print $2, $1}'| sort > func-reference.txt
# PF00001 91
sort gene-chr-start.txt | join -a1 - gene-func.txt > gene-chr-func.txt


echo "Making sliding windows on chromosomes..."
##--*--## MAKE SLIDING GENE WINDOWS ####
# split gene chrom file into each chromosome
for i in $(awk '{print $2}' gene-chr-start.txt | sort -u); do grep $i gene-chr-start.txt > $i.txt; done
# create for each sliding window a file with defined number of genes (eg. 100) 
# and calculate for each function the number of annotated genes, to make a fisher test table
# ** DEFINE YOUR WINDOW AND SLIDING SIZE IN sliding_window.sh **
mkdir slidingChr/
for i in $(awk '{print $1}' chr-length.txt); do ./1_Sliding_Window.sh $i.txt $1 $2; done
## delete empty table files
find . -size  0 -print0 |xargs -0 rm --


echo "Testing enrichment for each sliding window..."
##--*--## FISHER TEST FOR EACH WINDOW AND COMPILE RESULTS ####
for i in *.table; do Rscript 2-Fishers_Exact_Test.R $i; done
cat FisherResults-* | grep -v anno_block | sort -k2,2 -k1,1 > fisher_enriched_raw.txt
#Block   func     anno_block      anno_all        nonanno_block   nonanno_all     p_value Bonferroni      Holm    FDR
#SM_V7_1.txt-110	PF00001:7tm_1	5	91	95	9708	0.0028	0.2212	0.21	0.04424
rm -f *.table 


echo "Compiling test results and plotting clusters..."
# select clusters with most genes, need manual check for duplicates and multiple clusters
## get the key string in the chr name, eg. HMN, SM
KEYP=$(awk '(NR==1){print substr($1, 1, 3)}' chr-length.txt)
cat fisher_enriched_raw.txt | sort -k2,2 -k3,3n -k1,1| sed 's/\./ /' | awk '{print $3 "-" $1, $1 "." $2 "#" $4 "#" $11}' | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | awk '{print $1, $NF}' | sed "s/-$KEYP/ $KEYP/" |sed 's/#/ /g' | awk '{print $3, $1, $4, $5}'| sort > fisher_enriched_sliding.txt
# block_file function genes fdr
#SM_V7_1.txt-0 PF00169:PH 3 3e-04


##--*--## PLOT CHROMOSOMES AND CLUSTERS ####
# get coordinates of first and last genes in each cluster
cat fisher_enriched_sliding.txt | sed 's/:/ /'| awk -v OFS='' '{print "cat slidingChr/", $1, " | sort | join - gene-chr-func.txt | grep ", $2, "| awk \047{print $3}\047 | sort -nk1 | awk  \047NR==1;END{print}\047 | tr \047\134n\047 \047 \047 | sed \047s/ $//\047", "|awk \047{print \042", $1, "\042, \042", $2, "\042, \042", $3, "\042, ", $4, ", ", $5, ", ","$1, $2}\047 > ", $1, "-", $2, ".tab"}' > sigPoints-cmds.txt  #Use octal escape sequences ('\047') or printf ('printf "%c", 39') to print single quote under print: \047 single quote; \042 double quote; \134 backslash
sh -e sigPoints-cmds.txt
cat *.tab | sed 's/-/ /; s/\.txt//'|sort | join - chr-length.txt > plot_func-clusters_sliding.txt
# chr block func name genes fdr start end chrlen
# SM_V7_1 110 PF00001 7tm_1 5 0.001 42207659 45304828 88881357

# get genes in each cluster
cat fisher_enriched_sliding.txt | sed 's/:/ /'| awk -v OFS='' '{print "cat slidingChr/", $1, " | sort | join - gene-chr-func.txt | grep ", $2, "| awk \047{print $1}\047 | tr \047\134n\047 \047,\047 > ", $1, "-", $2, ".genes"}' > get_clusterGenes.cmds
sh -e get_clusterGenes.cmds

# plot with gene counts in each cluster
Rscript 3-cluster_geneCounts.R plot_func-clusters_sliding.txt
# plot with gene coordinates of each cluster
#Rscript 4-cluster_geneCoord.R plot_func-clusters_sliding.txt

rm -rf *.tab sigPoints-cmds.txt get_clusterGenes.cmds
#rm -rf slidingChr/

echo "Done!"

# check if multiple clusters on the same chromosome are omitted
echo "Check multiple clusters from the non-sliding approach:"
cat plot_func-clusters_nonsliding.txt | awk '{print $1 "_" $3}'| sort | uniq -c | awk '$1>1'
#cat fisher_enriched_raw.txt | sort -k2,2 -k3,3n -k1,1 | grep [func] | grep [chr] ... > sliding_2.txt
