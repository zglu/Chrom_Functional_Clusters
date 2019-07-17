#!/bin/bash

## find enriched gene clusters on chromosomes
## eg. make 50-gene blocks, and test in each block functional enrichment
## usage: ./getClusters_nonSliding.sh [BLOCKSIZE]


##--*--## INPUT FILES ####
#1 gene-chr-start.txt	gene chromosome start_coord (only for chromosomes)
# !!!sorted by chr and coord
# Smp_329140 SM_V7_1 88327
## grep gene ../Sm_v7.3_Basic-noGaps.gff | cut -d \; -f1 | awk '{print $9, $1, $4}'|sed 's/ID=//'| sort | grep -v 'H\|U\|SM_V7_100\|SM_V7_300\|W0' | sort -k2,2 -k3,3n> gene-chr-start.txt

#2 gene-func.txt	gene and function annotations
# Smp_000040	PF13374,PF13424

#3 func-names.txt	function id and name SHOULD NOT CONTAIN : ' or "
# PF00001	7tm_1

#4 chr-length.txt	chromosome length
# SM_V7_1 88881357

# CHECK ARGUMENTS
if [[ $# -lt 1 ]] ; then
    echo 'Usage: ./getClusters_nonSliding.sh [BLOCKSIZE]'
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


echo "Splitting chromosomes..."
##--*--## SPLIT CHROMOSOMES INTO GENE BLOCKS ####
## total number of genes to test ####
TOTAL="$(wc -l gene-chr-start.txt | awk '{print $1}')"
# gene block size
declare -i BLOCK=$1

# create for each gene block a file with genes inside
mkdir splitChr
for i in $(awk '{print $1}' chr-length.txt); do Rscript 1-splitChr.R $i "$BLOCK"; done

# calculate for each function number of genes and make Fisher test table
cd splitChr/
for i in *; do sort $i | join - ../gene-func.txt | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2}'| sort | uniq -c | awk '{print $2 " " $1}'|join - ../func-reference.txt | join ../func-names.txt - | awk '{print $1 ":" $2 " " $3 " " $4 " " "'$BLOCK'"-$3 " " "'$TOTAL'"-$4}' > $i-sum.txt; done

# print file name and columns
for i in *sum.txt; do awk '{print FILENAME $0}' $i | sed 's/\-sum\.txt/ /' > $i-table.txt; done

## delete empty table files
find . -size  0 -print0 |xargs -0 rm --


echo "Testing enrichment..."
##--*--## FISHER TEST ####
for i in *-table.txt; do Rscript ../2-Fishers_Exact_Test.R $i; done
cat FisherResults*.txt | grep -v anno_block |awk '{print $1,$2,$3,$10}'> ../fisher_enriched_nonsliding.txt
cd ../


echo "Plotting clusters..."
##--*--## PLOT CHROMOSOMES AND CLUSTERS ####
# get coordinates of first and last genes in each cluster
cat fisher_enriched_nonsliding.txt | sed 's/:/ /'| awk -v OFS='' '{print "cat splitChr/", $1, " | sort | join - gene-chr-func.txt | grep ", $2, "| awk \047{print $3}\047 | sort -nk1 | awk  \047NR==1;END{print}\047 | tr \047\134n\047 \047 \047 | sed \047s/ $//\047", "|awk \047{print \042", $1, "\042, \042", $2, "\042, \042", $3, "\042, ", $4, ", ", $5, ", ","$1, $2}\047 > ", $1, "-", $2, ".tab"}' > sigPoints-cmds.txt  #Use octal escape sequences ('\047') or printf ('printf "%c", 39') to print single quote under print: \047 single quote; \042 double quote; \134 backslash
sh -e sigPoints-cmds.txt
cat *.tab | sed 's/-/ /' | sort | join - chr-length.txt > plot_func-clusters_nonsliding.txt

# plot with gene counts in each cluster
#Rscript 3-cluster_geneCounts.R plot_func-clusters_sliding.txt
# plot with gene coordinates of each cluster
Rscript 4-cluster_geneCoord.R plot_func-clusters_nonsliding.txt

rm -rf *.tab sigPoints-cmds.txt
rm -rf splitChr/

echo "Done!"