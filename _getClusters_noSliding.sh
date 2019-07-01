#!/bin/bash

## find enriched gene clusters on chromosomes
## eg. fist make 50-gene blocks, and test in each block functional enrichment

#### need the following input files

## grep gene ../Sm_v7.3_Basic-noGaps.gff | cut -d \; -f1 | awk '{print $9, $1, $4}'|sed 's/ID=//'| sort | grep -v 'H\|U\|SM_V7_100\|SM_V7_300\|W0' | sort -k2,2 -k3,3n> gene-chr-start.txt

#1 gene-chr-start.txt	gene chromosome start_coord (only for chromosomes)
# !!!sorted by chr and coord
# Smp_329140 SM_V7_1 88327
# Smp_315690 SM_V7_1 103403
# Smp_317470 SM_V7_1 256087

#2 gene-func.txt	gene and domain ids separated by ,
# Smp_000020	PF07555
# Smp_000040	PF13374,PF13424
# Smp_000050	PF00520

#3 func-names.txt	domain id and name
# PF00001	7tm_1:7_transmembrane_receptor_(rhodopsin_family)
# PF00002	7tm_2:7_transmembrane_receptor_(Secretin_family)

#4 chr-length.txt	chromosome length
# SM_V7_1 88881357
# SM_V7_2 48130368

#### produce files to use ####
cat gene-func.txt | sed 's/,/ /g'|awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}'| awk '{print $2}'| sort | uniq -c | awk '{print $2, $1}'| sort > func-reference.txt
# PF00001 91
# PF00002 7
# PF00003 3

sort gene-chr-start.txt | join -a1 - gene-func.txt > gene-chr-func.txt

#### total number of genes to test ####
BLOCK=50
TOTAL="$(wc -l gene-chr-start.txt | awk '{print $1}')"

#### split chromosome and generate Fisher Table
# gene id from each block
rm -rf splitChr/
mkdir splitChr
# change according to the chromosome names
for i in {1..7} ZW; do Rscript 1-splitChr.R $i; done
# check the number of genes in the last block of each chromosome

# make Fisher table
cd splitChr/
for i in SM*.txt; do cat $i | sort | join - ../gene-func.txt | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2}'| sort | uniq -c | awk '{print $2 " " $1}'|join - ../func-reference.txt | join ../func-names.txt - | awk '{print $1 ":" $2 " " $3 " " $4 " " "'$BLOCK'"-$3 " " "'$TOTAL'"-$4}' > $i-iprsum.txt; done

# print file name and columns
for i in SM*iprsum.txt; do awk '{print FILENAME $0}' $i | sed 's/\.txt-iprsum\.txt/ /' > $i-show.txt; done

# make big Fisher Test table; modify the last block if necessary
cat *show.txt > ../FisherTable-"$BLOCK"gene-blocks.txt
rm -rf *iprsum*.txt

#### Fisher test
cd ../
Rscript 2-Fishers_Exact_Test.R

rm -f FisherTable-"$BLOCK"gene-blocks.txt

#### Plot chromosomes and clusters
cat FisherResults.txt | awk '$10<0.05 && $3>=3' > fisher_enriched.txt

# plot start and end coordinate of clusters
cat fisher_enriched.txt | sed 's/:/ /'| awk '{print $1 " " $2 " " $4}'| awk -v OFS='' '{print "cat splitChr/", $1, ".txt | grep Smp| sort | join - gene-chr-func.txt | grep ", $2, "| awk \047{print $3}\047 | sort -nk1 | awk  \047NR==1;END{print}\047 | tr \047\134n\047 \047 \047 | sed \047s/ $//\047", "|awk \047{print \042", $1, "\042, \042", $2, "\042, ", $3, ", ","$1, $2}\047 > ", $1, "-", $2, ".tab"}' > sigPoints-cmds.txt  #Use octal escape sequences ('\047') or printf ('printf "%c", 39') to print single quote under print: \047 single quote; \042 double quote; \134 backslash

sh -e sigPoints-cmds.txt

cat *.tab | sed 's/-/ /'> plot_func-clusters.txt

rm -rf *.tab sigPoints-cmds.txt

Rscript 4-plotClusters.R
