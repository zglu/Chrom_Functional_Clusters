#!/bin/bash

## find enriched gene clusters on chromosomes
## making 50-gene sliding windows and test in each block functional enrichment

#### need the following input files

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

##--*--##
#### produce files to use ####
# total number of genes for each function (as a reference for the test)
cat gene-func.txt | sed 's/,/ /g'|awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}'|sort -u| awk '{print $2}'| sort | uniq -c | awk '{print $2, $1}'| sort > func-reference.txt
# PF00001 91
# PF00002 7
# PF00003 3
sort gene-chr-start.txt | join -a1 - gene-func.txt > gene-chr-func.txt


##--*--##
#### make sliding gene blocks and do Fisher test ####
# split gene chrom file into each chromosome
for i in $(awk '{print $2}' gene-chr-start.txt | sort -u); do grep $i gene-chr-start.txt > $i.txt; done
# construct the fisher table for each sliding block, and test significance
rm -rf slidingChr/
mkdir slidingChr/
for i in $(cut -f1 chr-length.txt); do ./sliding_window.sh $i.txt; done

## delete empty table files
find . -size  0 -print0 |xargs -0 rm --

## Fisher test
for i in *.table; do Rscript 2-Fishers_Exact_Test.R $i; done

# select clusters with most genes
cat FisherResults-*.table | grep -v anno_block | sort -k2,2 -k1,1 > fisher_enriched_raw.txt
#Block   IPR     anno_block      anno_all        nonanno_block   nonanno_all     p_value Bonferroni      Holm    FDR
#SM_V7_1.txt-110	PF00001:7tm_1	5	91	95	9708	0.0028	0.2212	0.21	0.04424
rm -f *.table 

# get clusters with most genes
cat fisher_enriched_raw.txt | sort -k2,2 -k3,3n -k1,1| sed 's/\./ /' | awk '{print $3 "-" $1, $1 "." $2 "#" $4 "#" $11}' | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | awk '{print $1, $NF}' | sed 's/-HMN/ HMN/; s/#/ /g' | awk '{print $3, $1, $4, $5}'| sort > fisher_enriched_sliding.txt # needs manual-check for multiple clusters on the same chromosome and duplicated clusters

#SM_V7_1.txt-0 PF00169:PH 3 3e-04
#SM_V7_1.txt-110 PF00001:7tm_1	5 0.011

##--*--##
#### Plot chromosomes and clusters ####
# plot start and end coordinate of clusters
cat fisher_enriched_sliding.txt | sed 's/:/ /'| awk -v OFS='' '{print "cat slidingChr/", $1, " | sort | join - gene-chr-func.txt | grep ", $2, "| awk \047{print $3}\047 | sort -nk1 | awk  \047NR==1;END{print}\047 | tr \047\134n\047 \047 \047 | sed \047s/ $//\047", "|awk \047{print \042", $1, "\042, \042", $2, "\042, \042", $3, "\042, ", $4, ", ", $5, ", ","$1, $2}\047 > ", $1, "-", $2, ".tab"}' > sigPoints-cmds.txt  #Use octal escape sequences ('\047') or printf ('printf "%c", 39') to print single quote under print: \047 single quote; \042 double quote; \134 backslash

sh -e sigPoints-cmds.txt

# func-names without :
# cat *.tab | sed 's/-/ /' | awk '{print $3, $0}'| sort | join func-names.txt - | awk '{print $3, $4, $2, $6, $7, $8}'| sort |sed 's/\.txt//'> plot_func-clusters.txt
cat *.tab | sed 's/-/ /; s/\.txt//'|sort > plot_func-clusters_sliding.txt # file used to plot
# chr block func-id func-name genes fdr start end
#SM_V7_1 110 PF00001 7tm_1 5 0.001 42207659 45304828

rm -rf *.tab sigPoints-cmds.txt fisher_enriched_sliding.txt

Rscript 3-plotClusters.R plot_func-clusters_sliding.txt
