#!/bin/bash

declare -i WINSIZE=100
declare -i SLIDING=10
NGENES="$(wc -l $1 | awk '{print $1}')"
NSLIDES=$(((NGENES - WINSIZE) / SLIDING))
TOTAL="$(wc -l gene-chr-start.txt | awk '{print $1}')"

for ((i=0; i<=NSLIDES;i++)); do
	START=$((1+$i*SLIDING))
	END=$((START+(WINSIZE-1))) # not considering the last 4 or less genes
#	echo $TMPV $START $END
	awk -v start="$START" -v end="$END" 'NR>=start && NR<=end' $1 > slidingChr/$1-$i
	sort slidingChr/$1-$i | join - gene-func.txt | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2}'| sort | uniq -c | awk '{print $2 " " $1}'|join - func-reference.txt | join func-names.txt - | awk '{print $1 ":" $2 " " $3 " " $4 " " "'$WINSIZE'"-$3 " " "'$TOTAL'"-$4}' > $1-$i-sum.txt
	awk '{print FILENAME $0}' $1-$i-sum.txt | sed 's/-sum\.txt/ /' > $1-$i-table.txt
	Rscript 2-Fishers_Exact_Test.R $1-$i-table.txt
done

rm -f *sum.txt SM*-table.txt
