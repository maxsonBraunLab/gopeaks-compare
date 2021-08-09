#!/bin/bash

while getopts "i:s:" op
do
	case "$op" in
		i)  sample="$OPTARG";;
		s)  standard="$OPTARG";;
		\?) exit 1;;
	esac
done

# variables and functions ---------------------------------------------------------------

calculate_positive_rates () {
	awk -v true_term=$1 -v false_term=$2 'BEGIN {print true_term / (true_term + false_term)}' | cut -b1-5
}

calculate_f1_score () {
	awk -v true_term=$1 -v false_term=$2 'BEGIN {print 2 * (true_term * false_term) / (true_term + false_term)}' | cut -b1-5
}

# file I/O ------------------------------------------------------------------------------

file=$(find data/macs2 -name "${sample}*Peak") # this workaround accommodates for narrowPeak and broadPeak inputs

# identify method,condition,replicate,mark
method="macs2"
condition=$(basename $file | cut -d_ -f1)
replicate=$(basename $file | cut -d_ -f2)
mark=$(basename $file | cut -d_ -f3)

echo -e "method\tcondition\treplicate\tmark\tpval\tTP\tFP\tFN\tTN\tprecision\trecall\tfpr" # header

# file prep -----------------------------------------------------------------------------

# create a new column with pval in decimal form. Sort by that new column.
# identify all pvals within a sample.
sample_file=$(awk -v OFS='\t' '{print $0,10**-$9}' $file | sort -rg -k11)
sample_pvals=$(echo "$sample_file" | cut -f11 | uniq)

# pvalue threshold walk -----------------------------------------------------------------

for pval in $sample_pvals; do

	# subset peaks above and below the walking pvalue threshold.
	predicted_truth=$(echo "$sample_file" | awk -v p=$pval '$11 <= p')
	predicted_false=$(echo "$sample_file" | awk -v p=$pval '$11 > p')
	TP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -u -f 0.50 -a - -b $standard | wc -l)
	FP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)
	FN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -u -f 0.50 -a - -b $standard | wc -l)
	TN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)

	# calculate precision, recall, F1. export results.
	precision=$(calculate_positive_rates $TP $FP)
	recall=$(calculate_positive_rates $TP $FN)
	fpr=$(calculate_positive_rates $FP $TN)

	echo -e "$method\t$condition\t$replicate\t$mark\t$pval\t$TP\t$FP\t$FN\t$TN\t$precision\t$recall\t$fpr"

done

