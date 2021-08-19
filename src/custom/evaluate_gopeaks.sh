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

file=$(find $sample)

# identify method,condition,replicate,mark
method="gopeaks"
condition=$(basename $sample | cut -d_ -f2)
replicate=$(basename $sample | cut -d_ -f3)
mark=$(basename $sample | cut -d_ -f4 | cut -d. -f1)

echo -e "method\tcondition\treplicate\tmark\tsignal\tTP\tFP\tFN\tTN\tprecision\trecall\tfpr" # header

# file prep -----------------------------------------------------------------------------

# identify all pvals within a sample.
sample_file=$(sort -n -k4 $file)
sample_signals=$(echo "$sample_file" | cut -f4 | uniq)

for signal in $sample_signals; do

	predicted_truth=$(echo "$sample_file" | awk -v s=$signal '$4 >= s')
	predicted_false=$(echo "$sample_file" | awk -v s=$signal '$4 < s')
	TP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -u -a - -b $standard | wc -l)
	FP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)
	FN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -u -a - -b $standard | wc -l)
	TN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)

	# calculate precision, recall, F1. export results.
	precision=$(calculate_positive_rates $TP $FP)
	recall=$(calculate_positive_rates $TP $FN)
	fpr=$(calculate_positive_rates $FP $TN)

	echo -e "$method\t$condition\t$replicate\t$mark\t$signal\t$TP\t$FP\t$FN\t$TN\t$precision\t$recall\t$fpr"

done

