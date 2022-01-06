#!/bin/bash

while getopts "i:s:t:" op
do
	case "$op" in
		i)  sample="$OPTARG";;
		s)  standard="$OPTARG";;
		t)  threads="$OPTARG";;
		\?) exit 1;;
	esac
done

# variables and functions ---------------------------------------------------------------

calculate_rate() {
	awk -v true_term=$1 -v false_term=$2 'BEGIN {print true_term / (true_term + false_term)}' | cut -b1-5
}

calculate_f1_score() {
	awk -v true_term=$1 -v false_term=$2 'BEGIN {print 2 * (true_term * false_term) / (true_term + false_term)}' | cut -b1-5
}

threshold() {
	counts=$1
	predicted_truth=$(cat $sample_file_handle | awk -v s=$1 '$11 <= s')
	predicted_false=$(cat $sample_file_handle | awk -v s=$1 '$11 > s')
	TP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -u -a - -b $standard | wc -l)
	FP=$(echo "$predicted_truth" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)
	FN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -u -a - -b $standard | wc -l)
	TN=$(echo "$predicted_false" | cut -f1-3 | bedtools intersect -v -a - -b $standard | wc -l)

	# calculate precision, recall, F1. export results.
	precision=$(calculate_rate $TP $FP)
	recall=$(calculate_rate $TP $FN)
	fpr=$(calculate_rate $FP $TN)
	f1=$(calculate_f1_score $precision $recall)

	echo -e "$method\t$condition\t$replicate\t$mark\t$counts\t$TP\t$FP\t$FN\t$TN\t$precision\t$recall\t$fpr\t$f1"
}

# file I/O ------------------------------------------------------------------------------

file=$(find data/macs2 -name "${sample}*Peak") # this workaround accommodates for narrowPeak and broadPeak inputs

# identify method,condition,replicate,mark
method="macs2"
condition=$(basename $file | cut -d_ -f1)
replicate=$(basename $file | cut -d_ -f2)
mark=$(basename $file | cut -d_ -f3)

echo -e "method\tcondition\treplicate\tmark\tsignal\tTP\tFP\tFN\tTN\tprecision\trecall\tfpr\tf1" # header

# file prep -----------------------------------------------------------------------------

# create a new column with pval in decimal form. Sort by that new column.
# identify all pvals within a sample.
sample_file=$(awk -v OFS='\t' '{print $0,10**-$9}' $file | grep -vi "inf" | sort -rg -k11)
sample_pvals=$(echo "$sample_file" | cut -f11 | uniq)
sample_file_handle="data/evaluate_models/${method}_${condition}_${replicate}_${mark}.tmp.txt" # create tmp file to threshold from
echo "$sample_file" > $sample_file_handle

# export vars and fx to child processes
export -f calculate_rate
export -f threshold
export -f calculate_f1_score
if [ $method != "" ]; then export method; else echo "method $method is undefined" && exit 1; fi
if [ $condition != "" ]; then export condition; else echo "condition $condition is undefined" && exit 1; fi
if [ $replicate != "" ]; then export replicate; else echo "replicate $replicate is undefined" && exit 1; fi
if [ $mark != "" ]; then export mark; else echo "mark $mark is undefined" && exit 1; fi
if [ $sample != "" ]; then export sample; else echo "sample $sample is undefined" && exit 1; fi
if [ $standard != "" ]; then export standard; else echo "standard $standard is undefined" && exit 1; fi
export sample_file_handle

parallel --keep-order --trim lr -j $threads threshold ::: $sample_pvals # process peaks in parallel!

rm $sample_file_handle