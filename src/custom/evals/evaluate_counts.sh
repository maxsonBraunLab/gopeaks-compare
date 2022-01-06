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

calculate_rate () {
	awk -v term_1=$1 -v term_2=$2 'BEGIN {print term_1 / (term_1 + term_2)}' | cut -b1-5
}

calculate_f1_score () {
	awk -v term_1=$1 -v term_2=$2 'BEGIN {print 2 * (term_1 * term_2) / (term_1 + term_2)}' | cut -b1-5
}

threshold() {
	counts=$1
	predicted_truth=$(sort -n -k4 $sample | awk -v s=$1 '$4 >= s')
	predicted_false=$(sort -n -k4 $sample | awk -v s=$1 '$4 < s')
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

# identify method,condition,replicate,mark
method=$(basename $sample | cut -d_ -f1)
condition=$(basename $sample | cut -d_ -f2)
replicate=$(basename $sample | cut -d_ -f3)
mark=$(basename $sample | cut -d_ -f4 | sed 's/.bed//')

echo -e "method\tcondition\treplicate\tmark\tcounts\tTP\tFP\tFN\tTN\tprecision\trecall\tfpr\tf1" # header

# file prep -----------------------------------------------------------------------------

# identify all counts within a sample.
sample_signals=$(cut -f4 $sample | sort -n | uniq)

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



parallel --keep-order --trim lr -j $threads threshold ::: $sample_signals # process peaks in parallel!