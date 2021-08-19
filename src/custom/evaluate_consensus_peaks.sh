#!/bin/bash

while getopts "i:s:" op
do
	case "$op" in
		i)  sample="$OPTARG";;
		s)  standard="$OPTARG";;
		\?) exit 1;;
	esac
done

mkdir -p data/evaluate_consensus_peaks/peaks

# subset sample peaks to consensus peaks ----------------------------------------------------------

# gopeaks
method="gopeaks"
for consensus_peak in $(find data/consensus -name "${method}*.bed" | sort); do 

	# get metadata
	condition=$(basename $consensus_peak | cut -d_ -f2)
	mark=$(basename $consensus_peak | cut -d_ -f3 | cut -d. -f1)

	# locate all replicates
	all_replicates=$(find data/$method -name "${condition}*${mark}*.bed")

	# extract consensus peak for each relevant sample
	for sample in $all_replicates; do

		replicate=$(basename $sample | cut -d_ -f2)
		input_bam="data/ban/${condition}_${replicate}_${mark}.ban.sorted.markd.bam"
		outfile="data/evaluate_consensus_peaks/peaks/${method}_${condition}_${replicate}_${mark}.consensus.bed"

		# specific to gopeaks: grab the counts at the peaks
		bedtools intersect -wb -a $consensus_peak -b $sample | cut -f4- | bedtools intersect -C -a - -b $input_bam > $outfile

	done

done

# macs2
method="macs2"
for consensus_peak in $(find data/consensus -name "${method}*.bed" | sort); do 

	# get metadata
	condition=$(basename $consensus_peak | cut -d_ -f2)
	mark=$(basename $consensus_peak | cut -d_ -f3 | cut -d. -f1)

	# locate all replicates
	all_replicates=$(find data/$method -name "${condition}*${mark}*Peak")

	# extract consensus peak for each relevant sample
	for sample in $all_replicates; do

		replicate=$(basename $sample | cut -d_ -f2)
		outfile="data/evaluate_consensus_peaks/peaks/${method}_${condition}_${replicate}_${mark}.consensus.bed"
		bedtools intersect -wb -a $consensus_peak -b $sample | cut -f4- > $outfile

	done

done

# seacr-relaxed
method="seacr-relaxed"
for consensus_peak in $(find data/consensus -name "${method}*.bed" | sort); do 

	# get metadata
	condition=$(basename $consensus_peak | cut -d_ -f2)
	mark=$(basename $consensus_peak | cut -d_ -f3 | cut -d. -f1)

	# locate all replicates
	all_replicates=$(find data/seacr -name "${condition}*${mark}*.relaxed.bed")

	# extract consensus peak for each relevant sample
	for sample in $all_replicates; do

		replicate=$(basename $sample | cut -d_ -f2)
		outfile="data/evaluate_consensus_peaks/peaks/${method}_${condition}_${replicate}_${mark}.consensus.bed"
		bedtools intersect -wb -a $consensus_peak -b $sample | cut -f4- > $outfile

	done

done

method="seacr-stringent"
for consensus_peak in $(find data/consensus -name "${method}*.bed" | sort); do 

	# get metadata
	condition=$(basename $consensus_peak | cut -d_ -f2)
	mark=$(basename $consensus_peak | cut -d_ -f3 | cut -d. -f1)

	# locate all replicates
	all_replicates=$(find data/seacr -name "${condition}*${mark}*.stringent.bed")

	# extract consensus peak for each relevant sample
	for sample in $all_replicates; do

		replicate=$(basename $sample | cut -d_ -f2)
		outfile="data/evaluate_consensus_peaks/peaks/${method}_${condition}_${replicate}_${mark}.consensus.bed"
		bedtools intersect -wb -a $consensus_peak -b $sample | cut -f4- > $outfile

	done

done

# edit the 'evaluate_{methods}' rules to the above outfile directory, and output to a new folder.

