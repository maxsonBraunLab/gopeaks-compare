#!/bin/bash

# reset current consensus results
if [ -d "data/consensus" ]; then
	rm -r data/consensus
	mkdir data/consensus
else
	mkdir data/consensus
fi

log="data/logs/consensus.log"
if [ -f $log ]; then
	rm -r $log
else
	touch $log
fi

consensus_peak () {
	cat "$@" | sort -k1,1 -k2,2n | \
	bedtools merge | bedtools intersect -a - -b "$@" -c | \
	awk -v OFS='\t' '$4 >= 2 {print}' | cut -f1-3
}

# define all unique condition and marks
cut -f4,5 samplesheet.tsv | tail -n+2 | grep -v "IgG" | sort | uniq | awk -v OFS='\t' '{print $2,$1}' > data/consensus/all_groups.txt

# loop through condition+mark, take consensus peak split by method.
while IFS=$'\t' read condition mark
do

	# identify all input files
	echo finding consensus peaks by method for $condition $mark
	gopeaks=$(find data/gopeaks -name "${condition}*${mark}.bed")
	macs2=$(find data/macs2 -name "${condition}*${mark}*Peak")
	seacr_relaxed=$(find data/seacr -name "${condition}*${mark}*.relaxed.bed")
	seacr_stringent=$(find data/seacr -name "${condition}*${mark}*.stringent.bed")

	# take consensus peak by method
	consensus_peak $gopeaks > "data/consensus/gopeaks_${condition}_${mark}.bed"
	consensus_peak $macs2 > "data/consensus/macs2_${condition}_${mark}.bed"
	consensus_peak $seacr_relaxed > "data/consensus/seacr-relaxed_${condition}_${mark}.bed"
	consensus_peak $seacr_stringent > "data/consensus/seacr-stringent_${condition}_${mark}.bed"

	# export those cmd to the log file

done < data/consensus/all_groups.txt

# we do not create consensus peak files for each condition bc condition in this analysis
# do not refer to one data set. E.g. don't use K562 peaks to build consensus for Kasumi.