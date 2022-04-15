#!/usr/bin/bash

# one reviewer asked if our peaks intersected with igg peaks.
# i will intersect consensus peaks (mark) with consensus peaks (igg) and only report back if there are overlaps.

mkdir data/consensus_igg

# create consensus igg peaks per method -----------------------------------------------------------
all_methods=( gopeaks macs2 seacr-relaxed seacr-stringent )

# loop over all methods
for method in ${all_methods[@]}; do

	# loop over all conditions
	cut -f5 samplesheet.tsv | tail -n+2 | sort | grep -v "Pt176" | uniq | while read condition; do

		# define I/O
		echo $method $condition
		if [ "$method" == "seacr-relaxed" ]; then
			igg_beds=$(find data/seacr -name "*${condition}*IgG*relaxed.bed" | sort | tr '\n' ' ')
		elif [ "$method" == "seacr-stringent" ]; then
			igg_beds=$(find data/seacr -name "*${condition}*IgG*stringent.bed" | sort | tr '\n' ' ')
		else
			igg_beds=$(find data/$method -name "*${condition}*IgG*.bed" | sort | tr '\n' ' ')
		fi
		outfile="data/consensus_igg/${method}_${condition}_IgG.bed"

		# get consensus peaks for igg
		cat $igg_beds | sort -k1,1 -k2,2n | bedtools merge | \
		bedtools intersect -a stdin -b $igg_beds -c | \
		awk -v OFS='\t' '$4 >= 2 {print}' | cut -f1-3 > $outfile

	done
done

# There are a couple of IgG peaks, but Will doesn't need this type of analysis.

# intersect sample and igg peaks at consensus level -----------------------------------------------
# for consensus in $(find data/consensus -name "*.bed" | grep -v "Pt176" | sort); do

# 	# define I/O
# 	method=$(basename $consensus | cut -d_ -f1)
# 	condition=$(basename $consensus | cut -d_ -f2)
# 	mark=$(basename $consensus | cut -d_ -f3 | cut -d. -f1)
# 	igg_bed="data/consensus_igg/${method}_${condition}_IgG.bed"
# 	outfile="data/consensus_igg/intersection"

# 	echo $method $condition $mark

# 	bedtools intersect -u -a $consensus -b $igg_bed > $

# done