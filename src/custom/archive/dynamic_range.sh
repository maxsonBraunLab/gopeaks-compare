#!/usr/bin/bash

if [ -d "data/dynamic_range" ]; then
	rm -r "data/dynamic_range"
	mkdir "data/dynamic_range"
else
	mkdir "data/dynamic_range"
fi

# simply loop through every consensus file, and intersect with all the bam files.
for consensus_file in $(ls data/consensus/*.bed); do

	consensus=$(basename $consensus_file | cut -d. -f1)

	echo $consensus
	method=$(echo $consensus | cut -d_ -f1)
	condition=$(echo $consensus | cut -d_ -f2)
	mark=$(echo $consensus | cut -d_ -f3)
	input_bams=$(find data/ban -name "${condition}*${mark}*.bam" | sort)
	input_bams_entries=$(echo $input_bams | tr ' ' '\n' | wc -l)

	if [[ $input_bams_entries -ne 0 ]]; then

		for bam in $input_bams; do
			echo $bam
			replicate=$(basename $bam | cut -d. -f1 | cut -d_ -f2)
			outfile="data/dynamic_range/${method}_${condition}_${replicate}_${mark}.bed"
			cmd="bedtools intersect -C -a $consensus_file -b $bam > $outfile"
			eval $cmd
		done

	fi
done
