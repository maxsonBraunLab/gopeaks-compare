#!/usr/bin/bash

# reset current consensus_signal results
if [ -d "data/consensus_signal" ]; then
	rm -r "data/consensus_signal"
	mkdir "data/consensus_signal"
else
	mkdir "data/consensus_signal"
fi

# count reads at the consensus peak level.

for consensus in $(find data/consensus -name "*.bed" | sort); do

	file=$(basename $consensus | cut -d. -f1)
	method=$(basename $file | cut -d_ -f1)
	condition=$(basename $file | cut -d_ -f2)
	mark=$(basename $file | cut -d_ -f3)
	echo "${file} ${method} ${condition} ${mark} ------------------------------"

	for bam in $(find data/ban -name "${condition}*${mark}*.bam" | sort); do

		echo "$bam"
		replicate=$(basename $bam | cut -d. -f1 | cut -d_ -f2)
		outfile="data/consensus_signal/${method}_${condition}_${replicate}_${mark}.bed"

		# cmd="sbatch -t 01:00:00 --wait --wrap='bedtools intersect -C -a $consensus -b $bam > ${outfile}' &"
		cmd="bedtools intersect -C -a $consensus -b $bam > ${outfile} &"
		echo $cmd
		eval $cmd
	done
	wait

done