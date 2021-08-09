#!/usr/bin/bash

# reset current sample_signal results
if [ -d "data/sample_signal" ]; then
	rm -r "data/sample_signal"
	mkdir "data/sample_signal"
else
	mkdir "data/sample_signal"
fi

for bed in $(find data/gopeaks -name "*.bed" | grep -v "IgG" | sort); do

	file=$(basename $bed | cut -d. -f1)
	method="gopeaks"
	condition=$(basename $file | cut -d_ -f1)
	replicate=$(basename $file | cut -d_ -f2)
	mark=$(basename $file | cut -d_ -f3)

	for bam in $(find data/ban -name "*${condition}_${replicate}_${mark}*.bam" | grep -v "IgG"); do

		outfile="data/sample_signal/${method}_${condition}_${replicate}_${mark}.bed"

		cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
		echo $cmd
		eval $cmd
	done
done
wait

for bed in $(find data/macs2 -name "*.narrowPeak" | grep -v "IgG" | sort); do

	file=$(basename $bed | cut -d. -f1)
	method="macs2"
	condition=$(basename $file | cut -d_ -f1)
	replicate=$(basename $file | cut -d_ -f2)
	mark=$(basename $file | cut -d_ -f3)
	echo $bed

	for bam in $(find data/ban -name "*${condition}_${replicate}_${mark}*.bam" | grep -v "IgG"); do

		echo $bam
		outfile="data/sample_signal/${method}_${condition}_${replicate}_${mark}.bed"

		cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
		echo $cmd
		eval $cmd
	done
done
wait

for bed in $(find data/seacr -name "*.relaxed.bed" | grep -v "IgG" | sort); do

	file=$(basename $bed | cut -d. -f1)
	method="seacr-relaxed"
	condition=$(basename $file | cut -d_ -f1)
	replicate=$(basename $file | cut -d_ -f2)
	mark=$(basename $file | cut -d_ -f3)
	echo $bed

	for bam in $(find data/ban -name "*${condition}_${replicate}_${mark}*.bam" | grep -v "IgG"); do

		echo $bam
		outfile="data/sample_signal/${method}_${condition}_${replicate}_${mark}.bed"

		cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
		echo $cmd
		eval $cmd
	done
done
wait

for bed in $(find data/seacr -name "*.stringent.bed" | grep -v "IgG" | sort); do

	file=$(basename $bed | cut -d. -f1)
	method="seacr-stringent"
	condition=$(basename $file | cut -d_ -f1)
	replicate=$(basename $file | cut -d_ -f2)
	mark=$(basename $file | cut -d_ -f3)
	echo $bed

	for bam in $(find data/ban -name "*${condition}_${replicate}_${mark}*.bam" | grep -v "IgG"); do

		echo $bam
		outfile="data/sample_signal/${method}_${condition}_${replicate}_${mark}.bed"

		cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
		echo $cmd
		eval $cmd
	done
done
wait
