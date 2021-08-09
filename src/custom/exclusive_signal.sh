#!/usr/bin/bash

# reset current exclusive_signal results
if [ -d "data/exclusive_signal" ]; then
	rm -r "data/exclusive_signal"
	mkdir "data/exclusive_signal"
else
	mkdir "data/exclusive_signal"
fi

# exclusive gopeaks signal
echo "gopeaks"
for bed in $(find data/intervene/ -name "1000_gopeaks.bed"); do

	method="gopeaks"
	condition=$(echo $bed | cut -d/ -f3 | cut -d_ -f1)
	mark=$(echo $bed | cut -d/ -f3 | cut -d_ -f2)
	input_bams=$(find data/ban -name "${condition}*${mark}*.bam" | sort)
	input_bams_entries=$(echo $input_bams | tr ' ' '\n' | wc -l)

	if [[ $input_bams_entries -ne 0 ]]; then
		for bam in $input_bams; do
			replicate=$(basename $bam | cut -d. -f1 | cut -d_ -f2)
			outfile="data/exclusive_signal/${method}_${condition}_${replicate}_${mark}.bed"
			cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
			echo "${condition}_${replicate}_${mark}"
			eval $cmd
		done
		wait
	fi
done

echo "seacr_relaxed"
for bed in $(find data/intervene/ -name "0010_seacr-relaxed.bed"); do

	method="seacr-relaxed"
	condition=$(echo $bed | cut -d/ -f3 | cut -d_ -f1)
	mark=$(echo $bed | cut -d/ -f3 | cut -d_ -f2)
	input_bams=$(find data/ban -name "${condition}*${mark}*.bam" | sort)
	input_bams_entries=$(echo $input_bams | tr ' ' '\n' | wc -l)

	if [[ $input_bams_entries -ne 0 ]]; then
		for bam in $input_bams; do
			replicate=$(basename $bam | cut -d. -f1 | cut -d_ -f2)
			outfile="data/exclusive_signal/${method}_${condition}_${replicate}_${mark}.bed"
			cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
			echo "${condition}_${replicate}_${mark}"
			eval $cmd
		done
		wait
	fi
done

echo "seacr_stringent"
for bed in $(find data/intervene/ -name "0001_seacr-stringent.bed"); do

	method="seacr-stringent"
	condition=$(echo $bed | cut -d/ -f3 | cut -d_ -f1)
	mark=$(echo $bed | cut -d/ -f3 | cut -d_ -f2)
	input_bams=$(find data/ban -name "${condition}*${mark}*.bam" | sort)
	input_bams_entries=$(echo $input_bams | tr ' ' '\n' | wc -l)

	if [[ $input_bams_entries -ne 0 ]]; then
		for bam in $input_bams; do
			replicate=$(basename $bam | cut -d. -f1 | cut -d_ -f2)
			outfile="data/exclusive_signal/${method}_${condition}_${replicate}_${mark}.bed"
			cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
			echo "${condition}_${replicate}_${mark}"
			eval $cmd
		done
		wait
	fi
done

echo "macs2"
for bed in $(find data/intervene/ -name "0100_macs2.bed"); do

	method="macs2"
	condition=$(echo $bed | cut -d/ -f3 | cut -d_ -f1)
	mark=$(echo $bed | cut -d/ -f3 | cut -d_ -f2)
	input_bams=$(find data/ban -name "${condition}*${mark}*.bam" | sort)
	input_bams_entries=$(echo $input_bams | tr ' ' '\n' | wc -l)

	if [[ $input_bams_entries -ne 0 ]]; then
		for bam in $input_bams; do
			replicate=$(basename $bam | cut -d. -f1 | cut -d_ -f2)
			outfile="data/exclusive_signal/${method}_${condition}_${replicate}_${mark}.bed"
			cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > $outfile &"
			echo "${condition}_${replicate}_${mark}"
			eval $cmd
		done
		wait
	fi
done

echo "Exclusive signal finished."