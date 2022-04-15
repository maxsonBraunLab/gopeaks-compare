#!/usr/bin/bash

# take heatmap of exclusive peaks

a=3000
b=3000
p=8

mkdir -p data/exclusive_heatmaps/computeMatrix

get_file() {
	if [ -f "$1" ]; then
		echo "$1"
	else
		echo ""
	fi
}

awk '{print $5,$4}' samplesheet.tsv | tail -n+2 | grep -vE "Pt176|IgG" | sort | uniq | while read condition mark; do

	echo $condition $mark

	input_bw=$(find "data/tracks" -name "*${condition}*${mark}.bw" | sort | tr '\n' ' ')
	input_igg=$(find "data/tracks" -name "*${condition}*IgG.bw" | sort | tr '\n' ' ')
	input_gopeaks=$(get_file "data/intervene/${condition}_${mark}/sets/1000_gopeaks.bed")
	input_macs2=$(get_file "data/intervene/${condition}_${mark}/sets/0100_macs2.bed")
	input_seacr_relaxed=$(get_file "data/intervene/${condition}_${mark}/sets/0010_seacr-relaxed.bed")
	input_seacr_stringent=$(get_file "data/intervene/${condition}_${mark}/sets/0001_seacr-stringent.bed")
	output_mtx="data/exclusive_heatmaps/computeMatrix/${condition}_${mark}.gz"
	log_file="data/logs/exclusive_heatmaps_${condition}_${mark}.log"

	# computeMatrix
	cmd="computeMatrix reference-point --referencePoint center -S $input_bw $input_igg -R $input_gopeaks $input_macs2 $input_seacr_relaxed $input_seacr_stringent -a $a -b $b -p $p --smartLabel -o $output_mtx"
	echo $cmd
	sbatch -J "${condition}_${mark}" -c $p --mem 8G -t 01:00:00 --wait --wrap="computeMatrix reference-point --referencePoint center -S $input_bw $input_igg -R $input_gopeaks $input_macs2 $input_seacr_relaxed $input_seacr_stringent -a $a -b $b -p $p --smartLabel -o $output_mtx" &
	echo $cmd > $log_file # echo cmd into log file

done

sleep 15m # somehow wait doesn't work

echo "computeMatrix finished"

# plot results ------------------------------------------------------------------------------------

for input_mtx in $(find data/exclusive_heatmaps/computeMatrix/ -name "*.gz" | sort); do

	sample=$(basename $input_mtx | cut -d. -f1)
	output_heatmap="data/exclusive_heatmaps/${sample}.pdf"

	cmd=$(echo "plotHeatmap -m $input_mtx -o $output_heatmap --colorMap viridis")
	echo "$cmd"
	eval $cmd &
	echo "$cmd" >> $log_file

done

wait
sleep 10m