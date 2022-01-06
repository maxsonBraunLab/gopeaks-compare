#!/usr/bin/bash

mkdir -p data/scaledHeatmap
mkdir -p data/computeMatrix

cut -f4,5 samplesheet.tsv | tail -n+2 | grep -v "IgG" | sort | uniq | awk -v OFS="\t" '{print $2,$1}' > all_groups.txt

p=8
a=5000
b=5000

while IFS=$'\t' read condition mark
do

	# define I/O
	input_bw=$(find "data/tracks" -name "*${condition}*${mark}.bw" | sort | tr '\n' ' ')
	input_gopeaks="data/consensus/gopeaks_${condition}_${mark}.bed"
	input_macs2="data/consensus/macs2_${condition}_${mark}.bed"
	input_seacr_relaxed="data/consensus/seacr-relaxed_${condition}_${mark}.bed"
	input_seacr_stringent="data/consensus/seacr-stringent_${condition}_${mark}.bed"
	output_mtx="data/computeMatrix/${condition}_${mark}.gz"
	log_file="data/logs/heatmaps_${condition}_${mark}.log"

	# computeMatrix
	cmd=$(echo "computeMatrix scale-regions -S $input_bw -R $input_gopeaks $input_macs2 $input_seacr_relaxed $input_seacr_stringent -a $a -b $b -p $p --smartLabel -o $output_mtx")
	echo $cmd
	# eval $cmd
	eval "sbatch -c $p -t 06:00:00 --wait --wrap='$cmd'" >> $log_file 2>&1 &
	# echo $cmd > $log_file # echo cmd into log file

done < all_groups.txt

wait

# plot things

for input_mtx in $(ls data/computeMatrix/*.gz); do

	sample=$(basename $input_mtx | cut -d"." -f1)
	output_heatmap="data/scaledHeatmap/${sample}.pdf"

	cmd=$(echo "plotHeatmap -m $input_mtx -o $output_heatmap")
	echo $cmd
	eval $cmd &

done

wait
