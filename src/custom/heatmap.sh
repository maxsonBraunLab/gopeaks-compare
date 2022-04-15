#!/usr/bin/bash

mkdir -p data/plotHeatmap
mkdir -p data/computeMatrix

# computeMatrix parameters
p=8
a=5000
b=5000

cut -f4,5 samplesheet.tsv | tail -n+2 | grep -v "IgG" | awk -v OFS="\t" '{print $2,$1}' | sort | uniq | while read condition mark; do

	# define I/O
	input_bw=$(find "data/tracks" -name "*${condition}*${mark}.bw" | sort | tr '\n' ' ')
	input_igg=$(find "data/tracks" -name "*${condition}*IgG.bw" | sort | tr '\n' ' ')
	input_gopeaks="data/consensus/gopeaks_${condition}_${mark}.bed"
	input_macs2="data/consensus/macs2_${condition}_${mark}.bed"
	input_seacr_relaxed="data/consensus/seacr-relaxed_${condition}_${mark}.bed"
	input_seacr_stringent="data/consensus/seacr-stringent_${condition}_${mark}.bed"
	output_mtx="data/computeMatrix/${condition}_${mark}.gz"
	log_file="data/logs/heatmaps_${condition}_${mark}.log"

	# computeMatrix
	cmd="computeMatrix reference-point --referencePoint center -S $input_bw $input_igg -R $input_gopeaks $input_macs2 $input_seacr_relaxed $input_seacr_stringent -a $a -b $b -p $p --smartLabel -o $output_mtx"
	echo $cmd
	sbatch -J "${condition}_${mark}" -c $p --mem 8G -t 01:00:00 --wait --wrap="computeMatrix reference-point --referencePoint center -S $input_bw $input_igg -R $input_gopeaks $input_macs2 $input_seacr_relaxed $input_seacr_stringent -a $a -b $b -p $p --smartLabel -o $output_mtx" &
	echo $cmd > $log_file # echo cmd into log file

done

sleep 1h # somehow wait doesn't work

echo "computeMatrix finished"

# plot results ------------------------------------------------------------------------------------

for input_mtx in $(find data/computeMatrix/ -name "*.gz" | sort); do

	sample=$(basename $input_mtx | cut -d. -f1)
	output_heatmap="data/plotHeatmap/${sample}.pdf"

	cmd=$(echo "plotHeatmap -m $input_mtx -o $output_heatmap --colorMap viridis")
	echo "$cmd"
	eval $cmd &
	echo "$cmd" >> $log_file

done

wait
sleep 10m