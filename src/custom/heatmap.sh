#!/usr/bin/bash

while getopts "i:" op
do
	case "$op" in
		i)  i="$OPTARG";;
		\?) exit 1;;
	esac
done

# reset current computeMatrix results
if [ -d "data/computeMatrix" ]; then
	rm -r "data/computeMatrix"
	mkdir "data/computeMatrix"
else
	mkdir "data/computeMatrix"
fi

# reset current plotHeatmap results
if [ -d "data/plotHeatmap" ]; then
	rm -r "data/plotHeatmap"
	mkdir "data/plotHeatmap"
else
	mkdir "data/plotHeatmap"
fi

# check if input files exists.
if [ -f $i ];
then
	echo "$i file exists"
else
	echo -e "$i does not exist. Exiting program."
	exit
fi

# computeMatrix parameters
p=8
a=5000
b=5000

# put consensus peaks jobs to SLURM, process exclusive peaks in a compute node.

# consensus peaks ---------------------------------------------------------------------------------
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

	echo $input_bw

	echo "computeMatrix | mark,condition: $mark,$condition"

	# computeMatrix
	cmd=$(echo "computeMatrix reference-point --referencePoint center -S $input_bw -R $input_gopeaks $input_macs2 $input_seacr_relaxed $input_seacr_stringent -a $a -b $b -p $p --smartLabel -o $output_mtx")
	echo $cmd
	eval "sbatch -c $p -t 06:00:00 --wait --wrap='$cmd'" >> $log_file 2>&1 &
	echo $cmd > $log_file # echo cmd into log file

done < $i

# exclusive peak ----------------------------------------------------------------------------------

# for exclusive peaks, we use the consensus but method-exclusive peaks to plot at the replicate level
# to show 

# reset current computeMatrix/exclusive results
if [ -d "data/computeMatrix/exclusive" ]; then
	rm -r "data/computeMatrix/exclusive"
	mkdir "data/computeMatrix/exclusive"
else
	mkdir "data/computeMatrix/exclusive"
fi

# reset current plotHeatmap results
if [ -d "data/plotHeatmap/exclusive" ]; then
	rm -r "data/plotHeatmap/exclusive"
	mkdir "data/plotHeatmap/exclusive"
else
	mkdir "data/plotHeatmap/exclusive"
fi

# gopeaks exclusive peaks
method="gopeaks"
for file in $(find data/intervene/ -name "1000_gopeaks.bed" | sort); do

	# only calculate profile plots if a method calls > 100 peaks
	exclusive_peak_count=$(wc -l $file | cut -d" " -f1)
	if [ $exclusive_peak_count -gt 100 ]; then

		# identify sample
		condition_mark=$(echo $file | cut -d/ -f3)
		condition=$(echo $condition_mark | cut -d_ -f1)
		mark=$(echo $condition_mark | cut -d_ -f2)

		# file I/O
		outfile="data/computeMatrix/exclusive/${method}_${condition_mark}.gz"
		input_bw=$(find "data/tracks" -name "*${condition}*${mark}.bw" | sort | tr '\n' ' ') # outputs all replicates

		cmd="computeMatrix reference-point --referencePoint center -S $input_bw -R $file -a $a -b $b -p $p --smartLabel -o $outfile"
		eval "sbatch -c $p -t 00:30:00 --wait --wrap='$cmd'" &
		# eval $cmd & # no sbatch bc exclusive regions are small

	fi
done

# macs2 exclusive peaks
method="macs2"
for file in $(find data/intervene/ -name "0100_macs2.bed" | sort); do

	# only calculate profile plots if a method calls > 100 peaks
	exclusive_peak_count=$(wc -l $file | cut -d" " -f1)
	if [ $exclusive_peak_count -gt 100 ]; then

		# identify sample
		condition_mark=$(echo $file | cut -d/ -f3)
		condition=$(echo $condition_mark | cut -d_ -f1)
		mark=$(echo $condition_mark | cut -d_ -f2)

		# file I/O
		outfile="data/computeMatrix/exclusive/${method}_${condition_mark}.gz"
		input_bw=$(find "data/tracks" -name "*${condition}*${mark}.bw" | sort | tr '\n' ' ') # outputs all replicates

		cmd="computeMatrix reference-point --referencePoint center -S $input_bw -R $file -a $a -b $b -p $p --smartLabel -o $outfile"
		eval "sbatch -c $p -t 00:30:00 --wait --wrap='$cmd'" &

	fi
done

# seacr-relaxed exclusive peaks
method="seacr-relaxed"
for file in $(find data/intervene/ -name "0010_seacr-relaxed.bed" | sort); do

	# only calculate profile plots if a method calls > 100 peaks
	exclusive_peak_count=$(wc -l $file | cut -d" " -f1)
	if [ $exclusive_peak_count -gt 100 ]; then

		# identify sample
		condition_mark=$(echo $file | cut -d/ -f3)
		condition=$(echo $condition_mark | cut -d_ -f1)
		mark=$(echo $condition_mark | cut -d_ -f2)

		# file I/O
		outfile="data/computeMatrix/exclusive/${method}_${condition_mark}.gz"
		input_bw=$(find "data/tracks" -name "*${condition}*${mark}.bw" | sort | tr '\n' ' ') # outputs all replicates

		cmd="computeMatrix reference-point --referencePoint center -S $input_bw -R $file -a $a -b $b -p $p --smartLabel -o $outfile"
		eval "sbatch -c $p -t 00:30:00 --wait --wrap='$cmd'" &

	fi
done
wait # for exclusive peaks

# seacr-stringent does not have any unique peaks

# plot results ------------------------------------------------------------------------------------

wait # for consensus peaks

# plotHeatmap for consensus peaks
echo "Plotting heatmaps for consensus peaks"
for input_mtx in $(ls data/computeMatrix/*.gz); do

	sample=$(basename $input_mtx | cut -d"." -f1)
	output_heatmap="data/plotHeatmap/${sample}.pdf"

	cmd=$(echo "plotHeatmap -m $input_mtx -o $output_heatmap")
	echo "$cmd" >> $log_file
		eval "sbatch -c $p -t 00:30:00 --wait --wrap='$cmd'"

done

# plotProfile for exclusive peaks
echo "Plotting profiles for exclusive peaks"
for input_mtx in $(ls data/computeMatrix/exclusive/*.gz); do

	sample=$(basename $input_mtx | cut -d"." -f1)
	output="data/plotHeatmap/${sample}.pdf"

	cmd=$(echo "plotProfile -m $input_mtx -o $output")
	echo "$cmd" >> $log_file
	eval $cmd &

done

wait

echo "heatmaps made"