#!/bin/bash

# reset current intervene results
if [ -d "data/intervene" ]; then
	rm -r "data/intervene"
	mkdir "data/intervene"
else
	mkdir "data/intervene"
fi

while getopts "i:" op
do
	case "$op" in
		i)  i="$OPTARG";;
		\?) exit 1;;
	esac
done

# check if input files exists.
if [ -f $i ];
then
	echo "$i file exists"
else
	echo -e "$i does not exist. Exiting program."
	exit
fi

# intervene to find 4-way venn diagrams of peak sets
while IFS=$'\t' read condition mark
do

	# define I/O
	input_gopeaks="data/consensus/gopeaks_${condition}_${mark}.bed"
	input_macs2="data/consensus/macs2_${condition}_${mark}.bed"
	input_seacr_relaxed="data/consensus/seacr-relaxed_${condition}_${mark}.bed"
	input_seacr_stringent="data/consensus/seacr-stringent_${condition}_${mark}.bed"
	condition_mark=$(echo ${condition}_${mark})
	output_venn="data/intervene/$condition_mark"
	log_file="data/logs/intervene_${condition}_${mark}.log"

	echo "intervene | condition,mark: $condition,$mark"

	# define parameters
	title=$(echo "$condition_mark Consensus Peaks")

	# intervene
	cmd="intervene venn -i $input_gopeaks $input_macs2 $input_seacr_relaxed $input_seacr_stringent \
	--title '$title' \
	--names='gopeaks','macs2','seacr-relaxed','seacr-stringent' \
	--save-overlaps \
	--colors='#2C7BB6','#ABD9E9','#FDAE61','#D7191C' \
	--fontsize 12 \
	-o $output_venn"
	echo "$cmd" > $log_file 2>&1
	eval "$cmd" &

done < $i

echo "Taking intersection of consensus peaks..."
wait