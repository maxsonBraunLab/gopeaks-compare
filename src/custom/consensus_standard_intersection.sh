#!/usr/bin/bash

if [ ! -d data/consensus/consensus-standard-intersection ]; then
	mkdir data/consensus/consensus-standard-intersection
fi

out="data/consensus/consensus-standard-intersection/intersections-summary.txt"
log="data/logs/consensus_standard_intersection.log"

echo -e "method\tcondition\tmark\tstandard\tintersections\tpercent" > $out
echo "" > $log

consensus_files=$(find data/consensus -name "*.bed" | grep -v "Pt176" | sort)

for file in $consensus_files; do
	# consensus info
	method=$(basename $file | cut -d_ -f1)
	condition=$(echo $file | cut -d_ -f2)
	mark=$(echo $file | cut -d_ -f3 | cut -d. -f1)

	# define standard from src/config.yaml
	standardFile=$(grep "${condition}.*${mark}" src/config.yml | head -1 | cut -d\" -f2)
	standardID=$(basename $standardFile)
	standardTotal=$(wc -l $standardFile | cut -d" " -f1)

	# i/o and bedtools cmd
	outfile="data/consensus/consensus-standard-intersection/${method}-${condition}-${mark}-${standardID}.bed"
	cmd="bedtools intersect -u -a $file -b $standardFile"

	# intersect consensus peaks and standards
	intersections=$(eval $cmd | wc -l | cut -d" " -f1)
	percent=$(echo "scale=3; $intersections / $standardTotal" | bc)
	echo $cmd >> $log

	echo -e "$method\t$condition\t$mark\t$standardFile\t$intersections\t$percent"
done
