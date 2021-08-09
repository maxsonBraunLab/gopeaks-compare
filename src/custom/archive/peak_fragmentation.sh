#!/usr/bin/bash 

while getopts "g:d:o:" op
do
	case "$op" in
		g)  g="$OPTARG";;
		d)  d="$OPTARG";;
		o)  o="$OPTARG";;
		\?) exit 1;;
	esac
done

input_gtf=$g
downstream=$d
output_prefix=$o

# look only at promoters

# reset current promoter results
if [ -d "src/promoters" ]; then
	rm -r "src/promoters"
	mkdir "src/promoters"
else
	mkdir "src/promoters"
fi

if [ -d "data/fragments" ]; then
	rm -r "data/fragments"
	mkdir "data/fragments"
else
	mkdir "data/fragments"
fi

# convert gtf to bed format
output_bed=$(echo "src/${output_prefix}.bed")
zcat $input_gtf | awk '$3=="gene"' | grep "protein_coding" | awk '{print $1,$4,$5,$10"_"$14,0,$7}' | sed 's/"//g' | sed 's/;//g' | sed 's/^/chr/' | tr ' ' '\t' | sort -k1,1 -k2,2n > $output_bed

declare -a upstream_bp=(1000 2000 3000 4000 5000)

echo "Creating promoter buckets. One file per upstream bp."
# adjust promoter regions of genes upstream by specific intervals.
# should range from 1 kb into gene, to ${upstream_bp} amount of bp.
for bp in "${upstream_bp[@]}"; do
	outfile="src/promoters/${output_prefix}.promoters.upstream.${bp}.bed"
	awk -v OFS='\t' -v upstream=$bp -v downstream=$downstream '{
		if ($6 == "+") {
			print $1,$2-upstream,$2+downstream,$4,$5,$6;
		}
		else if ($6 == "-") {
			print $1,$3-downstream,$3+upstream,$4,$5,$6;
		}
		else {
			print "The following gene does not have a strand: ",$4
		}
	}' $output_bed | awk '$2 > 0' > $outfile # remove any intervals with negative start sites. usually in KI chromosomes.
done

# weirdly enough, this does the job for the most part.
# however, a minority of annotation does not match hg38 tracks in IGV. Clear example is CD72.
# its interval is wider in this GTF file than the IGV feature track.

# now intersect promoter regions with consensus peaks by method.
# count number of intersections, plot the distribution with R to assess
# how fragmented the peaks are.
echo "Intersecting promoter buckets with consensus peaks."
for file in $(find src/promoters -name "*.bed" | sort); do
	upstream_bp=$(basename $file | cut -d. -f4)
	echo $file
	for consensus in $(find data/consensus -name "*.bed" | sort); do
		method=$(basename $consensus | cut -d. -f1 | cut -d_ -f1)
		condition=$(basename $consensus | cut -d. -f1 | cut -d_ -f2)
		mark=$(basename $consensus | cut -d. -f1 | cut -d_ -f3)
		outfile="data/fragments/${method}_${condition}_${mark}.${upstream_bp}.bed"
		echo "$file x $consensus"
		bedtools intersect -c -a $file -b $consensus > $outfile
	done
done

