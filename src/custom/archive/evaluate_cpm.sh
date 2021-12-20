#!/bin/bash

# get CPM signal per peak by peak caller.
mkdir data/signal_cpm
mkdir data/signal_rpkm

# norm="cpm"
norm="rpkm"

# gopeaks
for sample in $(find data/gopeaks -name "*.bed" | grep -v "IgG" | sort); do

	# identify the sample
	method="gopeaks"
	condition=$(basename $sample | cut -d_ -f1)
	replicate=$(basename $sample | cut -d_ -f2)
	mark=$(basename $sample | cut -d_ -f3 | cut -d. -f1)
	samp="${condition}_${replicate}_${mark}"

	# file I/O
	bam="data/ban/${samp}.ban.sorted.markd.bam"
	N=$(samtools view -@ 4 -c $bam)

	# Normalization
	if [ $norm == "cpm" ]; then
		echo "Counting CPM signal for $samp"
		outfile="data/signal_cpm/${method}_${samp}.cpm.bed"
		bedtools intersect -C -a $sample -b $bam | awk -v N=$N -v OFS='\t' '{print $1,$2,$3,1e6*$4/N}' > $outfile
	else
		echo "Counting RPKM signal for $samp"
		outfile="data/signal_rpkm/${method}_${samp}.rpkm.bed"
		bedtools intersect -C -a $sample -b $bam | awk -v N=$N -v OFS='\t' '{print $1,$2,$3,$4/( $3-$2/1000 * N/1e6 )}' > $outfile
	fi

done

# macs2
for sample in $(find data/macs2 -name "*Peak" | grep -v "IgG" | sort); do

	# identify the sample
	method="macs2"
	condition=$(basename $sample | cut -d_ -f1)
	replicate=$(basename $sample | cut -d_ -f2)
	mark=$(basename $sample | cut -d_ -f3)
	samp="${condition}_${replicate}_${mark}"

	# file I/O
	bam="data/ban/${samp}.ban.sorted.markd.bam"
	N=$(samtools view -@ 4 -c $bam)

	# Normalization
	if [ $norm == "cpm" ]; then
		echo "Counting CPM signal for $samp"
		outfile="data/signal_cpm/${method}_${samp}.cpm.bed"
		cut -f1-3 $sample | bedtools intersect -C -a - -b $bam | awk -v N=$N -v OFS='\t' '{print $1,$2,$3,1e6*$4/N}' > $outfile
	else
		echo "Counting RPKM signal for $samp"
		outfile="data/signal_rpkm/${method}_${samp}.rpkm.bed"
		cut -f1-3 $sample | bedtools intersect -C -a - -b $bam | awk -v N=$N -v OFS='\t' '{print $1,$2,$3,$4/( $3-$2/1000 * N/1e6 )}' > $outfile
	fi

done

# seacr
for sample in $(find data/seacr -name "*.bed" | grep -v "IgG" | sort); do

	# identify the sample
	method="seacr-$(basename $sample | cut -d_ -f3 | cut -d. -f2)" # "seacr-relaxed" or "seacr-stringent"
	condition=$(basename $sample | cut -d_ -f1)
	replicate=$(basename $sample | cut -d_ -f2)
	mark=$(basename $sample | cut -d_ -f3 | cut -d. -f1)
	samp="${condition}_${replicate}_${mark}"

	# file I/O
	bam="data/ban/${samp}.ban.sorted.markd.bam"
	N=$(samtools view -@ 4 -c $bam)

	# Normalization
	if [ $norm == "cpm" ]; then
		echo "Counting CPM signal for $samp"
		outfile="data/signal_cpm/${method}_${samp}.cpm.bed"
		cut -f1-3 $sample | bedtools intersect -C -a - -b $bam | awk -v N=$N -v OFS='\t' '{print $1,$2,$3,1e6*$4/N}' > $outfile
	else
		echo "Counting RPKM signal for $samp"
		outfile="data/signal_rpkm/${method}_${samp}.rpkm.bed"
		cut -f1-3 $sample | bedtools intersect -C -a - -b $bam | awk -v N=$N -v OFS='\t' '{print $1,$2,$3,$4/( $3-$2/1000 * N/1e6 )}' > $outfile
	fi

done