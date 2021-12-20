rule macs2:
	input:
		sample = "data/ban/{sample}.ban.sorted.markd.bam",
		index = "data/ban/{sample}.ban.sorted.markd.bam.bai",
		igg = get_igg
	output:
		"data/macs2/{sample}_peaks.narrowPeak",
		"data/macs2/{sample}_peaks.xls"
	params:
		genome_size=config["GSIZE"],
		peak = macs2_peak,
		igg = macs2_igg
	conda:
		"../envs/macs2.yml"
	log:
		"data/logs/macs2_{sample}.log"
	shell:
		"macs2 callpeak -t {input.sample} {params.igg} "
		"--name {wildcards.sample} --format BAMPE "
		"--gsize {params.genome_size} --outdir data/macs2 "
		"{params.peak} --tempdir data/macs2 > {log} 2>&1"
# {input.igg} makes sure the corresponding IgG is made before macs2 is run.
# Otherwise {params.igg} does not enforce making IgG.

rule get_seacr:
	output:
		"src/SEACR-1.3/SEACR_1.3.sh"
	shell:
		"bash src/wget_seacr.sh"

rule get_seacr_14:
	output:
		"src/SEACR-1.4/SEACR_1.4.sh"
	shell:
		"git clone --branch SEACRv1.4_dev git@github.com:FredHutch/SEACR.git src/SEACR-1.4"

rule name_sort:
	input:
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/seacr/bam/{sample}.name.sorted.bam"
	conda:
		"../envs/sambamba.yml"
	threads: 4
	shell:
		"sambamba sort -n -t {threads} --tmpdir data/seacr/bam -o {output} {input}"

rule bedgraph:
	input:
		rules.name_sort.output
	output:
		"data/seacr/bedgraph/{sample}.bedgraph"
	params:
		chrom_size = config["CSIZES"]
	conda:
		"../envs/bedtools.yml"
	shell:
		"bedtools bamtobed -bedpe -i {input} | "
		"awk '$1==$4 && $6-$2 < 1000 {{print $0}}' | "
		"cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n | "
		"bedtools genomecov -bg -i stdin -g {params.chrom_size} > {output}"

rule seacr_relaxed:
	input:
		seacr = rules.get_seacr.output, # v1.3
		# seacr = rules.get_seacr_14.output, # v1.4
		sample = rules.bedgraph.output,
		igg = get_igg
	output:
		"data/seacr/{sample}.relaxed.bed"
	params:
		igg = seacr_igg,
		norm = seacr_norm
	conda:
		"../envs/seacr.yml"
	log:
		"data/logs/seacr-relaxed_{sample}.log"
	shell:
		"bash {input.seacr} {input.sample} {params.igg} {params.norm} relaxed data/seacr/{wildcards.sample} > {log} 2>&1" # SEACR v1.3
		# "bash {input.seacr} -b {input.sample} -c {params.igg} -n {params.norm} -m relaxed -o data/seacr/{wildcards.sample} > {log} 2>&1" # SEACR v1.4

rule seacr_stringent:
	input:
		seacr = rules.get_seacr.output, # v1.3
		# seacr = rules.get_seacr_14.output,
		sample = rules.bedgraph.output,
		igg = get_igg
	output:
		"data/seacr/{sample}.stringent.bed"
	params:
		igg = seacr_igg,
		norm = seacr_norm
	conda:
		"../envs/seacr.yml"
	log:
		"data/logs/seacr-stringent_{sample}.log"
	shell:
		"bash {input.seacr} {input.sample} {params.igg} {params.norm} stringent data/seacr/{wildcards.sample} > {log} 2>&1" # SEACR v1.3
		# "bash {input.seacr} -b {input.sample} -c {params.igg} -n {params.norm} -m stringent -o data/seacr/{wildcards.sample} > {log} 2>&1" # SEACR v1.4

# FRiP --------------------------------------------------------------------------------------------

rule frip_gopeaks:
	input:
		"data/gopeaks/{sample}.bed",
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/plotEnrichment/gopeaks_{sample}.png",
		"data/plotEnrichment/gopeaks_{sample}.tsv"
	conda:
		"../envs/dtools.yml"
	log: 
		"data/logs/plotEnrichment_gopeaks_{sample}.log"
	shell:
		"plotEnrichment -b {input[1]} --BED {input[0]} --regionLabels 'frip' "
		"--outRawCounts {output[1]} -o {output[0]} > {log} 2>&1"

rule frip_macs2:
	input:
		"data/macs2/{sample}_peaks.narrowPeak",
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/plotEnrichment/macs2_{sample}.png",
		"data/plotEnrichment/macs2_{sample}.tsv"
	conda:
		"../envs/dtools.yml"
	log:
		"data/logs/plotEnrichment_macs2_{sample}.log"
	shell:
		"plotEnrichment -b {input[1]} --BED {input[0]} --regionLabels 'frip' "
		"--outRawCounts {output[1]} -o {output[0]} > {log} 2>&1"

rule frip_seacr_relaxed:
	input:
		rules.seacr_relaxed.output,
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/plotEnrichment/seacr-relaxed_{sample}.png",
		"data/plotEnrichment/seacr-relaxed_{sample}.tsv"
	conda:
		"../envs/dtools.yml"
	log:
		"data/logs/plotEnrichment_seacr-relaxed_{sample}.log"
	shell:
		"plotEnrichment -b {input[1]} --BED {input[0]} --regionLabels 'frip' "
		"--outRawCounts {output[1]} -o {output[0]} > {log} 2>&1"

rule frip_seacr_stringent:
	input:
		rules.seacr_stringent.output,
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/plotEnrichment/seacr-stringent_{sample}.png",
		"data/plotEnrichment/seacr-stringent_{sample}.tsv"
	conda:
		"../envs/dtools.yml"
	log:
		"data/logs/plotEnrichment_seacr-stringent_{sample}.log"
	shell:
		"plotEnrichment -b {input[1]} --BED {input[0]} --regionLabels 'frip' "
		"--outRawCounts {output[1]} -o {output[0]} > {log} 2>&1"