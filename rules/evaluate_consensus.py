st = pd.read_table('samplesheet.tsv').set_index('sample',drop=False)
all_methods = ["gopeaks", "macs2", "seacr-relaxed", "seacr-stringent"]
all_groups = get_groups() # df of method, cond, and mark, columns w/out igg.
all_consensus = expand("data/consensus/{method}_{condition}_{mark}.bed", zip,
					method = list(all_groups.method),
					condition = list(all_groups.condition),
					mark = list(all_groups.mark)),

rule gopeaks_consensus_counts:
	input:
		"data/gopeaks/{sample}_peaks.bed",
		"data/seacr/bam/{sample}.name.sorted.bam",
		all_consensus
	output:
		"data/consensus_counts/gopeaks_{sample}.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cond=$(echo {wildcards.sample} | cut -d_ -f1)
		mark=$(echo {wildcards.sample} | cut -d_ -f3)
		consensus=$(find data/consensus -name "gopeaks*${{cond}}*${{mark}}*.bed")
		# count read pair pileup at consensus peaks, rather than {input[0]} aka sample-specific peak set.

		bedtools bamtobed -bedpe -i {input[1]} | \
		awk '$1==$4 && $6-$2 < 1000 {{print $0}}' | \
		cut -f1,2,6 | \
		sort -k1,1 -k2,2n -k3,3n | \
		bedtools intersect -C -a $consensus -b stdin > {output}
		"""

rule macs2_consensus_counts:
	input:
		"data/macs2/{sample}_peaks.bed",
		"data/seacr/bam/{sample}.name.sorted.bam",
		all_consensus
	output:
		"data/consensus_counts/macs2_{sample}.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cond=$(echo {wildcards.sample} | cut -d_ -f1)
		mark=$(echo {wildcards.sample} | cut -d_ -f3)
		consensus=$(find data/consensus -name "macs2*${{cond}}*${{mark}}*.bed")
		# count read pair pileup at consensus peaks, rather than {input[0]} aka sample-specific peak set.

		bedtools bamtobed -bedpe -i {input[1]} | \
		awk '$1==$4 && $6-$2 < 1000 {{print $0}}' | \
		cut -f1,2,6 | \
		sort -k1,1 -k2,2n -k3,3n | \
		bedtools intersect -C -a $consensus -b stdin > {output}
		"""

rule seacr_relaxed_consensus_counts:
	input:
		"data/seacr/{sample}.relaxed.bed",
		"data/seacr/bam/{sample}.name.sorted.bam",
		all_consensus
	output:
		"data/consensus_counts/seacr-relaxed_{sample}.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cond=$(echo {wildcards.sample} | cut -d_ -f1)
		mark=$(echo {wildcards.sample} | cut -d_ -f3)
		consensus=$(find data/consensus -name "seacr-relaxed*${{cond}}*${{mark}}*.bed")
		# count read pair pileup at consensus peaks, rather than {input[0]} aka sample-specific peak set.

		bedtools bamtobed -bedpe -i {input[1]} | \
		awk '$1==$4 && $6-$2 < 1000 {{print $0}}' | \
		cut -f1,2,6 | \
		sort -k1,1 -k2,2n -k3,3n | \
		bedtools intersect -C -a $consensus -b stdin > {output}
		"""

rule seacr_stringent_consensus_counts:
	input:
		"data/seacr/{sample}.stringent.bed",
		"data/seacr/bam/{sample}.name.sorted.bam",
		all_consensus
	output:
		"data/consensus_counts/seacr-stringent_{sample}.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cond=$(echo {wildcards.sample} | cut -d_ -f1)
		mark=$(echo {wildcards.sample} | cut -d_ -f3)
		consensus=$(find data/consensus -name "seacr-stringent*${{cond}}*${{mark}}*.bed")
		# count read pair pileup at consensus peaks, rather than {input[0]} aka sample-specific peak set.

		bedtools bamtobed -bedpe -i {input[1]} | \
		awk '$1==$4 && $6-$2 < 1000 {{print $0}}' | \
		cut -f1,2,6 | \
		sort -k1,1 -k2,2n -k3,3n | \
		bedtools intersect -C -a $consensus -b stdin > {output}
		"""

rule evaluate_consensus_counts:
	input:
		"data/consensus_counts/{method}_{sample}.bed"
	output:
		"data/evaluate_consensus_counts/{method}_{sample}.txt"
	params:
		standard = get_standard
	wildcard_constraints:
		method = "gopeaks|macs2|seacr-relaxed|seacr-stringent"
	conda:
		"../envs/eval.yml"
	threads: 8
	shell:
		"bash src/custom/evals/evaluate_counts.sh -i {input} -s {params.standard} -t {threads} > {output}"

rule plot_consensus_counts:
	input:
		expand("data/evaluate_consensus_counts/{method}_{sample}.txt",
			method = all_methods,
			sample = list(config["STANDARDS"].keys()))
	output:
		directory("data/figures-evaluate-consensus-counts/roc_counts"),
		directory("data/figures-evaluate-consensus-counts/pr_counts")
	conda:
		"../envs/plot.yml"
	script:
		"../src/custom/evals/plot_consensus_counts.R"

rule plot_consensus_peak_plots:
	input:
		expand("data/evaluate_consensus_counts/{method}_{sample}.txt",
			method = all_methods,
			sample = list(config["STANDARDS"].keys()))
	output:
		directory("data/figures-evaluate-consensus-counts/consensus_peak_plots")
	conda:
		"../envs/plot.yml"
	script:
		"../src/custom/consensus_peak_plots.R"

rule consensus_standard_intersection:
	input:
		expand("data/consensus/{method}_{condition}_{mark}.bed", zip,
			method = list(all_groups.method),
			condition = list(all_groups.condition),
			mark = list(all_groups.mark))
	output:
		"data/consensus/consensus-standard-intersection/intersections-summary.txt"
	conda:
		"../envs/bedtools.yml"
	shell:
		"bash src/custom/consensus_standard_intersection.sh >> {output}"

rule consensus_standard_intersection_plot:
	input:
		"data/consensus/consensus-standard-intersection/intersections-summary.txt"
	output:
		"data/figures-evaluate-consensus-counts/consensus-peak-standard-intersections.pdf"
	conda:
		"../envs/plot.yml"
	script:
		"../src/custom/consensus_standard_intersection.R"