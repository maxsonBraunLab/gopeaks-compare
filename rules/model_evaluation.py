all_methods = ["gopeaks", "macs2", "seacr-relaxed", "seacr-stringent"]

rule evaluate_macs2:
	input:
		# "data/macs2/{sample}_peaks.xls"
		"data/evaluate_consensus_peaks/peaks/macs2_{sample}.consensus.bed"
	output:
		"data/model_evaluation_consensus_peaks/macs2_{sample}.txt"
	conda:
		"../envs/bedtools.yml"
	params:
		standard=get_standard
	shell:
		"bash src/custom/evaluate_macs2.sh -i {wildcards.sample} -s {params.standard} > {output}"
# {wildcards.sample} workaround will accommodate for "narrowPeak" and "broadPeak" files.

rule evaluate_seacr_relaxed:
	input:
		# "data/seacr/{sample}.relaxed.bed"
		"data/evaluate_consensus_peaks/peaks/seacr-relaxed_{sample}.consensus.bed"
	output:
		"data/model_evaluation_consensus_peaks/seacr-relaxed_{sample}.txt"
	conda:
		"../envs/bedtools.yml"
	params:
		standard=get_standard
	shell:
		"bash src/custom/evaluate_seacr.sh -i {input} -s {params.standard} > {output}"

rule evaluate_seacr_stringent:
	input:
		# "data/seacr/{sample}.stringent.bed"
		"data/evaluate_consensus_peaks/peaks/seacr-stringent_{sample}.consensus.bed"
	output:
		"data/model_evaluation_consensus_peaks/seacr-stringent_{sample}.txt"
	conda:
		"../envs/bedtools.yml"
	params:
		standard=get_standard
	shell:
		"bash src/custom/evaluate_seacr.sh -i {input} -s {params.standard} > {output}"

rule gopeaks_signal:
	input:
		"data/gopeaks/{sample}.bed",
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/gopeaks_signal/{sample}.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"bedtools intersect -C -a {input[0]} -b {input[1]} > {output}"

rule evaluate_gopeaks:
	input:
		# "data/gopeaks_signal/{sample}.bed"
		"data/evaluate_consensus_peaks/peaks/gopeaks_{sample}.consensus.bed"
	output:
		"data/model_evaluation_consensus_peaks/gopeaks_{sample}.txt"
	conda:
		"../envs/bedtools.yml"
	params:
		standard=get_standard
	shell:
		"bash src/custom/evaluate_gopeaks.sh -i {input} -s {params.standard} > {output}"

rule plot_evaluation:
	input:
		expand("data/model_evaluation_consensus_peaks/{method}_{sample}.txt", 
			method = all_methods,
			sample = list(config["STANDARDS"].keys()))
	output:
		directory("data/figures/roc"),
		directory("data/figures/pr")
	conda:
		"../envs/plot.yml"
	script:
		"../src/custom/plot_evaluation.R"

# deprecated --------------------------------------------------------------------------------------

# rule evaluate_normalized_signal:
# 	input:
# 		"data/signal_rpkm/{method}_{sample}.rpkm.bed"
# 	output:
# 		"data/evaluate_rpkm/{method}_{sample}.txt"
# 	conda:
# 		"../envs/bedtools.yml"
# 	params:
# 		standard = lambda wildcards: config["STANDARDS"][wildcards.sample]
# 	wildcard_constraints:
# 		method="gopeaks|macs2|seacr-relaxed|seacr-stringent"
# 	shell:
# 		"bash src/custom/evaluate_signal.sh -i {input} -s {params.standard} > {output}"

# rule evaluate_consensus_peaks:
# 	input:
# 		"data/evaluate_consensus_peaks/peaks/{sample}.bed"
# 	output:
# 		directory("data/model_evaluation_consensus_peaks_consensus")
# 	conda:
# 		"../envs/bedtools.yml"
# 	shell:
# 		"bash src/custom/evaluate_consensus_peaks.sh"
