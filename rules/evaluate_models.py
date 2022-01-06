all_methods = ["gopeaks", "macs2", "seacr-relaxed", "seacr-stringent"]

rule evaluate_macs2:
	input:
		"data/macs2/{sample}_peaks.xls"
	output:
		"data/evaluate_models/macs2_{sample}.txt"
	conda:
		"../envs/eval.yml"
	params:
		standard=get_standard
	threads: 4
	resources:
		mem_mb = 2000
	shell:
		"bash src/custom/evals/evaluate_macs2.sh -i {wildcards.sample} -s {params.standard} -t {threads} > {output}"
# {wildcards.sample} workaround will accommodate for "narrowPeak" and "broadPeak" files.

rule evaluate_seacr_relaxed:
	input:
		"data/seacr/{sample}.relaxed.bed"
	output:
		"data/evaluate_models/seacr-relaxed_{sample}.txt"
	conda:
		"../envs/eval.yml"
	params:
		standard=get_standard
	threads: 4
	resources:
		mem_mb = 2000
	shell:
		"bash src/custom/evals/evaluate_seacr.sh -i {input} -s {params.standard} -t {threads} > {output}"

rule evaluate_seacr_stringent:
	input:
		"data/seacr/{sample}.stringent.bed"
	output:
		"data/evaluate_models/seacr-stringent_{sample}.txt"
	conda:
		"../envs/eval.yml"
	params:
		standard=get_standard
	threads: 4
	resources:
		mem_mb = 2000
	shell:
		"bash src/custom/evals/evaluate_seacr.sh -i {input} -s {params.standard} -t {threads} > {output}"

rule evaluate_gopeaks:
	input:
		"data/peak_counts/gopeaks_{sample}.bed"
	output:
		"data/evaluate_models/gopeaks_{sample}.txt"
	conda:
		"../envs/eval.yml"
	params:
		standard=get_standard
	threads: 4
	resources:
		mem_mb = 2000
	shell:
		"bash src/custom/evals/evaluate_gopeaks.sh -i {input} -s {params.standard} -t {threads} > {output}"

rule plot_roc:
	input:
		expand("data/evaluate_models/{method}_{sample}.txt",
			method = all_methods,
			sample = list(config["STANDARDS"].keys()))
	output:
		directory("data/figures/roc"),
		directory("data/figures/pr")
	conda:
		"../envs/plot.yml"
	script:
		"../src/custom/evals/plot_evals.R"
