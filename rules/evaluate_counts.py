all_methods = ["gopeaks", "macs2", "seacr-relaxed", "seacr-stringent"]

rule evaluate_gopeaks_counts:
	input:
		"data/peak_counts/gopeaks_{sample}.bed"
	output:
		"data/evaluate_counts/gopeaks_{sample}.txt"
	params:
		standard=get_standard
	conda:
		"../envs/bedtools.yml"
	shell:
		"bash src/custom/evals/evaluate_counts.sh -i {input} -s {params.standard} > {output}"

rule evaluate_macs2_counts:
	input:
		"data/peak_counts/macs2_{sample}.bed"
	output:
		"data/evaluate_counts/macs2_{sample}.txt"
	params:
		standard=get_standard
	conda:
		"../envs/bedtools.yml"
	shell:
		"bash src/custom/evals/evaluate_counts.sh -i {input} -s {params.standard} > {output}"

rule evaluate_seacr_relaxed_counts:
	input:
		"data/peak_counts/seacr-relaxed_{sample}.bed"
	output:
		"data/evaluate_counts/seacr-relaxed_{sample}.txt"
	params:
		standard=get_standard
	conda:
		"../envs/bedtools.yml"
	shell:
		"bash src/custom/evals/evaluate_counts.sh -i {input} -s {params.standard} > {output}"

rule evaluate_seacr_stringent_counts:
	input:
		"data/peak_counts/seacr-stringent_{sample}.bed"
	output:
		"data/evaluate_counts/seacr-stringent_{sample}.txt"
	params:
		standard=get_standard
	conda:
		"../envs/bedtools.yml"
	shell:
		"bash src/custom/evals/evaluate_counts.sh -i {input} -s {params.standard} > {output}"

rule plot_roc_counts:
	input:
		expand("data/evaluate_counts/{method}_{sample}.txt",
			method = all_methods,
			sample = list(config["STANDARDS"].keys()))
	output:
		directory("data/figures/roc_counts"),
		directory("data/figures/pr_counts")
	conda:
		"../envs/plot.yml"
	script:
		"../src/custom/evals/plot_counts.R"