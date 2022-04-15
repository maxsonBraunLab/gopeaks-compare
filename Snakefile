# CutAndQC performs initial QC on CutAndTag projects
import glob
import os 
from pathlib import Path,PurePath,PurePosixPath
from collections import defaultdict
import pandas as pd
from snakemake.utils import validate, min_version
import plotly as plt
import plotly.graph_objects as go
from snakemake.utils import Paramspace

##### set minimum snakemake version #####
min_version("5.32")

include: "src/common.py"
configfile: "src/config.yml"
include: "rules/peaks.py"
include: "rules/evaluate_models.py"
include: "rules/evaluate_counts.py"
include: "rules/evaluate_consensus.py"

st = pd.read_table('samplesheet.tsv').set_index('sample',drop=False)
validate(st, schema="schemas/samples.schema.yml")

samps = get_samples()
reads= get_reads()
marks=get_marks()

sample_noigg = [k for k in samps if config["IGG"] not in k]
marks_noigg = [m for m in marks if config["IGG"] not in m]

fastqScreenDict = {
'database': {
   'hg38': {
     'bowtie2': config["BOWTIE2"]["HG38"][0]},
   'mm10': {
     'bowtie2': config["BOWTIE2"]["MM10"][0]}, 
   'ecoli': {
     'bowtie2': config["BOWTIE2"]["ECOLI"][0]}, 
   'myco': {
     'bowtie2': config["BOWTIE2"]["MYCO"][0]}, 
 },
 'aligner_paths': {'bowtie2': 'bowtie2'}
}

all_methods = ["gopeaks", "macs2", "seacr-relaxed", "seacr-stringent"]
all_groups = get_groups() # df of method, cond, and mark, columns w/out igg.

rule all:
    input:
        # quality control ---------------------------------------------------------------
        expand("data/fastqc/{read}.html", read=reads),
        expand("data/fastq_screen/{read}_screen.txt", read=reads),
        expand(["data/ban/{sample}.ban.sorted.markd.bam",
                "data/tracks/{sample}.bw",
                "data/preseq/lcextrap_{sample}.txt"],
                sample=samps),
        # peak calling ------------------------------------------------------------------
        "src/gopeaks",
        "src/SEACR-1.3/SEACR_1.3.sh",
        "src/SEACR-1.4/SEACR_1.4.sh",
        expand("data/gopeaks/{sample}_peaks.bed", sample = samps),
        expand("data/macs2/{sample}_peaks.xls", sample = samps),
        expand("data/seacr/{sample}.{type}.bed",
            sample = samps,
            type = ["relaxed", "stringent"]),
        # custom analyses ---------------------------------------------------------------
        "data/exclusive_signal",
        "data/computeMatrix",
        "data/plotHeatmap",
        "data/intervene",
        # "data/consensus_igg",
        # run homer on CUT&RUN TFs
        expand("data/zip_HOMER/{method}_{condition}_{mark}.zip",
            method = all_methods,
            condition = ["H1D1"],
            mark = ["CTCF", "Sox2"]),
        # custom figures ----------------------------------------------------------------
        "data/figures-evaluate-consensus-counts/peak-distances",
        "data/figures-evaluate-consensus-counts/peak-counts",
        "data/figures-evaluate-consensus-counts/FRiP",
        "data/figures-evaluate-consensus-counts/roc_counts",
        "data/figures-evaluate-consensus-counts/pr_counts",
        "data/figures-evaluate-consensus-counts/consensus_peak_plots",
        "data/figures-evaluate-consensus-counts/exclusive-peaks",
        "data/figures-evaluate-consensus-counts/consensus-peak-standard-intersections.pdf",
        "data/consensus/consensus-standard-intersection/intersections-summary.txt",
        "data/figures-evaluate-consensus-counts/homer",
        # model evaluation --------------------------------------------------------------
        expand("data/evaluate_consensus_counts/{method}_{sample}.txt",
            method = all_methods,
            sample = list(config["STANDARDS"].keys()) )

# fastqc for each read
rule fastqc:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        html="data/fastqc/{read}.html",
        zip="data/fastqc/{read}_fastqc.zip"
    params:
        ""
    log:
        "data/logs/fastqc_{read}.log"
    threads: 4
    wrapper:
        "0.65.0/bio/fastqc"

rule fastq_screen:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        "data/fastq_screen/{read}_screen.txt"
    conda:
        "envs/fastq_screen.yml"
    log:
        "data/logs/fastq_screen_{read}.txt"
    threads: 8
    shell:
        "fastq_screen --aligner bowtie2 --threads {threads} --outdir data/fastq_screen "
        "--conf {config[FASTQ_SCREEN_CONFIG]} --force {input} > {log} 2>&1"

# align samples to genome
rule bowtie2:
    input:
        get_bowtie2_input
    output:
        "data/aligned/{sample}.bam"
    log:
        err="data/logs/bowtie2_{sample}.err"
    conda:
        "envs/align.yml"
    threads: 8
    shell:
        "bowtie2 --local --very-sensitive-local "
        "--no-unal --no-mixed --threads {threads} "
        "--no-discordant --phred33 "
        "-I 10 -X 700 -x {config[GENOME]} "
        "-1 {input[0]} -2 {input[1]} 2>{log.err} | samtools view -@ {threads} -Sbh - > {output}"

rule sort:
    input:
        "data/aligned/{sample}.bam"
    output: 
        temp("data/aligned/{sample}.sort.bam")
    conda:
        "envs/sambamba.yml"
    threads: 4
    log:
        "data/logs/sambamba_sort_{sample}.log"
    shell:
        "sambamba sort --tmpdir=data/aligned -t {threads} -o {output} {input} > {log} 2>&1"

rule markdup:
    input:
        rules.sort.output
    output:
        "data/markd/{sample}.sorted.markd.bam"
    conda:
        "envs/sambamba.yml"
    threads: 4
    log:
        "data/logs/sambamba_markdup_{sample}.log"
    shell:
        "sambamba markdup --tmpdir=data/markd -t {threads} {input} {output} > {log} 2>&1"

# remove reads in blacklist regions
rule banlist:
    input:
        rules.markdup.output
    output:
        "data/ban/{sample}.ban.sorted.markd.bam"
    conda:
        "envs/bedtools.yml"
    shell:
        "bedtools intersect -v -a {input} -b {config[BANLIST]} > {output}"

rule index:
    input:
        rules.banlist.output
    output:
        "data/ban/{sample}.ban.sorted.markd.bam.bai"
    conda:
        "envs/sambamba.yml"
    threads: 4
    log:
        "data/logs/samtools_index_{sample}.log"
    shell:
        "sambamba index -t {threads} {input} > {log} 2>&1"

rule tracks:
    input:
        rules.banlist.output,
        rules.index.output
    output:
        "data/tracks/{sample}.bw"
    conda:
        "envs/dtools.yml"
    threads:
        8
    shell:
        "bamCoverage -b {input[0]} -o {output} -p {threads} --binSize 10 --smoothLength 50 --normalizeUsing CPM"

rule fraglength:
    input:
        rules.banlist.output
    output:
        "data/ban/{sample}.sorted.markd.fraglen.tsv"
    conda:
        "envs/align.yml"
    shell:
        "src/fraglen-dist.sh {input} {output}"

rule fraglength_plot:
    input:
        expand("data/ban/{sample}.sorted.markd.fraglen.tsv", sample = samps)
    output:
        "data/ban/fraglen.html"
    run:
        pd.options.plotting.backend = "plotly"
        dfs = []
        for i in input:
            cond_marker = [os.path.basename(i).split(".")[0]]
            temp_df = pd.read_csv(i, sep = "\t", index_col = 0, names = cond_marker)
            dfs.append(temp_df)
        df = pd.concat(dfs, axis = 1)
        fraglen = df.plot()
        fraglen.update_layout( 
            title='Fragment Length Distribution', 
            xaxis_title='Fragment Length (bp)', 
            yaxis_title='Counts', 
            legend_title_text='Samples')
        fraglen.write_html(str(output))

rule preseq:
    input:
       rules.banlist.output
    output:
        "data/preseq/estimates_{sample}.txt"
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq c_curve -B -P -o {output} {input} > {log} 2>&1" 

rule preseq_lcextrap:
    input:
        rules.banlist.output
    output:
        "data/preseq/lcextrap_{sample}.txt"
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq lc_extrap -B -P -e 1000000000 -o {output} {input} > {log} 2>&1"

rule plotFinger:
    input:
        "data/ban/{sample}.ban.sorted.markd.bam", "data/ban/{sample}.ban.sorted.markd.bam.bai"
    output:
        "data/dtools/fingerprint_{sample}.tsv"
    conda:
        "envs/dtools.yml"
    log:
        "data/logs/fingerprint_{sample}.log"
    shell:
        "plotFingerprint -b {input[0]} --smartLabels --outRawCounts {output} > {log} 2>&1"

rule wget_gopeaks:
    output:
       "src/gopeaks"
    shell:
       "src/wget_gopeaks.sh {config[GOPEAKS_RELEASE_VERSION]}"

rule gopeaks:
    input:
        gopeaks = rules.wget_gopeaks.output,
        sample = "data/ban/{sample}.ban.sorted.markd.bam",
        index = "data/ban/{sample}.ban.sorted.markd.bam.bai",
        igg = get_igg
    output:
        "data/gopeaks/{sample}_peaks.bed"
    params:
        igg = gopeaks_igg,
        mindwidth = get_minwidth,
        step = get_step,
        slide = get_slide,
        mdist = get_mdist
    log:
        "data/logs/gopeaks_{sample}.log"
    shell:
        "{input.gopeaks} --bam {input.sample} {params.igg} "
        "{params.mdist} {params.mindwidth} "
        "{params.step} {params.slide} "
        "-o data/gopeaks/{wildcards.sample} > {log} 2>&1"
# input.igg requires the IgG bam file, even if treatment is IgG. however, params.igg will mask input.igg if treatment is IgG.
# so treatment file != control file for all samples.

rule consensus:
    input:
        group_reps
    output:
        "data/consensus/{method}_{condition}_{mark}.bed"
    params:
        consensus = consensus_params
    conda:
        "envs/bedtools.yml"
    shell:
        "cat {input} | sort -k1,1 -k2,2n | "
        "bedtools merge | bedtools intersect -a stdin -b {input} -c | "
        "awk -v OFS='\t' '$4 >= {params.consensus} {{print}}' | cut -f1-3 > {output} "

rule consensus_igg:
    input:
        expand("data/consensus/{method}_{condition}_{mark}.bed", zip,
            method = all_groups.method,
            condition = all_groups.condition,
            mark = all_groups.mark)
    output:
        directory("data/consensus_igg")
    conda:
        "envs/bedtools.yml"
    shell:
        "bash src/custom/consensus_igg.sh"

rule homer:
    input:
        "data/consensus/{method}_{condition}_{mark}.bed"
    output:
        directory("data/homer/{method}_{condition}_{mark}")
    wildcard_constraints:
        condition = "H1D.+"
    params:
        preparsed_dir = "/home/groups/MaxsonLab/indices/GRch38"
    conda:
        "envs/homer.yaml"
    threads: 8
    log:
        "data/logs/{method}_{condition}_{mark}.log"
    shell:
        "findMotifsGenome.pl {input} {config[FASTA]} {output} "
        "-size 200 -p {threads} -preparsedDir {params.preparsed_dir} > {log} 2>&1"

rule zip_homer:
    input:
        "data/homer/{method}_{condition}_{mark}"
    output:
        "data/zip_HOMER/{method}_{condition}_{mark}.zip"
    shell:
        "zip -r {output} {input}"

rule plot_homer:
    input:
        expand("data/zip_HOMER/{method}_{condition}_{mark}.zip",
            method = all_methods,
            condition = ["H1D1"],
            mark = ["CTCF", "Sox2"])
    output:
        directory("data/figures-evaluate-consensus-counts/homer")
    conda:
        "envs/plot.yml"
    script:
        "src/custom/plot_homer.R"

rule multiqc:
    input:
        directory("data/")
    output:
        "data/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yml"
    log:
        "data/logs/multiqc.log"
    shell:
        "multiqc -f -c src/multiqc_conf.yml -o data/multiqc data/ > {log} 2>&1"

# custom analyses ---------------------------------------------------------------------------------

rule peak_plot:
    input:
        expand("data/peak_counts/{method}_{sample}.bed",
            method = all_methods,
            sample = sample_noigg)
    output:
        directory("data/figures/peak_plots")
    conda:
        "envs/plot.yml"
    script:
        "src/custom/peak_plots.R"

# count basic peak statistics like peak counts, FRiP at the sample level.
# 'data/consensus' input ensures peak-calling is finished.
rule peak_characteristics:
    input:
        "data/consensus"
    output:
        directory("data/figures-evaluate-consensus-counts/peak-distances"),
        directory("data/figures-evaluate-consensus-counts/peak-counts"),
        directory("data/figures-evaluate-consensus-counts/FRiP"),
    conda:
        "envs/plot.yml"
    script:
        "src/custom/peak_characteristics.R"

# deeptools heatmap for all samples at consensus intervals.
rule heatmap:
    input:
        expand("data/tracks/{sample}.bw", sample = samps)
    output:
        directory("data/computeMatrix"),
        directory("data/plotHeatmap")
    conda:
        "envs/dtools.yml"
    shell:
        "bash src/custom/heatmap.sh"

rule exclusive_heatmaps:
    input:
        expand("data/tracks/{sample}.bw", sample = samps),
        "data/intervene"
    output:
        directory("data/exclusive-heatmaps")
    conda:
        "envs/dtools.yml"
    shell:
        "bash src/custom/exclusive_heatmaps.sh"

# venn diagram of consensus peaks
rule intervene:
    input:
        expand("data/consensus/{method}_{condition}_{mark}.bed", zip,
            method = list(all_groups.method),
            condition = list(all_groups.condition),
            mark = list(all_groups.mark)),
    output:
        directory("data/intervene")
    conda:
        "envs/intervene.yml"
    shell:
        "bash src/custom/intervene.sh"

# count reads at caller-exclusive peaks.
rule exclusive_counts:
    input:
        "data/intervene"
    output:
        directory("data/exclusive_signal")
    conda:
        "envs/bedtools.yml"
    shell:
        "bash src/custom/exclusive_signal.sh"

rule exclusive_counts_plot:
    input:
        "data/exclusive_signal"
    output:
        directory("data/figures-evaluate-consensus-counts/exclusive-counts")
    conda:
        "envs/bedtools.yml"
    shell:
        "bash src/custom/exclusive_signal.sh"

rule exclusive_annotation:
    input:
        "data/intervene"
    output:
        directory("data/figures-evaluate-consensus-counts/exclusive-peaks")
    conda:
        "envs/plot.yml"
    script:
        "src/custom/exclusive-peaks.R"

# count consensus peak intersections with promoters by method.
# promoter buckets usually 1kb downstream and 1-5kb upstream of a gene.
rule promoter_fragment:
    input:
        "data/consensus",
        config["GTF"]
    output:
        directory("data/fragments")
    conda:
        "envs/bedtools.yml"
    params:
        downstream = 1000,
        prefix = "hg38"
    shell:
        "bash src/custom/peak_fragmentation.sh -g {input[1]} -d {params.downstream} -o {params.prefix}"

rule promoter_fragment_plot:
    input:
        "data/fragments"
    output:
        "data/figures/promoter_fragment.png"
    conda:
        "envs/plot.yml"
    script:
        "src/custom/peak_fragmentation.R"

rule gene_ontology:
    input:
        "data/consensus",
        "data/intervene"
    output:
        directory("data/GO"),
        directory("data/figures/GO/exclusive")
    conda:
        "envs/gene_ontology.yml"
    script:
        "src/custom/go.R"

rule encode:
    input:
        "data/consensus"
    output:
        directory("data/encode")
    conda:
        "envs/bedtools.yml"
    shell:
        "bash src/custom/encode.sh"

rule plotEnrichment_consensus:
    input:
        "data/consensus/{method}_{condition}_{mark}.bed",
        "data/ban/{condition}_{replicate}_{mark}.ban.sorted.markd.bam"
    output:
        "data/plotEnrichment_consensus/{method}_{condition}_{replicate}_{mark}.png",
        "data/plotEnrichment_consensus/{method}_{condition}_{replicate}_{mark}.tsv"
    conda:
        "envs/dtools.yml"
    log:
        "data/logs/plotEnrichment_consensus_{method}_{condition}_{replicate}_{mark}.log"
    shell:
        "plotEnrichment -b {input[1]} --BED {input[0]} --regionLabels 'frip' "
        "--outRawCounts {output[1]} -o {output[0]} > {log} 2>&1"
