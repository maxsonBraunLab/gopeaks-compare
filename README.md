# cutTag-pipeline

Snakemake pipeline for Cut&amp;Tag analysis 

# Setup

1. Configure the project directory

```
#clone to your local called 'my-project'
git clone https://github.com/maxsonBraunLab/cutTag-pipeline.git my-project

#create directory for your fastq files
cd my-project
mkdir -p data/raw

# link fastqs to data/raw 
ln -s /path/to/fastq/files/* data/raw

# make scripts executable
chmod +x src/*.py src/*.sh
```


Rename all samples in data/raw to simplify sample information, for example:

This file:
LIB200706TB_M6Q3_RBP1_S93_L001_R1_001.fastq.gz

Is renamed to:
M6Q3_RBP1_S93_R1.fastq.gz


2. Edit configuration files 

Edit runtime configuration in the file src/config.yml:

- Specify the path to the bowtie2 index for the genome you are aligning to.
- Specify the path to the bowtie2 indices for the fastq_screen genomes.

The pipeline detects samples in the subdirectory data/raw with the following assumptions:

 - Paired end reads
 - Read 1 and 2 are designated using "_R1", and "_R2"
 - the epigenetic mark label is the second split of the sample name by _ delimeter. For example M6C3_4Me1_S12_R2.fastq.gz will have the {mark} wildcard set to 4me1. This effects the output files from the calculation of counts tables.

Edit the DESEq2 configuration file: src/deseq2_metadata.csv

 - should have two columns labeled "sample", and "condition"
 - the sample column corresponds to the replicates of the given condition, and should be the same as the first split of the raw file: e.g. M6C3_4Me1_S12_R2.fastq.gz will have "sample" equal to M6C3.
 - the condition should be the name for each sample condition, and doees not have to come from the file name.
 
 The file src/deseq2_metadata is populated with the following example data:
 
 ```
sample,condition
M6C1,GSKQUIZ
M6C2,GSKQUIZ
M6C3,GSKQUIZ
M6D1,DMSO
M6D2,DMSO
M6D3,DMSO
M6G1,GSK
M6G2,GSK
M6G3,GSK
M6Q1,QUIZ
M6Q2,QUIZ
M6Q3,QUIZ
 ```



# Execution

Setup snakemake profile to run on compute cluster:

SLURM: follow [these instructions](https://github.com/Snakemake-Profiles/slurm)

Example of how to run pipeline

`
snakemake --use-conda --profile slurm -j 60 --latency-wait 60 --cluster-config src/cluster.yml --latency-wait 200 
`

A "dry-run" can be accomplished to see what files would be generated by using the command:

`
snakemake -nrp
`

If the output of this command is not green, there is something wrong.


# Output


Below is an explanation of each output directory:

```
aligned - sorted and aligned sample bam files
callpeaks - the output of callpeaks.py with peak files and normalized bigwigs for each sample.
counts - the raw counts for each mark over a consensus peak list
deseq2 - the results of deseq2 fore each counts table (by mark)
dtools - fingerprint plot data for multiqc to use
fastqc - fastqc results
fastq_screen - fastq_screen results
logs - runtime logs for each snakemake rule
markd - duplicate marked bam files
multiqc - contains the file multiqc_report.html with a lot of QC info about the experiment.
plotEnrichment - FRIP statistics for each mark
preseq - library complexity data for multiqc report
```

## Deseq2 outputs

Each mark should have the following output files:

```
"data/deseq2/{mark}/{mark}-rld-pca.png" - PCA of counts after reguarlized log transformation. 
"data/deseq2/{mark}/{mark}-vsd-pca.png" - PCA of counts after variance stabilizing transformation.
"data/deseq2/{mark}/{mark}-normcounts.csv" - normalized count for each sample in each consensus peak.
"data/deseq2/{mark}/{mark}-lognormcounts.csv" - log2 normalized counts for each sample in each consensus peak.
"data/deseq2/{mark}/{mark}-rld.png", - the sdVsMean plot using regularized log transformation.
"data/deseq2/{mark}/{mark}-vsd.png" - the sdVsMean plot using variance stabilizing transformation.
"data/deseq2/{mark}/{mark}-vsd-dist.png" - the sample distance matrix after variance stabilizing transformation.
"data/deseq2/{mark}/{mark}-rld-dist.png" - the sample distance matrix using regularized log transformation.
"data/deseq2/{mark}/{mark}-dds.rds" - the R object with the results of running the DEseq() function.
```

For each contrast, the differentially expressed genes are written to a file ending in `-diffexp.tsv` as well as those with an adjusted p-value less than 0.05 with the extension `-sig05-diffexp.tsv`. A summary of the results usign an alpha of 0.05 is also written to a file with the extension `-sig05-diffexp-summary.txt`.

