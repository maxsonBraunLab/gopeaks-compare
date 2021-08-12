#!/bin/bash

# download & compress public datasets ------------------------------------------------------------------------

# Cut&Tag paper by Henikoff lab using K562 cells + histone mark.

# H3K4me1
fasterq-dump -O data/raw -p --split-files SRR8383512
fasterq-dump -O data/raw -p --split-files SRR8383513
# H3K4me2
fasterq-dump -O data/raw -p --split-files SRR8383514
fasterq-dump -O data/raw -p --split-files SRR8383515
# H3K4me3
fasterq-dump -O data/raw -p --split-files SRR8383516
fasterq-dump -O data/raw -p --split-files SRR8383517
# H3K27Ac
fasterq-dump -O data/raw -p --split-files SRR8383507
fasterq-dump -O data/raw -p --split-files SRR8383508
# IgG
fasterq-dump -O data/raw -p --split-files SRR8435051
fasterq-dump -O data/raw -p --split-files SRR8435052

ls data/raw/*.gz | parallel "gzip {}"

# rename files ------------------------------------------------------------------------------------

# H3K4me1
mv data/raw/SRR8383512_1.fastq.gz data/raw/K562_1_H3K4me1_R1.fastq.gz
mv data/raw/SRR8383512_2.fastq.gz data/raw/K562_1_H3K4me1_R2.fastq.gz

mv data/raw/SRR8383513_1.fastq.gz data/raw/K562_2_H3K4me1_R1.fastq.gz
mv data/raw/SRR8383513_2.fastq.gz data/raw/K562_2_H3K4me1_R2.fastq.gz

# H3K4me2
mv data/raw/SRR8383514_1.fastq.gz data/raw/K562_1_H3K4me2_R1.fastq.gz
mv data/raw/SRR8383514_2.fastq.gz data/raw/K562_1_H3K4me2_R2.fastq.gz

mv data/raw/SRR8383515_1.fastq.gz data/raw/K562_2_H3K4me2_R1.fastq.gz
mv data/raw/SRR8383515_2.fastq.gz data/raw/K562_2_H3K4me2_R2.fastq.gz

# H3K4me3
mv data/raw/SRR8383516_1.fastq.gz data/raw/K562_1_H3K4me3_R1.fastq.gz
mv data/raw/SRR8383516_2.fastq.gz data/raw/K562_1_H3K4me3_R2.fastq.gz

mv data/raw/SRR8383517_1.fastq.gz data/raw/K562_2_H3K4me3_R1.fastq.gz
mv data/raw/SRR8383517_2.fastq.gz data/raw/K562_2_H3K4me3_R2.fastq.gz

# H3K27Ac
mv data/raw/SRR8383507_1.fastq.gz data/raw/K562_1_H3K27Ac_R1.fastq.gz
mv data/raw/SRR8383507_2.fastq.gz data/raw/K562_1_H3K27Ac_R2.fastq.gz

mv data/raw/SRR8383508_1.fastq.gz data/raw/K562_2_H3K27Ac_R1.fastq.gz
mv data/raw/SRR8383508_2.fastq.gz data/raw/K562_2_H3K27Ac_R2.fastq.gz

# IgG
mv data/raw/SRR8435051_1.fastq.gz data/raw/K562_1_IgG_R1.fastq.gz
mv data/raw/SRR8435051_2.fastq.gz data/raw/K562_1_IgG_R2.fastq.gz

mv data/raw/SRR8435052_1.fastq.gz data/raw/K562_2_IgG_R1.fastq.gz
mv data/raw/SRR8435052_2.fastq.gz data/raw/K562_2_IgG_R2.fastq.gz