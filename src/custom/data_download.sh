#!/bin/bash

# download & compress public datasets ------------------------------------------------------------------------

# Cut&Tag paper by Henikoff lab using K562 cells + histone mark.

# K562 H3K4me1
fasterq-dump -O data/raw -p --split-files SRR8383512
fasterq-dump -O data/raw -p --split-files SRR8383513
# K562 H3K4me2
fasterq-dump -O data/raw -p --split-files SRR8383514
fasterq-dump -O data/raw -p --split-files SRR8383515
# K562 H3K4me3
fasterq-dump -O data/raw -p --split-files SRR8383516
fasterq-dump -O data/raw -p --split-files SRR8383517
# K562 H3K27Ac
fasterq-dump -O data/raw -p --split-files SRR8383507
fasterq-dump -O data/raw -p --split-files SRR8383508
# K562 IgG
fasterq-dump -O data/raw -p --split-files SRR8435051
fasterq-dump -O data/raw -p --split-files SRR8435052
# K562 H3K27me3
fasterq-dump -O data/raw -p --split-files SRR8435047
# K562 CTCF
fasterq-dump -O data/raw -p --split-files SRR8435039
fasterq-dump -O data/raw -p --split-files SRR8435040
fasterq-dump -O data/raw -p --split-files SRR8435041

ls data/raw/*.fastq | parallel -j 8 "gzip {}"

# rename files ------------------------------------------------------------------------------------

# K562 H3K4me1
mv data/raw/SRR8383512_1.fastq.gz data/raw/K562_1_H3K4me1_R1.fastq.gz
mv data/raw/SRR8383512_2.fastq.gz data/raw/K562_1_H3K4me1_R2.fastq.gz

mv data/raw/SRR8383513_1.fastq.gz data/raw/K562_2_H3K4me1_R1.fastq.gz
mv data/raw/SRR8383513_2.fastq.gz data/raw/K562_2_H3K4me1_R2.fastq.gz

# K562 H3K4me2
mv data/raw/SRR8383514_1.fastq.gz data/raw/K562_1_H3K4me2_R1.fastq.gz
mv data/raw/SRR8383514_2.fastq.gz data/raw/K562_1_H3K4me2_R2.fastq.gz

mv data/raw/SRR8383515_1.fastq.gz data/raw/K562_2_H3K4me2_R1.fastq.gz
mv data/raw/SRR8383515_2.fastq.gz data/raw/K562_2_H3K4me2_R2.fastq.gz

# K562 H3K4me3
mv data/raw/SRR8383516_1.fastq.gz data/raw/K562_1_H3K4me3_R1.fastq.gz
mv data/raw/SRR8383516_2.fastq.gz data/raw/K562_1_H3K4me3_R2.fastq.gz

mv data/raw/SRR8383517_1.fastq.gz data/raw/K562_2_H3K4me3_R1.fastq.gz
mv data/raw/SRR8383517_2.fastq.gz data/raw/K562_2_H3K4me3_R2.fastq.gz

# K562 H3K27Ac
mv data/raw/SRR8383507_1.fastq.gz data/raw/K562_1_H3K27Ac_R1.fastq.gz
mv data/raw/SRR8383507_2.fastq.gz data/raw/K562_1_H3K27Ac_R2.fastq.gz

mv data/raw/SRR8383508_1.fastq.gz data/raw/K562_2_H3K27Ac_R1.fastq.gz
mv data/raw/SRR8383508_2.fastq.gz data/raw/K562_2_H3K27Ac_R2.fastq.gz

# K562 IgG
mv data/raw/SRR8435051_1.fastq.gz data/raw/K562_1_IgG_R1.fastq.gz
mv data/raw/SRR8435051_2.fastq.gz data/raw/K562_1_IgG_R2.fastq.gz

mv data/raw/SRR8435052_1.fastq.gz data/raw/K562_2_IgG_R1.fastq.gz
mv data/raw/SRR8435052_2.fastq.gz data/raw/K562_2_IgG_R2.fastq.gz

# K562 H3K27me3
mv data/raw/SRR8435047_1.fastq.gz data/raw/K562_1_H3K27me3_R1.fastq.gz
mv data/raw/SRR8435047_2.fastq.gz data/raw/K562_1_H3K27me3_R2.fastq.gz

# K562 CTCF
mv data/raw/SRR8435039_1.fastq.gz data/raw/K562_1_CTCF_R1.fastq.gz
mv data/raw/SRR8435039_2.fastq.gz data/raw/K562_1_CTCF_R2.fastq.gz

mv data/raw/SRR8435040_1.fastq.gz data/raw/K562_2_CTCF_R1.fastq.gz
mv data/raw/SRR8435040_2.fastq.gz data/raw/K562_2_CTCF_R2.fastq.gz

mv data/raw/SRR8435041_1.fastq.gz data/raw/K562_3_CTCF_R1.fastq.gz
mv data/raw/SRR8435041_2.fastq.gz data/raw/K562_3_CTCF_R2.fastq.gz

# CUT&RUN -----------------------------------------------------------------------------------------

# source: Pioneer Factor-Nucleosome Binding Events during Differentiation Are Motif Encoded. Meers 2019.

# H1 Sox2 Day 1
fasterq-dump -O data/raw -p --split-files SRR8748855
fasterq-dump -O data/raw -p --split-files SRR8748856

# H1 CTCF Day 1
fasterq-dump -O data/raw -p --split-files SRR8748830
fasterq-dump -O data/raw -p --split-files SRR8748831

# H1 IgG Day 1
fasterq-dump -O data/raw -p --split-files SRR8748845

# H1 FoxA1 Day 5
fasterq-dump -O data/raw/ -p --split-files SRR8748811
fasterq-dump -O data/raw/ -p --split-files SRR8748812

# H1 IgG Day 5
fasterq-dump -O data/raw/ -p --split-files SRR8748821

ls data/raw/*.fastq | parallel -j 8 "gzip {}"

mv data/raw/SRR8748855_1.fastq.gz data/raw/H1D1_1_Sox2_R1.fastq.gz
mv data/raw/SRR8748855_2.fastq.gz data/raw/H1D1_1_Sox2_R2.fastq.gz

mv data/raw/SRR8748856_1.fastq.gz data/raw/H1D1_2_Sox2_R1.fastq.gz
mv data/raw/SRR8748856_2.fastq.gz data/raw/H1D1_2_Sox2_R2.fastq.gz

mv data/raw/SRR8748830_1.fastq.gz data/raw/H1D1_1_CTCF_R1.fastq.gz
mv data/raw/SRR8748830_2.fastq.gz data/raw/H1D1_1_CTCF_R2.fastq.gz

mv data/raw/SRR8748831_1.fastq.gz data/raw/H1D1_2_CTCF_R1.fastq.gz
mv data/raw/SRR8748831_2.fastq.gz data/raw/H1D1_2_CTCF_R2.fastq.gz

mv data/raw/SRR8748845_1.fastq.gz data/raw/H1D1_1_IgG_R1.fastq.gz
mv data/raw/SRR8748845_2.fastq.gz data/raw/H1D1_1_IgG_R2.fastq.gz

mv data/raw/SRR8748811_1.fastq.gz data/raw/H1D5_1_FoxA1_R1.fastq.gz
mv data/raw/SRR8748811_2.fastq.gz data/raw/H1D5_1_FoxA1_R2.fastq.gz

mv data/raw/SRR8748812_1.fastq.gz data/raw/H1D5_2_FoxA1_R1.fastq.gz
mv data/raw/SRR8748812_2.fastq.gz data/raw/H1D5_2_FoxA1_R2.fastq.gz

mv data/raw/SRR8748821_1.fastq.gz data/raw/H1D5_1_IgG_R1.fastq.gz
mv data/raw/SRR8748821_2.fastq.gz data/raw/H1D5_1_IgG_R2.fastq.gz