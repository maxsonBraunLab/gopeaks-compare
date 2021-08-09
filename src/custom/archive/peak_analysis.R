library(dplyr)
library(stringr)
library(conflicted)
library(GenomicRanges)
library(ggplot2)
library(stringr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")

if (!dir.exists("data/figures")) {
  dir.create("data/figures")
}

# REPLICATE-LEVEL ANALYSIS ------------------------------------------------------------------------

# COUNT PEAKS PER SAMPLE

macs2 <- snakemake@input[["macs2_peaks"]]
seacr_stringent <- snakemake@input[["seacr_stringent_peaks"]]
seacr_relaxed <- snakemake@input[["seacr_relaxed_peaks"]]
gopeaks <- snakemake@input[["gopeaks_peaks"]]

# import macs2
macs2_df <- data.frame(sample = str_replace(basename(macs2), "_peaks.*", ""),
                       method = "macs2",
                       file = macs2)

# import seacr
seacr_relaxed_df <- data.frame(sample = str_replace(basename(seacr_relaxed), ".relaxed.bed", ""),
                               method = "seacr_relaxed",
                               file = seacr_relaxed)
seacr_stringent_df <- data.frame(sample = str_replace(basename(seacr_stringent), ".stringent.bed", ""),
                               method = "seacr_stringent",
                               file = seacr_stringent)

# import gopeaks
gopeaks_df <- data.frame(sample = str_replace(basename(gopeaks), ".bed", ""),
                               method = "gopeaks",
                               file = gopeaks)
peak_counts_df <- bind_rows(macs2_df, seacr_relaxed_df, seacr_stringent_df, gopeaks_df)

# count peak files per sample
peak_counts <- lapply(peak_counts_df$file, function(x) {
  if (file.size(x) == 0) {
    return(0)
  } else {
    print(x)
    nrow(read.table(x, header = FALSE))
  }
})

peak_counts_df <- peak_counts_df %>% 
              mutate(counts = unlist(peak_counts)) %>%
              mutate(file = NULL) %>%
              mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", sample))

peak_counts_plot <- reshape2::melt(peak_counts_df) %>% 
  ggplot(aes(x = sample, y = value, fill = method)) + 
    geom_col(position = "dodge") + 
    facet_wrap(mark ~ ., scales = "free") +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x = element_text(hjust = 1)) +
    ggtitle("Peak counts with Kasumi_1_IgG")
peak_counts_plot

igg_peaks <- peak_counts_df %>% 
  filter(mark == "IgG") %>%
  reshape2::melt(.) %>%
    ggplot(aes(x = sample, y = log10(value), fill = method)) + 
      geom_col(position = "dodge") +
      ylab("log10(peak counts)")
igg_peaks

ggsave("data/figures/igg-peaks.png", igg_peaks, width = 16, height = 9, units = "in")
ggsave("data/figures/peak-counts.png", peak_counts_plot, width = 16, height = 9, units = "in")

# COUNT FRiP PER SAMPLE
macs2 <- snakemake@input[["macs2_frip"]]
seacr_stringent <- snakemake@input[["seacr_relaxed_frip"]]
seacr_relaxed <- snakemake@input[["seacr_stringent_frip"]]
gopeaks <- snakemake@input[["gopeaks_frip"]]

# import macs2
macs2_df <- data.frame(sample = str_replace(basename(macs2), "macs2_", ""),
                       method = "macs2",
                       file = macs2)
# import seacr
seacr_relaxed_df <- data.frame(sample = str_replace(basename(seacr_relaxed), "seacr_relaxed_", ""),
                               method = "seacr_relaxed",
                               file = seacr_relaxed)
seacr_stringent_df <- data.frame(sample = str_replace(basename(seacr_stringent), "seacr_stringent_", ""),
                                 method = "seacr_stringent",
                                 file = seacr_stringent)
# import gopeaks
gopeaks_df <- data.frame(sample = str_replace(basename(gopeaks), "gopeaks_", ""),
                         method = "gopeaks",
                         file = gopeaks)
frip_df <- bind_rows(macs2_df, seacr_relaxed_df, seacr_stringent_df, gopeaks_df)

# extract FRiP % per sample.
sample_FRiP <- lapply(frip_df$file, function(x) {
  if (file.size(x) == 0) {
    return(0)
  } else {
    print(x)
    read.table(x, header = TRUE)[["percent"]]
  }
})

# add frip column, emove file extensions for samples.
frip_df <- frip_df %>% 
  mutate(file = NULL) %>%
  mutate(sample = str_replace(sample, ".tsv", "")) %>%
  mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", sample)) %>%
  mutate(FRiP = unlist(sample_FRiP))

head(frip_df)

frip <- reshape2::melt(frip_df) %>%
  ggplot(aes(x = sample, y = value, fill = method)) +
    geom_col(position = "dodge") +
    facet_wrap(mark ~ ., scales = "free", nrow = 2) +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x = element_text(hjust = 1))
frip

ggsave("data/figures/FRiP.png", frip, width = 16, height = 9, units = "in")

# ANNOTATION PER SAMPLE
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
macs2 <- snakemake@input[["macs2_peaks"]]
seacr_stringent <- snakemake@input[["seacr_stringent_peaks"]]
seacr_relaxed <- snakemake@input[["seacr_relaxed_peaks"]]
gopeaks <- snakemake@input[["gopeaks_peaks"]]

# import macs2
macs2_df <- data.frame(sample = str_replace(basename(macs2), "_peaks.*", ""),
                       method = "macs2",
                       file = macs2)
# import seacr
seacr_relaxed_df <- data.frame(sample = str_replace(basename(seacr_relaxed), ".relaxed.bed", ""),
                               method = "seacr-relaxed",
                               file = seacr_relaxed)
seacr_stringent_df <- data.frame(sample = str_replace(basename(seacr_stringent), ".stringent.bed", ""),
                                 method = "seacr-stringent",
                                 file = seacr_stringent)

# import gopeaks
gopeaks_df <- data.frame(sample = str_replace(basename(gopeaks), ".bed", ""),
                         method = "gopeaks",
                         file = gopeaks)

# bind rows, filter away IgG peaks, arrange by sample and not method
anno_df <- bind_rows(macs2_df, seacr_relaxed_df, seacr_stringent_df, gopeaks_df) %>% 
              filter(str_detect(sample, "IgG", negate = TRUE)) %>%
              arrange(sample)

# filter bed files without peaks.
is_empty_file <- lapply(anno_df$file, function(x) {
  if (file.size(x) == 0) {
    return(1)
  } else {
    return(0)
  }
})

anno_df <- anno_df %>% 
              mutate(is_empty = unlist(is_empty_file)) %>% 
              filter(is_empty == 0)

# read peak file and annotate peaks
peak_list <- lapply(anno_df$file, readPeakFile)
names(peak_list) <- paste(anno_df$sample, anno_df$method, sep = "_")

peak_list_annotation <- lapply(peak_list, annotatePeak,
                                  TxDb = txdb, tssRegion = c(-3e3, 3e3))

# plot cleaner results.
peak_anno_df <- bind_rows( lapply(names(peak_list_annotation), function(x) {
  temp_sample <- peak_list_annotation[[x]]
  temp_sample@annoStat %>% 
    mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
    mutate(peak_counts = temp_sample@peakNum) %>%
    mutate(sample_method = x) %>%
    mutate(sample = gsub("(.*)_(.*)_(.*)_(.*)", "\\1_\\2_\\3", sample_method)) %>%
    mutate(samp = gsub("(.*)_(.*)_(.*)_(.*)", "\\1_\\2", sample_method)) %>%
    mutate(method = gsub("(.*)_(.*)_(.*)_(.*)", "\\4", sample_method)) %>%
    mutate(mark = gsub("(.*)_(.*)_(.*)_(.*)", "\\3", sample_method))
}) )
peak_anno_df

peak_anno_plot <- ggplot(peak_anno_df, aes(x = method, y = Frequency, fill = Feature, group = method, alpha = log10(peak_counts))) +
  geom_col(position = "stack") +
  geom_text(aes(x = method, y = 25, label = peak_counts), color = "white") +
  facet_wrap(mark + sample ~ ., scales = "free", nrow = 4) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x = element_text(hjust = 1)) +
  ggtitle("Peak Annotation by Sample and Method (Proportions)") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Peak annotation per sample at called")
peak_anno_plot

ggsave("data/figures/peak-annotation.png", peak_anno_plot, width = 18, height = 24, units = "in")

# CONSENSUS-LEVEL ANALYSIS ------------------------------------------------------------------------

# COUNT CONSENSUS PEAKS
consensus_counts_df <- data.frame(file = snakemake@input[["all_consensus"]]) %>%
  mutate(sample = str_replace(basename(file), ".bed", "")) %>%
  mutate(method = gsub("(.*)_(.*)_(.*)", "\\1", sample)) %>%
  mutate(condition = gsub("(.*)_(.*)_(.*)", "\\2", sample)) %>%
  mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", sample)) %>%
  mutate(sample = gsub("(.*)_(.*)_(.*)", "\\2_\\3", sample)) %>%
  mutate(counts = unlist(
    lapply(file, function(x) {nrow(read.table(x, header = FALSE))})
  ))
head(consensus_counts_df)

consensus_counts_plot <- ggplot(consensus_counts_df, aes(x = method, y = counts, fill = mark)) +
  geom_col(position = "dodge") +
  facet_wrap(condition ~ .) +
  ggtitle("Peak Counts with Kasumi_1_IgG at Consensus Peaks")

ggsave("data/figures/consensus-counts.png", consensus_counts_plot, width = 16, height = 9, units = "in", dpi = 600)

# ANNOTATION
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
macs2 <- snakemake@input[["macs2_consensus"]]
seacr_relaxed <- snakemake@input[["seacr_relaxed_consensus"]]
seacr_stringent <- snakemake@input[["seacr_stringent_consensus"]]
gopeaks <- snakemake@input[["gopeaks_consensus"]]

# import macs2
macs2_df <- data.frame(sample = str_replace(basename(macs2), ".bed", ""),
                       method = "macs2",
                       file = macs2)
# import seacr
seacr_relaxed_df <- data.frame(sample = str_replace(basename(seacr_relaxed), ".bed", ""),
                               method = "seacr-relaxed",
                               file = seacr_relaxed)
seacr_stringent_df <- data.frame(sample = str_replace(basename(seacr_stringent), ".bed", ""),
                                 method = "seacr-stringent",
                                 file = seacr_stringent)

# import gopeaks
gopeaks_df <- data.frame(sample = str_replace(basename(gopeaks), ".bed", ""),
                         method = "gopeaks",
                         file = gopeaks)

# bind rows, filter away IgG peaks, arrange by sample and not method
consensus_df <- bind_rows(macs2_df, seacr_relaxed_df, seacr_stringent_df, gopeaks_df) %>%
  filter(str_detect(sample, "IgG", negate = TRUE)) %>%
  arrange(sample) %>%
  mutate(sample = str_replace(sample, "seacr_relaxed", "seacr-relaxed")) %>%
  mutate(sample = str_replace(sample, "seacr_stringent", "seacr-stringent"))

# filter bed files without peaks.
is_empty_file <- lapply(consensus_df$file, function(x) {
  if (file.size(x) == 0) {
    return(1)
  } else {
    return(0)
  }
})

consensus_df <- consensus_df %>%
  mutate(is_empty = unlist(is_empty_file)) %>%
  filter(is_empty == 0)

peak_list <- lapply(consensus_df$file, readPeakFile)
names(peak_list) <- consensus_df$sample

peak_list_annotation <- lapply(peak_list, annotatePeak,
                               TxDb = txdb, tssRegion = c(-3e3, 3e3)) # mclapply for linux OS

peak_anno_df <- bind_rows( lapply(names(peak_list_annotation), function(x) {
  temp_sample <- peak_list_annotation[[x]]
  temp_sample@annoStat %>% 
    mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
    mutate(peak_counts = temp_sample@peakNum) %>%
    mutate(meta = x) %>%
    mutate(method = gsub("(.*)_(.*)_(.*)", "\\1", meta)) %>%
    mutate(condition = gsub("(.*)_(.*)_(.*)", "\\2", meta)) %>%
    mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", meta)) %>%
    mutate(sample = gsub("(.*)_(.*)_(.*)", "\\2_\\3", meta))
}) )

head(peak_anno_df)

p <- ggplot(peak_anno_df, aes(x = method, y = Frequency, fill = Feature, group = method, alpha = peak_counts)) +
  geom_col(position = "stack") +
  geom_text(aes(x = method, y = 25, label = peak_counts), color = "white") +
  facet_wrap(mark + sample ~ ., scales = "free", nrow = 4) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x = element_text(hjust = 1)) +
  ggtitle("Peak Annotation per sample at consensus peaks") +
  scale_fill_brewer(palette = "Spectral")
p

ggsave("data/figures/consensus-annotation.png", p, width = 18, height = 24, units = "in", dpi = 600)