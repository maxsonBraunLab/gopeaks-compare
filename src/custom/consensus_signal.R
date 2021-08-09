library(dplyr)
library(stringr)
library(conflicted)
library(GenomicRanges)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(scales)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

if (!dir.exists("data/figures/consensus_signal")) {
	dir.create("data/figures/consensus_signal")
} else {
	unlink("data/figures/consensus_signal", recursive = TRUE)
	dir.create("data/figures/consensus_signal")
}

# Peak size at consensus peaks --------------------------------------------------------------------

# peak size refers to height (counts of reads in that interval) and width.

peak_size_list <- lapply(list.files("data/consensus_signal", pattern = "*.bed", full.names = TRUE), function(x) {

  file = str_replace(basename(x), ".bed", "")
  method = str_split(file, "_")[[1]][1]
  condition = str_split(file, "_")[[1]][2]
  replicate = str_split(file, "_")[[1]][3]
  mark = str_split(file, "_")[[1]][4]
  
  read.table(x, col.names = c("chr", "start", "end", "counts")) %>%
    mutate(method = method) %>%
    mutate(condition = condition) %>%
    mutate(replicate = replicate) %>%
    mutate(mark = mark) %>%
    mutate(sample = paste(condition, replicate, mark, sep = "_")) %>%
    mutate(peaks = nrow(.))
})

peak_size <- bind_rows(peak_size_list) %>%
  mutate(width = end - start) %>%
  mutate(width = as.numeric(width)) %>%
  mutate(counts = counts + 1) # add one pseudo-count

peak_counts <- bind_rows(peak_size_list) %>%
  select(method, condition, replicate, mark, peaks) %>%
  distinct()

peak_size <- peak_size %>%
  inner_join(peak_counts)

peak_groups <- peak_size %>% group_by(condition, mark) %>% group_keys()
print(peak_groups)

for (i in 1:nrow(peak_groups)) {
  
  temp_condition <- as.character(peak_groups[i, "condition"])
  temp_mark <- as.character(peak_groups[i, "mark"])
  print(paste("Exporting consensus peak size for", temp_condition, temp_mark))

  temp_df <- peak_size %>% 
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)
  
  out_file <- paste0("data/figures/consensus_signal/", temp_condition, "_", temp_mark, ".pdf")
  plot_title <- paste(temp_condition, temp_mark, "Consensus Peak Size")
  
  temp_plot <- ggplot(temp_df, aes_string(x = "width", y = "counts", col = "method")) +
    geom_point() +
    geom_density_2d(color = "white") +
    geom_text(aes(x = 300, y = 10000, label = comma(peaks), size = 8), hjust = 0, col = "black") +
    facet_grid(sample ~ method) +
    scale_x_log10(labels = comma) +
    scale_y_log10(labels = comma) +
    scale_color_brewer(palette = "RdYlBu", direction = -1) +
    scale_size_continuous(guide = FALSE) +
    ggtitle(plot_title) +
    ylab("Counts") +
    xlab("Peak Width (bp)") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"),
          axis.text = element_text(color = "black"))

  print(temp_plot)
  ggsave(out_file, temp_plot, width = 16, height = 9, dpi = 600)
}


# Annotation per consensus peaks ------------------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks_df <- lapply(list.files("data/consensus", pattern="*.bed", full.names = TRUE), function(x) {

  method_condition_mark = str_replace(basename(x), ".bed", "")
  method = str_split(method_condition_mark, "_")[[1]][1]
  condition = str_split(method_condition_mark, "_")[[1]][2]
  mark = str_split(method_condition_mark, "_")[[1]][3]

  data.frame(method = method, condition = condition,
             mark = mark, sample = paste(condition, mark, sep = "_"),
             file = x)
})

peaks_df <- bind_rows(peaks_df) %>%
  arrange(sample)

peak_contents <- lapply(peaks_df$file, readPeakFile)
names(peak_contents) <- paste(peaks_df$method, peaks_df$condition, peaks_df$mark, sep = "_")

peak_annotation_list <- lapply(peak_contents, annotatePeak, TxDb = txdb, tssRegion = c(-3e3, 3e3))

peak_anno_df <- bind_rows( lapply(names(peak_annotation_list), function(x) {
  temp_sample <- peak_annotation_list[[x]]
  temp_sample@annoStat %>% 
    mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
    mutate(peak_counts = temp_sample@peakNum) %>%
    mutate(meta = x) %>%
    mutate(method = gsub("(.*)_(.*)_(.*)", "\\1", meta)) %>%
    mutate(condition = gsub("(.*)_(.*)_(.*)", "\\2", meta)) %>%
    mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", meta)) %>%
    mutate(sample = gsub("(.*)_(.*)_(.*)", "\\2_\\3", meta))
}) )

peak_anno_df <- peak_anno_df %>% arrange(condition, mark)
save.image()
p <- ggplot(peak_anno_df, aes(x = method, y = Frequency, fill = Feature, group = method)) +
  geom_col(position = "stack") +
  geom_text(aes(x = method, y = 25, label = peak_counts), color = "white", size = 4) +
  facet_wrap(sample ~ .) +
  ggtitle("Peak Annotation at Consensus Peaks") +
  scale_fill_brewer(palette = "Spectral") +
  theme(axis.text = element_text(color = "black"))
p

ggsave("data/figures/consensus_signal/consensus_annotation.pdf", p, width = 16, height = 9, units = "in", dpi = 600)