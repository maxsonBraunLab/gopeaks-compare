library(dplyr)
library(stringr)
library(conflicted)
library(GenomicRanges)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

if (!dir.exists("data/figures/exclusive-peaks")) {
	dir.create("data/figures/exclusive-peaks")
} else {
	unlink("data/figures/exclusive-peaks", recursive = TRUE)
	dir.create("data/figures/exclusive-peaks")
}

# Distribution of caller-exclusive-peaks counts ---------------------------------------------------
peak_size_list <- lapply(list.files("data/exclusive-peaks", pattern = "*.bed", full.names = TRUE), function(x) {

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
	mutate(counts = counts + 1) %>%
	mutate(norm_counts = counts / width * 1000)

exclusive_counts <- ggplot(peak_size, aes(x = method, y = norm_counts, fill = method, col = method)) +
            geom_violin() +
            geom_text(aes(x = method, y = 1, label = peaks), size = 4, col = "black") +
            facet_wrap(. ~ sample) +
            scale_y_log10() + 
            theme(axis.text = element_text(color = "black"))
exclusive_counts

ggsave("data/figures/exclusive-peaks/count_violin.pdf", exclusive_counts, width = 9, height = 16, dpi = 600)

# exclusive peak plots -----------------------------------------------------------------------------
peak_groups <- peak_size %>% group_by(condition, mark) %>% group_keys()
print(peak_groups)

for (i in 1:nrow(peak_groups)) {
  
  temp_condition <- as.character(peak_groups[i, "condition"])
  temp_mark <- as.character(peak_groups[i, "mark"])
  print(paste("Exporting exclusive peak plots for", temp_condition, temp_mark))

  temp_df <- peak_size %>% 
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)
  
  out_file <- paste0("data/figures/exclusive-peaks/", temp_condition, "_", temp_mark, ".pdf")
  plot_title <- paste(temp_condition, temp_mark, "Exclusive Peak Size")
  
  temp_plot <- ggplot(temp_df, aes_string(x = "width", y = "counts", col = "method")) +
    geom_point() +
    geom_density_2d(color = "black") +
    geom_text(aes_string(x = 300, y = 10000, label = "peaks", size = 8), hjust = 0, col = "black") +
    facet_grid(sample ~ method) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c("gopeaks" = "#2C7BB6", "macs2" = "#ABD9E9", "seacr-relaxed" = "#D7191C")) +
    scale_size_continuous(guide = FALSE) +
    ggtitle(plot_title) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"), 
    axis.text = element_text(color = "black"))

  print(temp_plot)
  ggsave(out_file, temp_plot, width = 16, height = 9, dpi = 600)
}






















# because we reach into bed files with counts information and replicate-level info, this is built
# differently than consensus_signal. my bad for the discontinuity.
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# exclusive_peak_list <- lapply(list.files("data/exclusive-peaks", pattern="*.bed", full.names = TRUE), function(x) {

#   method_condition_replicate_mark = str_replace(basename(x), ".bed", "")
#   method = str_split(method_condition_replicate_mark, "_")[[1]][1]
#   condition = str_split(method_condition_replicate_mark, "_")[[1]][2]
#   replicate = str_split(method_condition_replicate_mark, "_")[[1]][3]
#   mark = str_split(method_condition_replicate_mark, "_")[[1]][4]

#   data.frame(method = method, condition = condition, replicate = replicate,
#              mark = mark, method_condition_mark = paste(method, condition, mark, sep = "_"),
#              file = x)
# })

# peaks_df <- bind_rows(exclusive_peak_list) %>%
# 	arrange(method_condition_mark) %>%
# 	filter(replicate == 1) # since these are shared exclusive peaks, only take from replicate 1.

# peak_contents <- lapply(peaks_df$file, readPeakFile)
# names(peak_contents) <- peaks_df$method_condition_mark

# peak_annotation_list <- lapply(peak_contents, annotatePeak, TxDb = txdb, tssRegion = c(-3e3, 3e3))

# peak_anno_df <- bind_rows( lapply(names(peak_annotation_list), function(x) {
#   temp_sample <- peak_annotation_list[[x]]
#   temp_sample@annoStat %>%
#     mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
#     mutate(peak_counts = temp_sample@peakNum) %>%
#     mutate(meta = x) %>%
#     mutate(method = gsub("(.*)_(.*)_(.*)", "\\1", meta)) %>%
#     mutate(condition = gsub("(.*)_(.*)_(.*)", "\\2", meta)) %>%
#     mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", meta)) %>%
#     mutate(sample = gsub("(.*)_(.*)_(.*)", "\\2_\\3", meta)) %>%
#     mutate(meta = NULL)
# }) )

# peak_anno_df <- peak_anno_df %>% arrange(condition, mark)

# p <- ggplot(peak_anno_df, aes(x = method, y = Frequency, fill = Feature, group = method)) +
#   geom_col(position = "stack") +
#   geom_text(aes(x = method, y = 25, label = peak_counts), color = "white", size = 4) +
#   facet_wrap(sample ~ .) +
#   ggtitle("Peak Annotation at Peak-Caller Exclusive Peaks") +
#   scale_fill_brewer(palette = "Spectral")
# p

# ggsave("data/figures/exclusive-peaks/exclusive_annotation.pdf", p, width = 16, height = 9, units = "in", dpi = 600)

