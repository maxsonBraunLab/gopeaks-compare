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

if (!dir.exists("data/figures-evaluate-consensus-counts/consensus_peak_plots")) {
	dir.create("data/figures-evaluate-consensus-counts/consensus_peak_plots")
} else {
	dir.create("data/figures-evaluate-consensus-counts/consensus_peak_plots")
}

# Peak size at consensus peaks --------------------------------------------------------------------

# peak size refers to height (counts of reads in that interval) and width.

peak_size_list <- lapply(list.files("data/consensus_counts", pattern = "*.bed", full.names = TRUE), function(x) {

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
  
  out_pdf <- paste0("data/figures-evaluate-consensus-counts/consensus_peak_plots/", temp_condition, "_", temp_mark, ".pdf")
  out_png <- paste0("data/figures-evaluate-consensus-counts/consensus_peak_plots/", temp_condition, "_", temp_mark, ".png")
  plot_title <- paste(temp_condition, temp_mark, "Consensus Peak Size")
  
  temp_plot <- ggplot(temp_df, aes_string(x = "width", y = "counts", col = "method")) +
    geom_point() +
    geom_density_2d(color = "white") +
    geom_text(aes(x = 300, y = 10000, label = comma(peaks), size = 8), hjust = 0, col = "black", check_overlap = TRUE) +
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
  ggsave(out_pdf, temp_plot, width = 16, height = 9, dpi = 600)
  ggsave(out_png, temp_plot, width = 16, height = 9, dpi = 600)
}

