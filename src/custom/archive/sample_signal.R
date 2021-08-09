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

peak_size_list <- lapply(list.files("data/sample_signal", pattern = "*.bed", full.names = TRUE), function(x) {
  
  method_condition_replicate_mark = str_replace(basename(x), ".bed", "")
  
  method = str_split(method_condition_replicate_mark, "_")[[1]][1]
  condition = str_split(method_condition_replicate_mark, "_")[[1]][2]
  replicate = str_split(method_condition_replicate_mark, "_")[[1]][3]
  mark = str_split(method_condition_replicate_mark, "_")[[1]][4]
  
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
  mutate(counts = counts + 1) %>% # add one pseudo-count
  select(method, condition, replicate, mark, width, counts) %>%
  mutate(sample = paste(condition, replicate, mark, sep = "_"))

peak_groups <- peak_size %>% group_by(condition, mark) %>% group_keys()
print(peak_groups)

peak_counts <- bind_rows(peak_size_list) %>%
  select(method, condition, replicate, mark, peaks) %>%
  distinct()

peak_size <- peak_size %>%
  inner_join(peak_counts)

if (!dir.exists("data/figures/sample_signal")) {
  dir.create('data/figures/sample_signal')
}

for (i in 1:nrow(peak_groups)) {
  
  temp_condition <- as.character(peak_groups[i, "condition"])
  temp_mark <- as.character(peak_groups[i, "mark"])
  
  temp_df <- peak_size %>% 
    filter(condition == eval(parse(text = "temp_condition"))) %>%
    filter(mark == eval(parse(text = "temp_mark")))
  
  out_file <- paste0("data/figures/sample_signal/", temp_condition, "_", temp_mark, ".png")
  plot_title <- paste(temp_condition, temp_mark, "peak sizes")
  
  temp_plot <- ggplot(temp_df, aes_string(x = "width", y = "counts", col = "method")) +
    geom_point() +
    geom_density_2d(color = "white") +
    facet_grid(sample ~ method) +
    scale_x_log10() +
    scale_y_log10() +
    geom_text(aes(x = 300, y = 10000, label = peaks), size = 6, col = "black") +
    labs(col = "# of Counts") +
    ggtitle(plot_title) +
    scale_color_brewer(palette = "Spectral")
  
  print(temp_plot)
  ggsave(out_file, temp_plot, width = 16, height = 9, dpi = 600)
  
}
