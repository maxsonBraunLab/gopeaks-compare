library(dplyr)
library(stringr)
library(conflicted)
library(GenomicRanges)
library(ggplot2)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

fragment_counts_list <- lapply(list.files("data/fragments", pattern = "*.bed", full.names = TRUE), function(x) {
  
  file = basename(x)
  upstream = str_split(file, "\\.")[[1]][2]

  read.table(x, col.names = c("chr", "start", "stop", "name", "score", "strand", "counts")) %>%
    mutate(file = file) %>%
    mutate(upstream = upstream) %>%
    mutate(meta = str_replace(file, ".\\d+.bed", "")) %>%
    tidyr::separate(meta, into = c("method", "condition", "mark"), sep = "_")
})

fragment_counts <- bind_rows(fragment_counts_list)
fragment_counts <- fragment_counts %>% 
  select(upstream, method, condition, mark, counts) %>%
  group_by(upstream, method, condition, mark, counts) %>%
  summarise(n = n()) %>%
  mutate(sample = paste(condition, mark, sep = "_"))

fragment_counts <- fragment_counts %>%
  filter(upstream %in% c(1000, 3000, 5000))

promoter_fragments <- ggplot(fragment_counts, aes(x = counts, y = n, color = method)) +
  geom_line(size = 2) +
  geom_point(size = 4) +
  facet_grid(sample ~ upstream, scales = "free_y") +
  scale_y_log10() +
  ggtitle("Counts of peaks between 1-kbp downstream and n-bp upstream of hg38 promoters") +
  xlab("Number of intersections at promoters") +
  ylab("Frequency of intersections")
promoter_fragments

ggsave(snakemake@output[[1]], promoter_fragments, width = 16, height = 16, units = "in", dpi = 600)
