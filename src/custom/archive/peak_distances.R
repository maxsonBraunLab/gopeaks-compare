library(dplyr)
library(stringr)
library(conflicted)
library(GenomicRanges)
library(ggplot2)
# library(plyranges)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

if (!dir.exists("data/figures/peak_distances")) {
    dir.create("data/figures/peak_distances")
} else {
    unlink("data/figures/peak_distances", recursive = TRUE)
    dir.create("data/figures/peak_distances")
}

all_peaks <- lapply(list.files("data/consensus", pattern = "*.bed", full.names = TRUE), function(x) {

  file = str_replace(basename(x), ".bed", "")
  method = str_split(file, "_")[[1]][1]
  condition = str_split(file, "_")[[1]][2]
  mark = str_split(file, "_")[[1]][3]
  
  read.table(x, col.names = c("seqnames", "start", "end")) %>%
    mutate(method = method) %>%
    mutate(condition = condition) %>%
    mutate(mark = mark) %>%
    mutate(sample = paste(condition, mark, sep = "_")) %>%
    mutate(peaks = nrow(.))
})

all_peaks <- bind_rows(all_peaks)

peak_groups <- all_peaks %>% group_by(condition, mark) %>% group_keys()
all_methods <- c("gopeaks", "macs2", "seacr-relaxed", "seacr-stringent")

for (i in 1:nrow(peak_groups)) {
  
  temp_condition = as.character(peak_groups[i, "condition"])
  temp_mark = as.character(peak_groups[i, "mark"])
  
  peaks_by_condition_method <- all_peaks %>%
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark) %>%
    GRanges %>%
    GenomicRanges::split(x = ., .$method)
  
  distance_by_methods <- lapply(all_methods, function(x) {
    temp_distances <- distanceToNearest(peaks_by_condition_method[[x]])
    elementMetadata(temp_distances) %>% as.data.frame %>% mutate(method = x)
  })
  
  peak_distance_plot <- bind_rows(distance_by_methods) %>%
    ggplot(aes_string(x = "distance", fill = "method", color = "method")) +
      geom_density(alpha = 0.5) +
      scale_x_log10() +
      ggtitle(paste("Distance between closest peaks by method", temp_condition, temp_mark)) +
      xlab("Distance (bp)")
  
  peak_distance_dist <- paste0("data/figures/peak_distances/", temp_condition, "_", temp_mark, ".png")
  ggsave(peak_distance_dist, peak_distance_plot, width = 16, height = 9, dpi = 600)

  peak_distance_stat <- paste0("data/figures/peak_distances/", temp_condition, "_", temp_mark, ".txt")
  bind_rows(distance_by_methods) %>%
    group_by(method) %>%
    summarise(avg_dist = mean(distance), stdev_dist = sd(distance), med_dist = median(distance)) %>%
    mutate(condition = temp_condition) %>% mutate(mark = temp_mark) %>%
    write.table(., peak_distance_stat, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}

# Violin view
# bind_rows(distance_by_methods) %>%
#     ggplot(aes_string(x = "method", y = "distance", fill = "method", color = "method")) +
#     geom_violin(alpha = 0.5) +
#     ggtitle(paste("Distance between peaks", temp_condition, temp_mark)) +
#     scale_y_log10()