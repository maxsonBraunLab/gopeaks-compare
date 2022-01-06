library(dplyr)
library(stringr)
library(conflicted)
library(GenomicRanges)
library(ggplot2)
library(tidyr)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

if (!dir.exists("data/figures/peak-counts")) {
    dir.create("data/figures/peak-counts")
} else {
    unlink("data/figures/peak-counts", recursive = TRUE)
    dir.create("data/figures/peak-counts")
}

if (!dir.exists("data/figures/FRiP")) {
    dir.create("data/figures/FRiP")
} else {
    unlink("data/figures/FRiP", recursive = TRUE)
    dir.create("data/figures/FRiP")
}

if (!dir.exists("data/figures/peak-distances")) {
    dir.create("data/figures/peak-distances")
} else {
    unlink("data/figures/peak-distances", recursive = TRUE)
    dir.create("data/figures/peak-distances")
}

method_finder = function(filename) {
  if (str_detect(filename, "gopeaks")) {
    return("gopeaks")
  } else if (str_detect(filename, "macs2")) {
    return("macs2")
  } else if (str_detect(filename, "relaxed")) {
    return("seacr-relaxed")
  } else if (str_detect(filename, "stringent")) {
    return("seacr-stringent")
  }
}

# count peaks -------------------------------------------------------------------------------------

print("Counting peaks at the replicate level")

gopeaks <- list.files("data/gopeaks", pattern = "*.bed", full.names = TRUE, recursive = TRUE)
macs2 <- list.files("data/macs2", pattern = "Peak", full.names = TRUE, recursive = TRUE)
seacr_relaxed <- list.files("data/seacr", pattern = "*relaxed*.bed", full.names = TRUE, recursive = TRUE)
seacr_stringent <- list.files("data/seacr", pattern = "*stringent*.bed", full.names = TRUE, recursive = TRUE)

peaks_df = bind_rows(lapply(c(gopeaks, macs2, seacr_relaxed, seacr_stringent), function(x) {
  # file i/o
  sample = strsplit(basename(x), "\\.")[[1]][1]
  condition = strsplit(sample, "_")[[1]][1]
  replicate = strsplit(sample, "_")[[1]][2]
  mark = strsplit(sample, "_")[[1]][3]
  peaks = nrow(read.table(x)) # assume each line is a peak.
  method = method_finder(x)

  data.frame(method = method,
             sample = sample,
             condition = condition,
             replicate = replicate,
             mark = mark,
             counts = peaks)
})) %>% mutate(sample = str_replace(sample, "_peaks", "")) # fix macs2 _peaks suffix

# define condition,mark groupings
peak_groups <- peaks_df %>% group_by(condition, mark) %>% group_keys()

# plot peak counts for all condition,mark.
for (i in 1:nrow(peak_groups)) {

  # define condition,mark
  temp_condition = as.character(peak_groups[i, "condition"])
  temp_mark = as.character(peak_groups[i, "mark"])

  # subset DF by condition,mark
  peaks_by_condition_method <- peaks_df %>%
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)

  # file I/O
  outfile = paste0("data/figures/peak-counts/", temp_condition, "_", temp_mark, ".pdf")

  # plot peak counts
  p <- peaks_by_condition_method %>%
    ggplot(aes_string(x = "replicate", y = "counts", fill = "method")) +
      geom_col(position = "dodge") +
      geom_text(aes_string(label = "counts"), vjust = 2, color = "white", check_overlap = TRUE) +
      xlab("Replicates") +
      ylab("Number of Peaks") +
      scale_y_continuous(labels = scales::comma) +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      ggtitle("Peak Counts") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
}

# count consensus peaks ---------------------------------------------------------------------------

print("Counting peaks at the consensus level")

peaks_df = bind_rows(lapply(list.files("data/consensus", pattern = "*.bed", full.names=TRUE), function(x) {
  # file i/o
  sample = strsplit(basename(x), "\\.")[[1]][1]
  method = strsplit(sample, "_")[[1]][1]
  condition = strsplit(sample, "_")[[1]][2]
  mark = strsplit(sample, "_")[[1]][3]
  peaks = nrow(read.table(x)) # assume each line is a peak.

  data.frame(method = method,
             sample = sample,
             condition = condition,
             mark = mark,
             counts = peaks)
}))

peak_groups <- peaks_df %>% group_by(condition, mark) %>% group_keys

for (i in 1:nrow(peak_groups)) {
  
  # define condition,mark
  temp_condition = as.character(peak_groups[i, "condition"])
  temp_mark = as.character(peak_groups[i, "mark"])
  
  # subset DF by condition,mark
  peaks_by_condition_method <- peaks_df %>%
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)
  
  # file I/O
  outfile = paste0("data/figures/peak-counts/", temp_condition, "_", temp_mark, "-consensus", ".pdf")
  
  # plot peak counts
  p <- peaks_by_condition_method %>%
    ggplot(aes_string(x = "method", y = "counts", fill = "method")) +
    geom_col(position = "dodge") +
    geom_text(aes_string(label = "counts"), vjust = 2, color = "white", check_overlap = TRUE) +
    xlab("Method") +
    ylab("Number of Peaks") +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    ggtitle(paste(temp_condition, temp_mark, "Consensus Peak Counts")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"))
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
}

# FRiP --------------------------------------------------------------------------------------------

print("Plotting FRiP at the replicate level")

frip_df = bind_rows(lapply(list.files("data/plotEnrichment", pattern="*.tsv", full.names=TRUE), function(x) {
  # file i/o
  sample = strsplit(basename(x), "\\.")[[1]][1]
  method = strsplit(sample, "_")[[1]][1]
  condition = strsplit(sample, "_")[[1]][2]
  replicate = strsplit(sample, "_")[[1]][3]
  mark = strsplit(sample, "_")[[1]][4]
  frip = read.table(x, header = TRUE)[1, "percent"]

  data.frame(method = method,
             sample = sample,
             condition = condition,
             replicate = replicate,
             mark = mark,
             frip = frip)
}))

# define condition,mark groupings
peak_groups <- frip_df %>% group_by(condition, mark) %>% group_keys()

# plot peak-peak distances for all condition,mark.
for (i in 1:nrow(peak_groups)) {

  # define condition,mark
  temp_condition = as.character(peak_groups[i, "condition"])
  temp_mark = as.character(peak_groups[i, "mark"])

  # subset DF by condition,mark
  peaks_by_condition_method <- frip_df %>%
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)

  # file I/O
  outfile = paste0("data/figures/FRiP/", temp_condition, "_", temp_mark, ".pdf")

  # plot peak counts
  p <- peaks_by_condition_method %>%
    ggplot(aes(x = replicate, y = frip, fill = method)) +
      geom_col(position = "dodge") +
      ylab("Fraction of Reads in Peaks") +
      xlab("Replicate") +
      theme(plot.background = element_rect(fill = "white")) +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
}

# FRiP consensus peaks ----------------------------------------------------------------------------

print("Plotting FRiP at the consensus level")

frip_df = bind_rows(lapply(list.files("data/plotEnrichment_consensus", pattern="*.tsv", full.names=TRUE), function(x) {
  # file i/o
  sample = strsplit(basename(x), "\\.")[[1]][1]
  method = strsplit(sample, "_")[[1]][1]
  condition = strsplit(sample, "_")[[1]][2]
  replicate = strsplit(sample, "_")[[1]][3]
  mark = strsplit(sample, "_")[[1]][4]
  frip = read.table(x, header = TRUE)[1, "percent"]

  data.frame(method = method,
             sample = sample,
             condition = condition,
             replicate = replicate,
             mark = mark,
             frip = frip)
}))

peak_groups <- frip_df %>% group_by(condition, mark) %>% group_keys

for (i in 1:nrow(peak_groups)) {
  
  # define condition,mark
  temp_condition = as.character(peak_groups[i, "condition"])
  temp_mark = as.character(peak_groups[i, "mark"])
  
  # subset DF by condition,mark
  peaks_by_condition_mark <- frip_df %>%
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)
  
  # file I/O
  outfile = paste0("data/figures/FRiP/", temp_condition, "_", temp_mark, "-consensus", ".pdf")
  
  # plot FRiP
  p <- peaks_by_condition_mark %>%
    ggplot(aes_string(x = "replicate", y = "frip", fill = "method")) +
    geom_bar(position = "dodge", stat = 'identity') +
    geom_text(aes_string(x = "replicate", y = "frip", label = "frip"), position = position_dodge(width = 0.9), vjust = 2, color = "white") +
    ylab("Fraction of Reads in Peaks") +
    xlab("Replicate") +
    ggtitle(paste(temp_condition, temp_mark, "FRiP at Consensus Peaks")) +
    theme(plot.background = element_rect(fill = "white")) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"))
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
}

# peak-peak distances -----------------------------------------------------------------------------

print("Plotting peak-peak distances")

# import all consensus peaks
all_peaks <- bind_rows(lapply(list.files("data/consensus", pattern = "*.bed", full.names = TRUE), function(x) {
  file = str_replace(basename(x), ".bed", "")
  method = str_split(file, "_")[[1]][1]
  condition = str_split(file, "_")[[1]][2]
  mark = str_split(file, "_")[[1]][3]

  read.table(x, col.names = c("seqnames", "start", "end")) %>%
    mutate(method = method) %>%
    mutate(condition = condition) %>%
    mutate(mark = mark) %>%
    mutate(sample = paste(condition, mark, sep = "_"))
}))

# define condition,mark groupings
peak_groups <- all_peaks %>% group_by(condition, mark) %>% group_keys()
all_methods <- c("gopeaks", "macs2", "seacr-relaxed", "seacr-stringent")

# plot peak-peak distances for all condition,mark.
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
    elementMetadata(temp_distances) %>%
        as.data.frame %>%
        mutate(method = x)
  })

  peak_distance_plot <- bind_rows(distance_by_methods) %>%
    ggplot(aes_string(x = "distance", fill = "method", color = "method")) +
      geom_density(alpha = 0.5) +
      scale_x_log10() +
      ggtitle(paste(temp_condition, temp_mark, "peak-peak distances")) +
      xlab("Distance (bp)") +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      scale_color_brewer(palette = "RdYlBu", direction = -1) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))

  peak_distance_dist <- paste0("data/figures/peak-distances/", temp_condition, "_", temp_mark, ".pdf")
  ggsave(peak_distance_dist, peak_distance_plot, width = 8, height = 6, units = "in", dpi = 600)

  peak_distance_stat <- paste0("data/figures/peak-distances/", temp_condition, "_", temp_mark, ".txt")
  bind_rows(distance_by_methods) %>%
    group_by(method) %>%
    summarise(avg_dist = mean(distance), stdev_dist = sd(distance), med_dist = median(distance)) %>%
    mutate(condition = temp_condition) %>% mutate(mark = temp_mark) %>%
    write.table(., peak_distance_stat, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}