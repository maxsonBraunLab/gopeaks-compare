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

if (!dir.exists("data/figures/peak_distances")) {
    dir.create("data/figures/peak_distances")
} else {
    unlink("data/figures/peak_distances", recursive = TRUE)
    dir.create("data/figures/peak_distances")
}

# count peaks -------------------------------------------------------------------------------------

print("Counting peaks at the replicate level")

gopeaks <- list.files("data/gopeaks", pattern = "*.bed", full.names = TRUE, recursive = TRUE)
seacr_stringent <- list.files("data/seacr", pattern = "*stringent*.bed", full.names = TRUE, recursive = TRUE)
seacr_relaxed <- list.files("data/seacr", pattern = "*relaxed*.bed", full.names = TRUE, recursive = TRUE)
macs2 <- list.files("data/macs2", pattern = "*.narrowPeak", full.names = TRUE, recursive = TRUE)

macs2_df <- data.frame(sample = str_replace(basename(macs2), "_peaks.*", ""),
                       method = "macs2",
                       file = macs2)

seacr_relaxed_df <- data.frame(sample = str_replace(basename(seacr_relaxed), ".relaxed.bed", ""),
                               method = "seacr-relaxed",
                               file = seacr_relaxed)

seacr_stringent_df <- data.frame(sample = str_replace(basename(seacr_stringent), ".stringent.bed", ""),
                               method = "seacr-stringent",
                               file = seacr_stringent)

gopeaks_df <- data.frame(sample = str_replace(basename(gopeaks), ".bed", ""),
                               method = "gopeaks",
                               file = gopeaks)

peaks_df <- bind_rows(macs2_df, seacr_relaxed_df, seacr_stringent_df, gopeaks_df)

# count peak files per sample
peak_counts <- lapply(peaks_df$file, function(x) {
    nrow(read.table(x, header = FALSE))
  }
)

peaks_df <- peaks_df %>%
              mutate(counts = unlist(peak_counts)) %>%
              mutate(file = NULL) %>%
              mutate(condition = gsub("(.*)_(.*)_(.*)", "\\1", sample)) %>%
              mutate(replicate = gsub("(.*)_(.*)_(.*)", "\\2", sample)) %>%
              mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", sample))

# define condition,mark groupings
peak_groups <- peaks_df %>% group_by(condition, mark) %>% group_keys()

# plot peak-peak distances for all condition,mark.
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
      xlab("Replicates") +
      ylab("Number of Peaks") +
      scale_y_continuous(labels = scales::comma) +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      ggtitle("Peak Counts") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))

  print(p)
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
}

# count consensus peaks ---------------------------------------------------------------------------

peak_counts <- lapply(list.files("data/consensus", pattern = "*.bed", full.names = TRUE), function(x) {
  nrow(read.table(x, header = FALSE))
})

names(peak_counts) <- lapply(list.files("data/consensus", pattern = "*.bed", full.names = TRUE), function(x) {
  str_replace(basename(x), ".bed", "")
})

peaks_df <- data.frame(metadata = names(peak_counts), counts = unname(unlist(peak_counts))) %>% 
  separate(metadata, into = c("method", "condition", "mark"), sep = "_")

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
    geom_text(aes_string(label = "counts"), vjust = 2, color = "white") +
    xlab("Method") +
    ylab("Number of Peaks") +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    ggtitle(paste(temp_condition, temp_mark, "Consensus Peak Counts")) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"))
  
  print(p)
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
  
}

# FRiP --------------------------------------------------------------------------------------------

print("Plotting FRiP")

gopeaks <- list.files("data/plotEnrichment", pattern = "gopeaks.*.tsv", full.names = TRUE, recursive = TRUE)
seacr_stringent <- list.files("data/plotEnrichment", pattern = "seacr-stringent.*.tsv", full.names = TRUE, recursive = TRUE)
seacr_relaxed <- list.files("data/plotEnrichment", pattern = "seacr-relaxed.*.tsv", full.names = TRUE, recursive = TRUE)
macs2 <- list.files("data/plotEnrichment", pattern = "macs2.*.tsv", full.names = TRUE, recursive = TRUE)

macs2_df <- data.frame(sample = str_replace(basename(macs2), "macs2_", ""),
                       method = "macs2",
                       file = macs2)
seacr_relaxed_df <- data.frame(sample = str_replace(basename(seacr_relaxed), "seacr-relaxed_", ""),
                               method = "seacr-relaxed",
                               file = seacr_relaxed)
seacr_stringent_df <- data.frame(sample = str_replace(basename(seacr_stringent), "seacr-stringent_", ""),
                                 method = "seacr_stringent",
                                 file = seacr_stringent)
gopeaks_df <- data.frame(sample = str_replace(basename(gopeaks), "gopeaks_", ""),
                         method = "gopeaks",
                         file = gopeaks)
peaks_df <- bind_rows(macs2_df, seacr_relaxed_df, seacr_stringent_df, gopeaks_df) %>%
              mutate(sample = str_replace(sample, ".tsv", ""))

# extract FRiP % per sample.
sample_FRiP <- lapply(peaks_df$file, function(x) {
    read.table(x, header = TRUE)[["percent"]]
  }
)

# add frip column, emove file extensions for samples.
peaks_df <- peaks_df %>% 
  mutate(file = NULL) %>%
  mutate(sample = str_replace(sample, ".tsv", "")) %>%
  mutate(condition = gsub("(.*)_(.*)_(.*)", "\\1", sample)) %>%
  mutate(replicate = gsub("(.*)_(.*)_(.*)", "\\2", sample)) %>%
  mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", sample)) %>%
  mutate(FRiP = unlist(sample_FRiP))

print(peaks_df)

# define condition,mark groupings
peak_groups <- peaks_df %>% group_by(condition, mark) %>% group_keys()

# plot peak-peak distances for all condition,mark.
for (i in 1:nrow(peak_groups)) {

  # define condition,mark
  temp_condition = as.character(peak_groups[i, "condition"])
  temp_mark = as.character(peak_groups[i, "mark"])

  # subset DF by condition,mark
  peaks_by_condition_method <- peaks_df %>%
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)

  # file I/O
  outfile = paste0("data/figures/FRiP/", temp_condition, "_", temp_mark, ".pdf")

  # plot peak counts
  p <- peaks_by_condition_method %>%
    ggplot(aes(x = replicate, y = FRiP, fill = method)) +
      geom_col(position = "dodge") +
      ylab("Fraction of Reads in Peaks") +
      xlab("Replicate") +
      theme(plot.background = element_rect(fill = "white")) +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))

  print(p)
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
}

# FRiP consensus peaks ----------------------------------------------------------------------------

FRiP <- lapply(list.files("data/plotEnrichment_consensus/", pattern = "*.tsv", full.names = TRUE), function(x) {
  read.table(x, header = TRUE)[["percent"]]
})

names(FRiP) <- lapply(list.files("data/plotEnrichment_consensus", pattern = "*.tsv", full.names = TRUE), function(x) {
  str_replace(basename(x), ".tsv", "")
})

FRiP_df <- data.frame(metadata = names(FRiP), FRiP = unname(unlist(FRiP))) %>% 
  separate(metadata, into = c("method", "condition", "replicate", "mark"), sep = "_") %>%
  mutate(replicate = as.factor(replicate))

peak_groups <- FRiP_df %>% group_by(condition, mark) %>% group_keys

for (i in 1:nrow(peak_groups)) {
  
  # define condition,mark
  temp_condition = as.character(peak_groups[i, "condition"])
  temp_mark = as.character(peak_groups[i, "mark"])
  
  # subset DF by condition,mark
  peaks_by_condition_mark <- FRiP_df %>%
    filter(condition == temp_condition) %>%
    filter(mark == temp_mark)
  
  # file I/O
  outfile = paste0("data/figures/FRiP/", temp_condition, "_", temp_mark, "-consensus", ".pdf")
  
  # plot FRiP
  p <- peaks_by_condition_mark %>%
    ggplot(aes_string(x = "replicate", y = "FRiP", fill = "method")) +
    geom_bar(position = "dodge", stat = 'identity') +
    geom_text(aes_string(x = "replicate", y = "FRiP", label = "FRiP"), position = position_dodge(width = 0.9), vjust = 2, color = "white") +
    ylab("Fraction of Reads in Peaks") +
    xlab("Replicate") +
    ggtitle(paste(temp_condition, temp_mark, "FRiP at Consensus Peaks")) +
    theme(plot.background = element_rect(fill = "white")) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"))
  
  print(p)
  ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 600)
}

# peak-peak distances -----------------------------------------------------------------------------

print("Calculating peak-peak distances")

# import all consensus peaks
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
      ggtitle(paste("Distance between closest peaks by method", temp_condition, temp_mark)) +
      xlab("Distance (bp)") +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      scale_color_brewer(palette = "RdYlBu", direction = -1) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white"))

  peak_distance_dist <- paste0("data/figures/peak_distances/", temp_condition, "_", temp_mark, ".pdf")
  ggsave(peak_distance_dist, peak_distance_plot, width = 8, height = 6, units = "in", dpi = 600)

  peak_distance_stat <- paste0("data/figures/peak_distances/", temp_condition, "_", temp_mark, ".txt")
  bind_rows(distance_by_methods) %>%
    group_by(method) %>%
    summarise(avg_dist = mean(distance), stdev_dist = sd(distance), med_dist = median(distance)) %>%
    mutate(condition = temp_condition) %>% mutate(mark = temp_mark) %>%
    write.table(., peak_distance_stat, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}