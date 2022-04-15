library(dplyr)
library(stringr)
library(tibble)
library(stringr)
library(ggplot2)

if (!dir.exists("data/figures-evaluate-consensus-counts/homer")) {
  dir.create("data/figures-evaluate-consensus-counts/homer")
}

# plot HOMER results as dot plots from RNA-Seq data -----------------------------------------------

homer_files = sort(list.files("data/homer", pattern = "knownResults.txt", recursive = TRUE, full.names = TRUE))

homer_columns = c("Motif",
                  "Consensus",
                  "pval",
                  "logpval",
                  "qval",
                  "target_seqs_with_motifs",
                  "prop_targets_with_motifs",
                  "no_bg_seqs_with_motif",
                  "prop_bg_seqs_with_motif")

homer_contents = lapply(homer_files, function(x) {
  
  info = str_split(x, "/")[[1]][3]
  method = str_split(info, "_")[[1]][1]
  condition = str_split(info, "_")[[1]][2]
  mark = str_split(info, "_")[[1]][3]
  
  
  df <- read.table(x, header = FALSE, skip = 1, col.names = homer_columns, sep = "\t") %>%
    
    # take negative log pval. convert from decimal to phred scale.
    mutate(p_pval = -log10(pval)) %>%
    
    # TF name is everything before the first "("
    mutate(name = sub("\\(.*", "", Motif)) %>% 
    
    # remove "%" sign from proportion columns and then calculate enrichment of target over bg.
    mutate(prop_targets_with_motifs = as.numeric(str_remove(prop_targets_with_motifs, "%"))) %>%
    mutate(prop_bg_seqs_with_motif = as.numeric(str_remove(prop_bg_seqs_with_motif, "%"))) %>%
    mutate(motif_enrichment = prop_targets_with_motifs / prop_bg_seqs_with_motif) %>%
    
    # paste TF name and motif together
    mutate(TF_consensus = paste(name, Consensus, sep = "-")) %>%
    
    # sample identifiers
    mutate(method = method) %>%
    mutate(condition = condition) %>%
    mutate(mark = mark)
  
})

names(homer_contents) <- lapply(homer_files, function(x) {str_split(x, "/")[[1]][3]})
all_homer <- bind_rows(homer_contents) %>% arrange(pval)

# plot TFs of interest ----------------------------------------------------------------------------

# need to fix things when pval < 10e328 because that's below R's detection limit.
all_homer = all_homer %>%
  mutate(p_pval = ifelse(p_pval == Inf, 400, p_pval)) %>%
  mutate(star = ifelse(p_pval == 400, 1, 0))

# plot sox2 TF
sox2 = all_homer %>%
  filter(mark == "Sox2") %>%
  filter(grepl("Sox", name, ignore.case = TRUE)) %>%
  ggplot(aes(x = method, y = TF_consensus, col = motif_enrichment, size = p_pval)) +
  geom_point() +
  scale_size(range = c(2, 12)) +
  xlab("") +
  ylab("") +
  theme_minimal() +
  scale_color_viridis_c() +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(labels = c("GoPeaks\n3194", "MACS2\n6147", "SEACR-Relaxed\n5531", "SEACR-Stringent\n1928"))

# add asterisks
sox2 = sox2 + 
  geom_text(label = "*", nudge_x = 0.2,
            nudge_y = 0.2, color = "black",
            size = 8, data = all_homer %>%
              filter(mark == "Sox2") %>%
              filter(grepl("Sox", name, ignore.case = TRUE)) %>%
              filter(star == 1))

ggsave("data/figures-evaluate-consensus-counts/homer/Sox2.pdf", sox2, width = 16, height = 9, dpi = 300, units = "in")

# plot CTCF TF
ctcf = all_homer %>%
  filter(mark == "CTCF") %>%
  filter(grepl("CTCF", name, ignore.case = TRUE)) %>%
  ggplot(aes(x = method, y = TF_consensus, col = motif_enrichment, size = p_pval)) +
  geom_point() +
  scale_size(range = c(2, 12)) +
  xlab("") +
  ylab("") +
  theme_minimal() +
  scale_color_viridis_c() +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(labels = c("GoPeaks\n28973", "MACS2\n42713", "SEACR-Relaxed\n19569", "SEACR-Stringent\n9572"))

# add asterisks
ctcf = ctcf + 
  geom_text(label = "*", nudge_x = 0.2,
            nudge_y = 0.2, color = "black",
            size = 8, data = all_homer %>%
              filter(mark == "CTCF") %>%
              filter(grepl("CTCF", name, ignore.case = TRUE)) %>%
              filter(star == 1))

ggsave("data/figures-evaluate-consensus-counts/homer/CTCF.pdf", ctcf, width = 16, height = 9, dpi = 300, units = "in")

