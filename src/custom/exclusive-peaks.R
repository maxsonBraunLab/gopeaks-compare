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

if (!dir.exists("data/figures-evaluate-consensus-counts/exclusive-peaks")) {
	dir.create("data/figures-evaluate-consensus-counts/exclusive-peaks")
} else {
	unlink("data/figures-evaluate-consensus-counts/exclusive-peaks", recursive = TRUE)
	dir.create("data/figures-evaluate-consensus-counts/exclusive-peaks")
}

# Annotation per exclusive peaks ------------------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

import_peaks <- function(list_of_bed_files) {

	bind_rows(lapply(list_of_bed_files, function(x) {

		condition_mark = str_split(x, "/")[[1]][3]
		condition = str_split(condition_mark, "_")[[1]][1]
		mark = str_split(condition_mark, "_")[[1]][2]
		method = gsub(".*_(.*).bed", "\\1", x)
		data.frame(method = method, condition = condition,
					mark = mark, file = x)
	}))
}

exclusive_files <- bind_rows(
	import_peaks(list.files("data/intervene", pattern="1000_gopeaks.bed", full.names=TRUE, recursive=TRUE)),
	import_peaks(list.files("data/intervene", pattern="0100_macs2.bed", full.names=TRUE, recursive=TRUE)),
	import_peaks(list.files("data/intervene", pattern="0010_seacr-relaxed.bed", full.names=TRUE, recursive=TRUE)),
	import_peaks(list.files("data/intervene", pattern="0001_seacr-stringent.bed", full.names=TRUE, recursive=TRUE))
) %>% mutate(method_condition_mark = paste(method, condition, mark, sep = "_"))

peak_contents <- lapply(exclusive_files$file, readPeakFile)
names(peak_contents) <- exclusive_files$method_condition_mark

peak_annotation_list <- lapply(names(peak_contents), function(x) {
	# annotate peaks per method x cond x mark. if a sample only has non-canonical chr, this will error. Just return NA in that case.
	tryCatch({
		print(x)
		annotatePeak(peak_contents[[x]], TxDb = txdb, tssRegion = c(-3e3, 3e3))
	}, error = function(e) { return(NA) }
	)
})
names(peak_annotation_list) = exclusive_files$method_condition_mark
peak_annotation_list = peak_annotation_list[!unlist(lapply(peak_annotation_list, is.na))] # filter away bad annotations

peak_anno_df <- bind_rows( lapply(names(peak_annotation_list), function(x) {
  temp_sample <- peak_annotation_list[[x]]
  temp_sample@annoStat %>%
    mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
    mutate(peak_counts = temp_sample@peakNum) %>%
    mutate(meta = x) %>%
    mutate(method = gsub("(.*)_(.*)_(.*)", "\\1", meta)) %>%
    mutate(condition = gsub("(.*)_(.*)_(.*)", "\\2", meta)) %>%
    mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", meta)) %>%
    mutate(sample = gsub("(.*)_(.*)_(.*)", "\\2_\\3", meta)) %>%
    mutate(meta = NULL)
}) ) %>% arrange(condition, mark)

p <- ggplot(peak_anno_df, aes(x = method, y = Frequency, fill = Feature, group = method)) +
  geom_col(position = "stack") +
  geom_text(aes(x = method, y = 25, label = peak_counts), color = "white", size = 4, check_overlap = TRUE) +
  facet_wrap(sample ~ .) +
  ggtitle("Peak Annotation at Peak-Caller Exclusive Peaks") +
  scale_fill_brewer(palette = "Spectral") +
  theme(axis.text = element_text(color = "black"))
p

ggsave("data/figures-evaluate-consensus-counts/exclusive-peaks/exclusive_annotation.pdf", p, width = 16, height = 9, units = "in", dpi = 600)
write.table(peak_anno_df, "data/figures-evaluate-consensus-counts/exclusive-peaks/exclusive_annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)