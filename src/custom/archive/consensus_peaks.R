library(dplyr)
library(stringr)
library(GenomicRanges)
library(plyranges)
library(tidyr)
library(data.table)

if (dir.exists("data/consensus/")) {
	unlink("data/consensus", recursive = TRUE)
	dir.create("data/consensus")
} else {
	dir.create("data/consensus")
}

macs2_files <- snakemake@input[["macs2"]]
gopeaks_files <- snakemake@input[["gopeaks"]]
seacr_relaxed_files <- snakemake@input[["seacr_relaxed"]]
seacr_stringent_files <- snakemake@input[["seacr_stringent"]]

import_peaks <- function(peak_files, method) {

	# input: list of bed files
	# method: read in first 3 columns of each file. concatenate all bed files to DF. manually define method (peak caller) column.
	# output: DF of all peaks in a method (macs2, seacr, gopeaks)

	list_of_peaks <- lapply(peak_files, function(x) {
		if (file.size(x) == 0) {
			message(paste("MESSAGE:", x, "does not have any peaks. This file will be omitted from the analysis."))
		} else {
			temp_df <- fread(x, header = FALSE, select = c(1:3), col.names = c("seqnames", "start", "end"))
			temp_df$sample <- get_sample_info(x)
			temp_df$caller <- method
			temp_df$strand <- "*"
			temp_df
		}
	})

	all_peaks_df <- bind_rows(list_of_peaks) %>%
		mutate(condition = gsub("(.*)_(.*)_(.*)", "\\1", sample)) %>%
		mutate(replicate = gsub("(.*)_(.*)_(.*)", "\\2", sample)) %>%
		mutate(mark = gsub("(.*)_(.*)_(.*)", "\\3", sample)) %>%
		select(seqnames, start, end, strand, everything()) %>%
		filter(mark != "IgG")

	return(all_peaks_df)
}

get_sample_info <- function(bed_file) {

	# input: full name of one bed file
	# method: take basename of file, rm file extension.
	# output: name of individual sample in this format: {condition}_{replicate}_{mark}

	sample_info <- str_replace(basename(bed_file), "\\..*$", "")
	# adjust for macs2 filenames
	if (str_detect(sample_info, "_peaks")) {
		sample_info <- str_replace(sample_info, "_peaks", "")
	}
	return(sample_info)
}

# all_peaks contains all peak coordinates split by method,mark,condition,replicate.
all_peaks <- bind_rows(
					import_peaks(macs2_files, method = "macs2"),
					import_peaks(gopeaks_files, method = "gopeaks"),
					import_peaks(seacr_relaxed_files, method = "seacr-relaxed"),
					import_peaks(seacr_stringent_files, method = "seacr-stringent"),
)

write.table(all_peaks, "data/consensus/all_peaks.txt",
			sep = "\t", row.names = FALSE,
			col.names = FALSE, quote = FALSE)

# list all distinct caller,mark,condition groupings
all_groups <- all_peaks %>% group_by(caller, mark, condition) %>% summarise()
print("All groupings for all_peaks:")
print(all_groups)

# write out distinct mark,condition groupings to file.
all_groups %>%
	ungroup %>%
	select(mark, condition) %>%
	distinct() %>%
	write.table(., snakemake@output[["all_groups"]],
				sep = "\t", row.names = FALSE,
				col.names = FALSE, quote = FALSE)

all_peaks <- all_peaks %>% as_granges()

consensus_peaks <- data.frame()
for (i in 1:nrow(all_groups)) {

	# loop through and define individual terms of caller,mark,condition groups
	temp_caller = as.character(all_groups[i, "caller"])
	temp_mark = as.character(all_groups[i, "mark"])
	temp_condition = as.character(all_groups[i, "condition"])

	# define how many times a peak should appear across replicates
	# to be considered a consensus peak. per condition.
	n_intersects <- snakemake@config$CONSENSUS_COUNT[[temp_condition]]
	# n_intersects <- config_file$CONSENSUS_COUNT[[temp_condition]]

	# subset all_peaks using those groupings
	subset_peaks <- subset(all_peaks, caller == temp_caller & mark == temp_mark & condition == temp_condition)

	# intersect peaks and count number of intersections by biological replicate.
	# this is an "intermediate" set of high-integrity peaks to analyze downstream.
	# the user will tell us how many replicates they want this peak to appear in.
	subset_common <- subset_peaks %>% compute_coverage() %>% filter(score >= n_intersects)

	# find overlapping peaks in the common vs subset_peaks.
	# reduce ranges to widest interval (union) at each interval.
	temp_consensus <- findOverlapPairs(subset_common, subset_peaks) %>%
		punion %>%
		reduce %>%
		as.data.frame

	# # count intersections, filter by n_intersects by condition. turn into df.
	# temp_consensus <- subset_peaks %>%
	# 	reduce() %>%
	# 	filter(score >= n_intersects) %>%
	# 	as.data.frame()

	# if consensus peaks exist, add group vars and append to consensus df.
	if (nrow(temp_consensus) == 0) {
		print(paste("MESSAGE: This group has no consensus peaks. Please double check!", paste(temp_caller, temp_mark, temp_condition)))
		next
	} else {
		temp_consensus <- temp_consensus %>%
			mutate(caller = temp_caller) %>%
			mutate(mark = temp_mark) %>%
			mutate(condition = temp_condition) %>%
			filter(start < end) # sometimes start == end. odd.

		# # filter away short consensus peaks for gopeaks.
		# if (temp_caller == "gopeaks") {
		# 	temp_consensus <- temp_consensus %>%
		# 		filter(width >= snakemake@config$CONSENSUS_MIN_WIDTH)
		# }

	}
	consensus_peaks <- rbind(consensus_peaks, temp_consensus)
}

save.image()

# print out results
for (i in 1:nrow(all_groups)) {

	temp_caller = as.character(all_groups[i, "caller"])
	temp_mark = as.character(all_groups[i, "mark"])
	temp_condition = as.character(all_groups[i, "condition"])

	subset_peaks <- subset(consensus_peaks, caller == temp_caller & mark == temp_mark & condition == temp_condition)
	subset_peaks <- subset_peaks[,1:3] %>% 
		mutate(name = paste(temp_caller, temp_mark, temp_condition, row_number(), sep = "_"))

	outfile <- paste0("data/consensus/", temp_caller, "_", temp_condition, "_", temp_mark, ".bed")
	write.table(subset_peaks, outfile, sep = "\t", row.names = FALSE, col.names = FALSE , quote = FALSE)
}