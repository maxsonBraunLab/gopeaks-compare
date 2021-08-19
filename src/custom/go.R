library(dplyr)
library(stringr)
library(ChIPseeker)
library(topGO)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(data.table)
library(RColorBrewer)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

if (!dir.exists("data/GO/gopeaks_consensus_exclusive")) {
  dir.create("data/GO/gopeaks_consensus_exclusive")
} else {
  unlink("data/GO/gopeaks_consensus_exclusive", recursive = TRUE)
  dir.create("data/GO/gopeaks_consensus_exclusive")
}

if (!dir.exists("data/GO/consensus_venn")) {
  dir.create("data/GO/consensus_venn")
} else {
  unlink("data/GO/consensus_venn", recursive = TRUE)
  dir.create("data/GO/consensus_venn")
}

# functions -----------------------------------------------------------------------------
runGO <- function(geneList,xx=xx,otype,setName, Dir, printTree = 1) {
  setLength       <- sum(as.numeric(levels(geneList))[geneList]) 
  fname           <- paste(Dir, paste(setName, otype, "GO.txt", sep="_"), sep="/")
  GOData          <- new("topGOdata", ontology=otype, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
  x               <- GenTable(GOData, classicFisher=resultFisher, topNodes=length(names(resultFisher@score)))## make go table for all terms
  x               <- data.frame(x)
  pVal            <- data.frame(pval=signif(resultFisher@score, 6)) ## get unrounded pvalue
  x$enrich        <- x$Significant/x$Expected ## calculate enrichment based on what you expect by chance
  x$p.unround     <- pVal[x$GO.ID,"pval"]## put unrounded pvalue in the table
  x$p.adj         <- signif(p.adjust(x$p.unround, method="BH"), 6)## calculate the adjusted pvalue with Benjamini & Hochberg correction
  x$log.p.adj     <- -log10(x$p.adj) ## convert adjusted p value to -log10 for plot magnitude
  #x$Term.full     <- sapply(x$GO.ID, FUN=function(n){Term(xx[[n]])}) ## get the full term name
  x <- x[order(x$p.adj),]
  write.table(x, file=fname, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE) ## save the table
  ## you can print the tree if you want, but since I keep the list of all of them skip
  if(printTree>0){
    printGraph(GOData,## make the tree for the go data
               resultFisher,
               firstSigNodes = 5,
               fn.prefix = sub("_GO.txt$", "", fname),
               useInfo = "all",
               pdfSW = TRUE
    )
  }
  return(x)
}

# set up annotations --------------------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# load gene id to GO relationships as named list
geneID2GO <- get(load("anno/biomart/hg38.Ens_90.biomaRt.GO.external.geneID2GO.RData"))
xx <- get(load("anno/biomart/GO.db.Term.list.rda"))
assembly = "hg38"

# load ensembl gene names
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
if (assembly == "mm10") {
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
} else if (assembly == "hg38") {
  ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) 
} else {
  stop("Must use either mm10 or hg38 assembly")
}

gene_id_name <- getBM(attributes = c("entrezgene_id", "external_gene_name"), mart = ensembl)
gene_universe <- unique(names(geneID2GO)) # list of all possible genes

# exclusive peak gene ontology ----------------------------------------------------------

# list peak caller exclusive peaks
gopeaks_exclusive <- list.files("data/intervene", pattern = "1000_gopeaks.bed", recursive = TRUE, full.names = TRUE)
macs2_exclusive <- list.files("data/intervene", pattern = "0100_macs2.bed", recursive = TRUE, full.names = TRUE)
seacrr_exclusive <- list.files("data/intervene", pattern = "0010_seacr-relaxed.bed", recursive = TRUE, full.names = TRUE)
seacrs_exclusive <- list.files("data/intervene", pattern = "0001_seacr-stringent.bed", recursive = TRUE, full.names = TRUE)

# import peaks and annotate them with ChIPseekeR to get nearest gene for each peak.
exclusive_peak_files <- lapply(c(gopeaks_exclusive, macs2_exclusive, seacrr_exclusive), readPeakFile)
exclusive_peak_anno <- lapply(exclusive_peak_files, function(x) {
  annotatePeak(x, c(-3e3, 3e3), TxDb = txdb)
})

# GO for each exclusive peak set
for (i in names(exclusive_peak_anno)) {
  dir.create(paste0("data/GO/exclusive/", i), recursive = TRUE)
}

names(exclusive_peak_anno) <- lapply(c(gopeaks_exclusive, macs2_exclusive, seacrr_exclusive), function(x) {
  
  condition_mark <- gsub("data/intervene/(.*)/sets/(.*)_(.*).bed", replacement = "\\1", x)
  method <- gsub("data/intervene/(.*)/sets/(.*)_(.*).bed", replacement = "\\3", x)
  paste(method, condition_mark, sep = "_")
  
})

# gene ontology for exclusive peaks
go_results <- list()
for (x in names(exclusive_peak_anno)) {
  
  temp_anno <- exclusive_peak_anno[[x]]
  
  # do not run GO if there are < 10 exclusive peaks
  if (length(temp_anno@anno) < 10) {
    
    print(paste("WARNING:", x, "had < 10 exclusive peaks. topGO will not be run on this peak set."))
    go_results[[x]] <- "too short"
    
  } else {
    
    # convert entrez ID (from ChIPseekeR) to gene name
    # remove genes that have blanks as a name.
    method_exclusive_genes <- subset(gene_id_name, entrezgene_id %in% temp_anno@anno$geneId) %>%
      dplyr::select(external_gene_name) %>%
      filter(external_gene_name != "")
    method_exclusive_genes <- as.character(method_exclusive_genes$external_gene_name)
    
    # if we annotated a peak with a gene name and find it in the gene universe, set named list value to 1. Otherwise 0.
    method_exclusive_gene_list <- factor(as.integer(gene_universe %in% method_exclusive_genes))
    names(method_exclusive_gene_list) <- gene_universe
    
    # run gene ontology and export results to own directory
    temp_GO <- runGO(geneList = method_exclusive_gene_list, xx = xx,
                     otype = "BP", setName = x,
                     Dir = paste0("data/GO/exclusive/", x))
    go_results[[x]] <- temp_GO
  }
  
}

# plot exclusive results ----------------------------------------------------------------

all_exclusive_go <- list.files("data/GO/exclusive", pattern = "*.txt", recursive = TRUE, full.names = TRUE)

dir.create("data/figures/GO/exclusive", recursive = TRUE)

all_exclusive_content <- lapply(all_exclusive_go, function(x) {
  
  # define sample info 
  sample_file <- basename(x)
  method <- str_split(sample_file, "_")[[1]][1]
  condition <- str_split(sample_file, "_")[[1]][2]
  mark <- str_split(sample_file, "_")[[1]][3]
  
  # import content, add identifying cols
  fread(x, header = TRUE, sep = "\t") %>%
    mutate(method = method) %>%
    mutate(condition = condition) %>%
    mutate(mark = mark) %>%
    select(-classicFisher) %>%
    # mutate(classicFisher = as.double(classicFisher)) %>%
    filter(p.adj < 0.3)
  
})

exclusive_go <- bind_rows(all_exclusive_content)
group_keys <- exclusive_go %>% group_by(condition, mark) %>% group_keys() # group by condition,mark

for (i in 1:nrow(group_keys)) {
  
  # define condition,mark
  tmp_condition <- as.character(group_keys[i,"condition"])
  tmp_mark <- as.character(group_keys[i,"mark"])
  
  # find the top n GO.ID with the lowest pval
  top_hits <- exclusive_go %>%
    filter(condition == tmp_condition) %>%
    filter(mark == tmp_mark) %>% 
    select(GO.ID) %>% 
    head(10)
  
  # Use those GO.ID to subset the DF across methods.
  # Multiple methods can usually pick up the same term but different significance. 
  top_go_by_methods <- exclusive_go %>%
    filter(condition == tmp_condition) %>%
    filter(mark == tmp_mark) %>%
    filter(GO.ID %in% top_hits$GO.ID)
  
  # plot metadata
  tmp_title <- paste(tmp_condition, tmp_mark, "signficant GO by method")
  outfile <- paste0("data/figures/GO/exclusive/", tmp_condition, "_", tmp_mark, "-GO.pdf")
      
  # plot and export top n GO.ID across methods
  tmp_plot <- top_go_by_methods %>%
    ggplot(aes( x = method, y = Term, color = method, size = -log(p.adj) )) +
      geom_point() +
      theme_minimal() +
      scale_size(name = "-ln(p.adj)", breaks = seq(0, 20, 4), range = c(2,12)) +
      scale_color_manual(values = c("gopeaks" = "#2C7BB6", "macs2" = "#ABD9E9", "seacr-relaxed" = "#FDAE61")) +
      ggtitle(tmp_title)
  
  print(tmp_plot)
  ggsave(outfile, tmp_plot, width = 16, height = 9, units = "in")
  
}

# analyze consensus peaks ---------------------------------------------------------------

all_consensus_beds <- list.files("data/consensus/", pattern = "*.bed", full.names = TRUE)
consensus_peak_files <- lapply(all_consensus_beds, readPeakFile)
consensus_peak_anno <- lapply(consensus_peak_files, function(x) {
  annotatePeak(x, c(-3e3, 3e3), TxDb = txdb)
})

names(consensus_peak_anno) <- lapply(all_consensus_beds, function(x) {
  str_replace(basename(x), ".bed", "")  
})

for (i in names(consensus_peak_anno)) {
  dir.create(paste0("data/GO/consensus/", i), recursive = TRUE)
}

# gene ontology for consensus peaks
go_results <- list()
for (x in names(consensus_peak_anno)) {
  
  temp_anno <- consensus_peak_anno[[x]]
  
  # do not run GO if there are < 10 exclusive peaks
  if (length(temp_anno@anno) < 10) {
    
    print(paste("WARNING:", x, "had < 10 exclusive peaks. topGO will not be run on this peak set."))
    go_results[[x]] <- "too short"
    
  } else {
    
    # convert entrez ID (from ChIPseekeR) to gene name, then remove genes that have blanks as a name.
    method_consensus_genes <- subset(gene_id_name, entrezgene_id %in% temp_anno@anno$geneId) %>%
      dplyr::select(external_gene_name) %>%
      filter(external_gene_name != "")
    method_consensus_genes <- as.character(method_consensus_genes$external_gene_name)
    
    # if we annotated a peak with a gene name and find it in the gene universe, set named list value to 1. Otherwise 0.
    method_consensus_gene_list <- factor(as.integer(gene_universe %in% method_consensus_genes))
    names(method_consensus_gene_list) <- gene_universe
    
    # run gene ontology and export results to own directory
    temp_GO <- runGO(geneList = method_consensus_gene_list, xx = xx,
                     otype = "BP", setName = x,
                     Dir = paste0("data/GO/consensus/", x))
    go_results[[x]] <- temp_GO
  }
  
}

# set analysis for consensus GOs + unique GO's by method --------------------------------

all_consensus_go <- list.files("data/GO/consensus", pattern = "*.txt", recursive = TRUE, full.names = TRUE)

all_consensus_content_list <- lapply(all_consensus_go, function(x) {
  
  # define sample info 
  sample_file <- basename(x)
  method <- str_split(sample_file, "_")[[1]][1]
  condition <- str_split(sample_file, "_")[[1]][2]
  mark <- str_split(sample_file, "_")[[1]][3]
  
  # import content, add identifying cols
  fread(x, header = TRUE, sep = "\t") %>%
    mutate(method = method) %>%
    mutate(condition = condition) %>%
    mutate(mark = mark) %>%
    select(-classicFisher) %>%
    # mutate(classicFisher = as.double(classicFisher)) %>%
    filter(p.adj < 0.3)
  
})

unique_list_items <- function(x, var) {
  
  # input: x = named list. var = string of variable of interest, must be a name in x. 
  # method: take set difference of var and the rest of the other lists
  # output: list
  
  if (var %in% names(x)) {
    
    list_of_interest <- x[[var]]
    other_lists <- unlist(unname(  x[names(x)[!grepl(var, names(x))]]  )) # subset x by names != var. unname and unlist to mix all the GOs.
    unique_items_in_var <- setdiff(list_of_interest, other_lists)
    return(unique_items_in_var)

  }
}

all_methods <- c("gopeaks", "macs2", "seacr-relaxed", "seacr-stringent")

all_consensus_content %>% 
  group_by(condition, mark) %>% 
  group_walk(.keep = TRUE, .f = function(x, y) {
    
    # x = grouped dataframe
    # y = group keys in dataframe format
    
    # group GO DF by condition,mark then split DF by method
    # to get named list. key = method. value = GO terms.
    tmp_go_results_by_methods <- split(x, x$method)
    tmp_go_results_by_methods <- lapply(tmp_go_results_by_methods, function(x) {x$GO.ID})
    
    # file I/O
    outpdf=paste0("data/GO/consensus_venn/", y$condition, "_", y$mark , ".pdf")
    
    # venn diagram of GO terms split by method
    venn_title <- paste(y$condition, y$mark, "Consensus Peak GO terms by method")
    p <- ggvenn(tmp_go_results_by_methods, 
                fill_color = rev(brewer.pal(4, "RdYlBu")),
                fill_alpha = 0.75,
                stroke_size = 0,
                stroke_alpha = 0.75)
    
    pdf(outpdf, width = 8, height = 6)
    plot(p)
    dev.off()
    
    for (method in all_methods) {

      # unique GO term by method
      unique_terms <- unique_list_items(tmp_go_results_by_methods, method)

      # file I/O
      out_table <- paste0("data/GO/consensus_peak_exclusive_go/", method, "_", y$condition, "_", y$mark, ".txt")

      x %>%
        filter(GO.ID %in% unique_terms) %>%
        filter(method == method) %>%
        write.table(., out_table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    }
    
  })
