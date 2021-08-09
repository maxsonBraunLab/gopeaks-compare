library(dplyr)
library(stringr)
library(ChIPseeker)
library(topGO)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(data.table)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

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
