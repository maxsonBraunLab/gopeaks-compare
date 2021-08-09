library(dplyr)
library(stringr)
library(conflicted)
library(GenomicRanges)
library(ggplot2)
library(ChIPseeker)
library(topGO)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(doParallel)

# define standards
registerDoParallel(cores = 20)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# load GO sources
geneID2GO <- get(load("anno/biomart/hg38.Ens_90.biomaRt.GO.external.geneID2GO.RData"))
xx <- get(load("anno/biomart/GO.db.Term.list.rda"))

# load biomaRt resource
assembly = "hg38"
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
if (assembly == "mm10") {
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
} else if (assembly == "hg38") {
  ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) 
} else {
  stop("Must use either mm10 or hg38 assembly")
}

# get entrez gene number and name; define gene universe
gene_id_name <- getBM(attributes = c("entrezgene_id", "external_gene_name"), mart = ensembl)
gene_universe <- unique(names(geneID2GO))

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

# method-exclusive peaks --------------------------------------------------------------------------
gopeaks_exclusive <- list.files("data/intervene", pattern = "1000_gopeaks.bed", recursive = TRUE, full.names = TRUE)
macs2_exclusive <- list.files("data/intervene", pattern = "0100_macs2.bed", recursive = TRUE, full.names = TRUE)
seacrr_exclusive <- list.files("data/intervene", pattern = "0010_seacr-relaxed.bed", recursive = TRUE, full.names = TRUE)
seacrs_exclusive <- list.files("data/intervene", pattern = "0001_seacr-stringent.bed", recursive = TRUE, full.names = TRUE)

exclusive_peak_files <- lapply(c(gopeaks_exclusive, macs2_exclusive, seacrr_exclusive), readPeakFile)
exclusive_peak_anno <- lapply(exclusive_peak_files, function(x) {
  return(annotatePeak(x, c(-3e3, 3e3), TxDb = txdb))
})

names(exclusive_peak_anno) <- lapply(c(gopeaks_exclusive, macs2_exclusive, seacrr_exclusive), function(x) {

  condition_mark <- gsub("data/intervene/(.*)/sets/(.*)_(.*).bed", replacement = "\\1", x)
  method <- gsub("data/intervene/(.*)/sets/(.*)_(.*).bed", replacement = "\\3", x)
  paste(method, condition_mark, sep = "_")

})

# convert entrez ID to gene name, set important genes to 1 and bg genes to 0, and run GO.
exclusive_peak_go <- list()
# for (x in names(exclusive_peak_anno)) {
foreach(x = names(exclusive_peak_anno)) %dopar% {

  temp_anno <- exclusive_peak_anno[[x]]
  
  if (length(temp_anno@anno) < 10) {
    print(paste("WARNING:", x, "had < 10 exclusive peaks. topGO will not be run on this peak set."))
    exclusive_peak_go[[x]] <- "too short"
  } else {
    dir.create(paste0("data/GO/exclusive/", x), recursive = TRUE)
    method_exclusive_genes <- subset(gene_id_name, entrezgene_id %in% temp_anno@anno$geneId) %>%
      select(external_gene_name) %>%
      filter(external_gene_name != "") # remove genes that have blanks as a name.
    method_exclusive_genes <- as.character(method_exclusive_genes$external_gene_name)
    
    method_exclusive_gene_list <- factor(as.integer(gene_universe %in% method_exclusive_genes))
    names(method_exclusive_gene_list) <- gene_universe

    print(x)
    print(length(method_exclusive_gene_list[method_exclusive_gene_list == 1]))
    temp_GO <- runGO(geneList = method_exclusive_gene_list, xx = xx, otype = "BP", setName = x, Dir = paste0("data/GO/exclusive/", x))
    print(head(temp_GO))
    exclusive_peak_go[[x]] <- temp_GO
  }
}

# consensus peaks ---------------------------------------------------------------------------------
consensus_files <- list.files("data/consensus", pattern = "*.bed", full.names=TRUE)
consensus_peak_files <- lapply(consensus_files, readPeakFile)
consensus_peak_anno <- lapply(consensus_peak_files, function(x) {
  return(annotatePeak(x, c(-3e3, 3e3), TxDb = txdb))
})

names(consensus_peak_anno) <- lapply(consensus_files, function(x) {
  str_replace(basename(x), ".bed", "")
})

consensus_peak_go <- list()
foreach(x = names(consensus_peak_anno)) %dopar% {
# for (x in names(consensus_peak_anno)) {
  
  temp_anno <- consensus_peak_anno[[x]]
  
  if (length(temp_anno@anno) < 10) {
    print(paste("WARNING:", x, "had < 10 exclusive peaks. topGO will not be run on this peak set."))
    consensus_peak_go[[x]] <- "too short"
  } else {
    dir.create(paste0("data/GO/consensus/", x), recursive = TRUE)
    method_exclusive_genes <- subset(gene_id_name, entrezgene_id %in% temp_anno@anno$geneId) %>%
      select(external_gene_name) %>%
      filter(external_gene_name != "") # remove genes that have blanks as a name.
    method_exclusive_genes <- as.character(method_exclusive_genes$external_gene_name)
    
    method_exclusive_gene_list <- factor(as.integer(gene_universe %in% method_exclusive_genes))
    names(method_exclusive_gene_list) <- gene_universe

    print(x)
    print(length(method_exclusive_gene_list[method_exclusive_gene_list == 1]))
    temp_GO <- runGO(geneList = method_exclusive_gene_list, xx = xx, otype = "BP", setName = x, Dir = paste0("data/GO/consensus/", x))
    print(head(temp_GO))
    consensus_peak_go[[x]] <- temp_GO
  }
}

save.image()