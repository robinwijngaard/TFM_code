library(biomaRt)
library(tidyr)
library(dplyr)
library(readxl)

# mama, colon, melanoma o endocri

setwd("~/Dropbox/Master UOC/TFM/Report")
panells <- c("mama", "colon", "melanoma", "endocri")

reportDir <- "~/Dropbox/Master UOC/TFM/Report"

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
bedFile <- read.delim("~/Dropbox/Master UOC/TFM/Genes/BED files/hg38/ICR96_hg38_noSNP.txt", header=FALSE)

for (panell in panells){
  genes_panell <- read_excel("gens_nm_panell.xlsx", sheet = panell, col_names = FALSE)
  colnames(genes_panell) <- c("gene", "refseq_id")
  
  my_attribute <- c('chromosome_name',
                  'refseq_mrna',
                  'external_gene_name',
                  'strand',
                  'ensembl_transcript_id')
  
  #fetch gene info
  a <- getBM(attributes=my_attribute,
           filters = c('refseq_mrna'),
           values = list(refseq_mrna=genes_panell$refseq_id),
           mart = ensembl)

  #build attribute vector 2
  my_attribute2_exon <- c('exon_chrom_start',
                        'exon_chrom_end',
                        'rank',
                        'ensembl_transcript_id',
                        'cds_start', 'cds_end')

  #fetch exon info
  b <- getBM(attributes=my_attribute2_exon,
           filters = c('refseq_mrna'),
           values = list(refseq_mrna=genes_panell$refseq_id),
           mart = ensembl)

  genes_info <- merge(a,b)
  genes_info <- genes_info[, c(2, 6, 7, 4, 3, 5, 8:10)]
  colnames(genes_info) <- c("chr", "start", "end", "gene", "nm","strand", "rank", "cds_start", "cds_end")
  write.table(genes_info, file.path(reportDir, "genes_info.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
  
  # Load bedFile
  bedFilePanel <- subset(bedFile, bedFile$V4 %in% genes_panell$gene)
  colnames(bedFilePanel) <- c("chr", "start", "end", "gene")

  # Give ID to region in bed
  bedFilePanel$ID <- 1:nrow(bedFilePanel)
  write.table(bedFilePanel, file.path(reportDir, "bedFilePanel.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
  
  setwd(reportDir)
  system(paste("bedtools intersect -wao -a bedFilePanel.bed -b genes_info.bed >", paste0(panell, "_annotated.bed")))
  system(paste("bedtools intersect -wao -a bedFilePanel.bed -b genes_info.bed -v >", paste0(panell, "_notpresent.bed")))
}

colonFile <- read.delim("colon_annotated.bed", header = FALSE)[c(1:5,10:14)]
mamaFile <- read.delim("mama_annotated.bed", header = FALSE)[c(1:5,10:14)]
melanomaFile <- read.delim("melanoma_annotated.bed", header = FALSE)[c(1:5,10:14)]
endocriFile <- read.delim("endocri_annotated.bed", header = FALSE)[c(1:5,10:14)]

newnames <- c("chr", "start", "end", "gene","id", "nm","strand", "rank", "cds_start", "cds_end")
names(colonFile) <- names(mamaFile) <- names(melanomaFile) <- names(endocriFile) <- newnames

NewExonStart <- function(file){
  file <- file[order(file$gene, file$rank), ]
  file$exonpb <- file$end - file$start
  file$intropb <- NA
  file$newstart <- NA
  file$newend <- NA
  genes <- unique(file$gene)
  for(gene in genes){
    a <- which(file$gene == gene & file$strand == 1)
    if(length(a) > 0){
      file$intropb[a] <- abs(round((file$end[a] - lead(file$start)[a]) / 1000, 0))
    }
    b <- which(file$gene == gene & file$strand == -1)
    if(length(b) > 0){
      file$intropb[b] <- round((file$start[b] - lead(file$end)[b]) / 1000, 0)
    }
    n <- which(file$gene == gene)
    newcoor <- c(0, head(cumsum(c(rbind(file$exonpb[n], file$intropb[n]))), -1))
    file$newstart[n] <- newcoor[seq(from = 1, to = length(newcoor), by = 2)]
    file$newend[n] <- newcoor[seq(from = 2, to = length(newcoor), by = 2)]
  }
  return(file)
}

colonFile <- NewExonStart(colonFile)
mamaFile <- NewExonStart(mamaFile)
endocriFile <- NewExonStart(endocriFile)
melanomaFile <- NewExonStart(melanomaFile) 

annotatedFile <- data.frame(rbind(colonFile, mamaFile, melanomaFile, endocriFile), panell = c(rep("colon", nrow(colonFile)), rep("mama", nrow(mamaFile)), rep("melanoma", nrow(melanomaFile)), rep("endocri", nrow(endocriFile)))) 
write.table(annotatedFile, file.path(reportDir, "annotatedFile.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  

