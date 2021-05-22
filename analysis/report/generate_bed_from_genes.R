#R > 4.0
library(biomaRt)
library(tidyr)
library(dplyr) #Load here so it does not interfere with the other select function

args = commandArgs(TRUE)
genes.file	<- 	args[1]
output.file		<-	args[2]

#Read gene and transcript identifiers
list_genes <- read.table(genes.file, header = FALSE, sep = "\t")
keys = as.character(list_genes$V2) #Keep refseq transcript ID

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

#biomaRt does not support a full query including gene and exon info. It has to be split into two queryes and merge them
#build attribute vector 1
my_attribute <- c('chromosome_name',
                  'refseq_mrna',
                  'external_gene_name',
                  'strand',
                'ensembl_transcript_id')

#fetch gene info
a <- getBM(attributes=my_attribute,
                        filters = c('refseq_mrna'),
                        values = list(refseq_mrna=keys),
                        mart = ensembl)

#build attribute vector 2
my_attribute2_exon <- c('exon_chrom_start',
                  'exon_chrom_end',
                  'rank',
                'ensembl_transcript_id')

#fetch exon info
b <- getBM(attributes=my_attribute2_exon,
                        filters = c('refseq_mrna'),
                        values = list(refseq_mrna=keys),
                        mart = ensembl)

my_refseq_loci <- merge(a,b)

#convert the strand into '-' and '+'
my_refseq_loci$strand <- gsub(pattern='-1', replacement='-', my_refseq_loci$strand)
my_refseq_loci$strand <- gsub(pattern='1', replacement='+', my_refseq_loci$strand)

#add a 'chr' into the chromosome_name
my_refseq_loci$chromosome_name <- gsub(pattern="^",
                                       replacement='chr',
                                       my_refseq_loci$chromosome_name)

#Generate the bed file with the unified columns symbol, refseq transcript and exon number as name
bed_df <- my_refseq_loci %>%
unite("name", c("external_gene_name","refseq_mrna","ensembl_transcript_id","rank"), sep=";", remove = FALSE) %>%
dplyr::select(c("chromosome_name","exon_chrom_start","exon_chrom_end", "name", "rank", "strand"))

#Write to file
write.table(bed_df, file = output.file, quote = FALSE, sep='\t', row.names = FALSE, col.names = FALSE)
