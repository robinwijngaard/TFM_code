library(readr)

#setwd("~/Dropbox/Master UOC/TFM/Genes/BED files/hg38")

setwd("~/Documents/Project/Samples/bedfiles")

ICR96_hg38_bed <- read_delim("ICR96_hg38.bed", 
                                  "\t", escape_double = FALSE, col_names = FALSE, 
                                  col_types = cols(X1 = col_character()), 
                                  trim_ws = TRUE)

ICR96_hg38_bed <- transform(ICR96_hg38_bed, x= ave(X2,X4,FUN=function(x) order(x,decreasing=F)))


# Gene_n format
ICR96_hg38_bed$X6 <- paste0(ICR96_hg38_bed$X4, "_", ICR96_hg38_bed$x)

exons_hg38 <- ICR96_hg38_bed[, c(1:4,6)]

write.table(exons_hg38, file = "exons_hg38.bed", quote=FALSE, sep="\t", row.names=FALSE)


#Gene.En format

ICR96_hg38_bed$X7 <- paste0(ICR96_hg38_bed$X4, ".E", ICR96_hg38_bed$x)

exons.hg38 <- ICR96_hg38_bed[, c(1:4,7)]

write.table(exons.hg38, file = "exons.hg38.bed", quote=FALSE, sep="\t", row.names=FALSE)
