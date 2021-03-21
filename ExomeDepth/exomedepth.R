library(R.utils)

args=commandArgs(asValues = TRUE)

library(ExomeDepth)
library(GenomicRanges)


bam_file=args$bams
bedfile=args$bed
fasta=args$fasta
exonfile=args$exon
output=args$out

#bam_file="/home/robin/Documents/Project/Samples/example/bam_ex"                                                                   #location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
#bedfile="/home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38.bed"                                                                       #name of bed file
#fasta="/home/robin/Documents/Project/Samples/hg38/hg38.fa"                                                                     #name of fasta file
#exonfile="/home/robin/Documents/Project/TFM_code/Files/ExomeDepth/exons_hg38.bed"
#output="/home/robin/Documents/Project/Results/ExomeDepth"


# Create count data
bams<-list.files(bam_file,pattern=".bam",full.names=T)  
bais<-grep("bai",bams)
if(length(bais)>0){
  bams<-bams[-bais]
}

bed.file<-read.table(paste(bedfile))                                                    #reads in the bedfile and gives each column a name - expects 4 columns: chr, start, stop, name/gene.
colnames(bed.file)<-c("chromosome","start","end","name")

my.counts <- getBamCounts(bed.frame = bed.file,
                          bam.files = bams, 
                          include.chr = FALSE,
                          referenceFasta = fasta)

ExomeCount.dafr <- as(my.counts, "data.frame")


# Pipeline


#### get the annotation datasets to be used later

#data(Conrad.hg19)
#write.table(exons.hg19, "exons.hg19.txt",sep = "\t", row.names = FALSE, quote = FALSE)

exons.hg38 <- read.table(paste(exonfile))
colnames(exons.hg38)<-c("chromosome","start","end","gene", "name")

exons.hg38$chromosome <- as.character(exons.hg38$chromosome)

exons.hg38.GRanges <- GRanges(seqnames = exons.hg38$chromosome,
                              IRanges(start=exons.hg38$start,
                                      end=exons.hg38$end),
                              names = exons.hg38$name)

### prepare the main matrix of read count data
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = "X")])
nsamples <- ncol(ExomeCount.mat)

### start looping over each sample
for (i in 1:nsamples){
  
#### Create the aggregate reference set for this sample
  my.choice <- select.reference.set (test.counts =  ExomeCount.mat[,i],
                                     reference.counts = ExomeCount.mat[,-i],
                                     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                     n.bins.reduced = 10000)
  
  my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)
  
  message('Now creating the ExomeDepth object')
  all.exons <- new('ExomeDepth',test = ExomeCount.mat[,i],
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  
  
  ################ Now call the CNVs
  
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = ExomeCount.dafr$chromosome,
                        start = ExomeCount.dafr$start,
                        end = ExomeCount.dafr$end,
                        name = ExomeCount.dafr$exon)
  
  ########################### Now annotate the ExomeDepth object

  #all.exons <- AnnotateExtra(x = all.exons,
  #                           reference.annotation = Conrad.hg19.common.CNVs,
  #                           min.overlap = 0.5,
  #                           column.name = 'Conrad.hg19')
  
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.hg38.GRanges,
                             min.overlap = 0.0001,
                             column.name = 'exons.hg38')
  
  output.file <- paste('Exome_', i, '.csv', sep = '')
  
  write.csv(file = file.path(output, output.file), x = all.exons@CNV.calls, row.names = FALSE)
  
}
  
  
  
  
  
  
  
  