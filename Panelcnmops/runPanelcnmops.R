
library(R.utils)

args=commandArgs(asValues = TRUE)

library(panelcn.mops)
library(GenomicRanges)

#bam_file=args$bams
#bedfile=args$bed
#fasta=args$fasta
#exonfile=args$exon
#output=args$out

bam_file="/home/robin/Documents/Project/Samples/example/bam_ex"                                                                   #location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
bedfile="/home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38.bed"                                                                       #name of bed file
fasta="/home/robin/Documents/Project/Samples/hg38/hg38.fa"                                                                     #name of fasta file
exonfile="/home/robin/Documents/Project/Samples/bedfiles/exons.hg38.bed"
output="/home/robin/Documents/Project/Results/panelcnmops"

# Input
countWindow <- getWindows(exonfile)


# Count BamReads

bams<-list.files(bam_file,pattern=".bam",full.names=T)  
bais<-grep("bai",bams)
if(length(bais)>0){
  bams<-bams[-bais]
}

### Sample names
``
multi_strsplit<-function(x,splits,y){                                                  
  X<-x
  for(i in 1:length(splits)){X=strsplit(X,splits[i])[[1]][y[i]]}
  return(X)
}

a<-length(strsplit(bams[1],"/")[[1]])                                                  
sample.names<-sapply(bams,multi_strsplit,c("/",".bam"),c(a,1))                       
names(sample.names)<-NULL


bamCount <- countBamListInGRanges(countWindows = countWindow,
                                  bam.files = bams, read.width = FALSE)
