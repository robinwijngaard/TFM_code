
library(R.utils)
library(plyr)

args=commandArgs(asValues = TRUE)

library(panelcn.mops)
library(GenomicRanges)

bam_file=args$bams
bedfile=args$bed
fasta=args$fasta
exonfile=args$exon
output=args$out

#bam_file="/home/robin/Documents/Project/Samples/example/bam_ex"                                                                   #location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
#bedfile="/home/robin/Documents/Project/Samples/bedfiles/ICR96_hg38.bed"                                                                       #name of bed file
#fasta="/home/robin/Documents/Project/Samples/hg38/hg38.fa"                                                                     #name of fasta file
#exonfile="/home/robin/Documents/Project/Samples/bedfiles/exons.hg38.bed"
#output="/home/robin/Documents/Project/Results/Panelcnmops"

# Input
countWindow <- getWindows(exonfile)


# Count BamReads

bams<-list.files(bam_file,pattern=".bam",full.names=T)  
bais<-grep("bai",bams)
if(length(bais)>0){
  bams<-bams[-bais]
}

### Sample names

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



# Run Panelcnmops

XandCB <- bamCount
elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), elementMetadata(XandCB))
resultList <- runPanelcnMops(XandCB, 1:ncol(elementMetadata(bamCount)),countWindows = countWindow)

# Build results table
sampleNames <- colnames(elementMetadata(bamCount))

finalResultsTable <- createResultTable(resultlist = resultList, XandCB = XandCB, countWindows = countWindow,
                                       sampleNames = sampleNames)
allResults <- ldply(finalResultsTable, data.frame) # concat output from all samples

# Save results
colNames <- c("Sample", "Gene", "Chr", "Start", "End", "lowQual", "CN")

auxCNname <- function(x) {
  if (x %in% c("CN0", "CN1")) return("deletion") 
  else if (x %in% c("CN3", "CN4")) return("duplication")
}

outputFile <- file.path(output, "cnvFounds.txt")

filteredResults <- allResults[(allResults$CN != "CN2") & (allResults$lowQual != "lowQual"),colNames] # only deletions/duplications (CN2 means normal), and high quality
filteredResults$CNV.type <- lapply(filteredResults$CN, function(x) sapply(x, auxCNname)) # Add CNV.type column before storing file
filteredResults$CNV.type <- as.factor(unlist(filteredResults$CNV.type))  # R things...
write.table(filteredResults, outputFile, sep="\t", row.names=FALSE, quote = FALSE)  # write output file

# Save failed ROIs in a common format
failedROIs <- allResults[allResults$lowQual == "lowQual", colNames] # get all low qual
names(failedROIs)[1]<- "SampleID" # rename Sample column
failedROIs <- failedROIs[,c(1,3,4,5,2)] # reorder and filter columns
failedROIs[,1] <- unlist(strsplit(data.frame(lapply(failedROIs, as.character), stringsAsFactors=FALSE)[,1],"\\.bam"))  # remove .bam from sample names
write.table(failedROIs, file.path(output, "failedROIs.csv"), sep="\t", row.names=FALSE, quote = FALSE) # save


