# Obtain validated ROIs per sample from validated file
library(R.utils)
args=commandArgs(asValues = TRUE)

outputDir <- args$outputDir
bedFile <- args$bedFile
validatedFile <- args$validatedFile
EPCAM <- args$EPCAM

# make temp dir, results dir en graphs dir
tempDir <- file.path(outputDir, "temp")
resultDir <- file.path(outputDir, "results")
graphsDir <- file.path(outputDir, "graphs")

if(!dir.exists(tempDir)) {dir.create(tempDir)}
if(!dir.exists(resultDir)) {dir.create(resultDir)}
if(!dir.exists(graphsDir)) {dir.create(graphsDir)}

# Obtain ROIs for all validated regions
roiFile <- file.path(tempDir, "sample_ROI.bed")
system(paste("bedtools intersect -wb -wa -a", validatedFile, "-b", bedFile,  ">", roiFile))
roiData <- read.table(roiFile, sep = "\t", stringsAsFactors=FALSE)

# Order columsn
roiData <- roiData[, c(10:13, 5:9)]
names(roiData) <- c("chr", "start", "end", "gene", "sampleID", "cnv_description", "cnv", "cnv_type", "exon_type")

# Read bedData
bedData <- read.table(bedFile, sep = "\t", stringsAsFactors=FALSE)

# Fix positives
roiPositive <- subset(roiData, roiData$cnv == "ExonCNV")
roiNormal <- data.frame(matrix(ncol = 5, nrow = 0))

# iterate over samples
samples <- unique(roiPositive$sampleID)
for (sample in samples){
  sampleData <- subset(roiPositive, roiPositive$sampleID == sample)
  
  # iterate over genes
  genes <- sort(unique(sampleData$gene))
  if(EPCAM) {genes[genes !="EPCAM"]}
  for (gene in genes){
    
    # Generate specific sample-gene bedfiles for intersect function
    sampleGene <- subset(sampleData, sampleData$gene == gene)
    bedGene <- subset(bedData, bedData$V4 == gene)
    
    sampleGeneFile <- file.path(tempDir, "samplegene.bed")
    write.table(sampleGene, sampleGeneFile, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
    
    bedGeneFile <- file.path(tempDir, "bedGene.bed")
    write.table(bedGene, bedGeneFile, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # Intersect no overlapping to obtain NORMAL ROIs from the gene with CNV
    roigeneFile <- file.path(tempDir, "roinormal.bed")
    system(paste("bedtools intersect -wa -a", bedGeneFile, "-b", sampleGeneFile , "-v >", roigeneFile))
    
    if (file.size(roigeneFile) != 0) {roigenenormal <- read.table(roigeneFile, sep = "\t", stringsAsFactors=FALSE)} else {roigenenormal <- NULL}
    
    # Add sample ID and add to roiNormal dataframe
    if(!is.null(roigenenormal)){
      n <- nrow(roigenenormal)
      sampleID <- rep(sample, n)
      roigenenormal <- cbind(roigenenormal, sampleID)
      roiNormal <- rbind(roiNormal, roigenenormal)
    }
  }
}

# Add colnames and missing columns
names(roiNormal) <- c("chr", "start", "end", "gene", "sampleID")

# Add missing columns
roiNormal$cnv_description <- rep("Normal", nrow(roiNormal))
roiNormal$cnv <- rep("Normal", nrow(roiNormal))
roiNormal$cnv_type <- rep("", nrow(roiNormal))
roiNormal$exon_type <- rep("", nrow(roiNormal))

# Rbind
roiData <- rbind(roiData, roiNormal)

# Add roi length column
roiData$pb_length <- roiData$end - roiData$start

# Sort bedfile
roiData <- roiData[order(roiData$sampleID, roiData$chr, roiData$start), ]
write.table(roiData, file.path(outputDir, "validated_rois.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE) 