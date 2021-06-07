# Obtain ROI files for datasets

## Comment to myself: I have corrected the EZH2 gene in validated results since it is an exon 1-20 deletion however only one exon was included in the deletion region.

# Define dir
analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
bedDir <- file.path(analysisDir, "bedfiles")
validatedDir <- file.path(bedDir, "validated")
tempDir <- file.path(bedDir, "temp")

# Load used bedfile in CNV calls
bedFile <- file.path(bedDir, "ICR96_hg38_noSNP.bed")
bedData <- read.table(bedFile, sep = "\t", stringsAsFactors=FALSE)

# Iterate over both dataset to obtain ROIs for validated results
validatedFile <- file.path(validatedDir, "ICR96_validated_regions38_all.bed")
  
# Obtain ROIs for all validated regions
roiFile <- file.path(tempDir, "all_ROI.bed")
system(paste("bedtools intersect -wb -wa -a", validatedFile, "-b", bedFile,  ">", roiFile))
roiData <- read.table(roiFile, sep = "\t", stringsAsFactors=FALSE)
  
# Fix positive ROIs per sample and gene
roiPositive <- subset(roiData, roiData$V8 == "Deletion" | roiData$V8 == "Duplication")
roiNormal <- data.frame(matrix(ncol = 5, nrow = 0))
  
# iterate over samples
samples <- unique(roiPositive$V5)
for (sample in samples){
  sampleData <- subset(roiPositive, roiPositive$V5 == sample)
  
  # iterate over genes exclude EPCAM
  genes <- sort(unique(sampleData$V4))
  genes <- genes[genes != "EPCAM"]
  for (gene in genes){
      
    # Generate specific sample-gene bedfiles for intersect function
    sampleGene <- subset(sampleData, sampleData$V4 == gene)
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

 # Generate final ROI file, with roiData + roiNormal
roiData <- roiData[, c(10, 11, 12, 13, 5, 6, 7, 8, 9)]
colnames(roiData) <- c("chr", "start", "end", "gene", "sampleID", "cnv_description", "cnv", "cnv_type", "exon_type")
  
# Add information missing columns
roiNormal$cnv_description <- rep("Normal", nrow(roiNormal))
roiNormal$cnv <- rep("Normal", nrow(roiNormal))
roiNormal$cnv_type <- rep("", nrow(roiNormal))
roiNormal$exon_type <- rep("", nrow(roiNormal))
colnames(roiNormal) <- colnames(roiData)
  
# rbind bot dataframes and export result
roiData <- rbind(roiData, roiNormal)
  
# Add roi length column
roiData$pb_length <- roiData$end - roiData$start
  
# Sort bedfile
roiData <- roiData[order(roiData$sampleID, roiData$chr, roiData$start), ]
write.table(roiData, file.path(bedDir, "all_rois.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE) 




  
   

 
  
