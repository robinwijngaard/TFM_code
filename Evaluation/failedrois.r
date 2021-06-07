# Failed ROIs
library(R.utils)
library(dplyr)
library(reshape)

args=commandArgs(asValues = TRUE)

# Define dir
outputDir <- args$outputDir
failedRois <- args$failedRois
roisFile <- args$roisFile

# make temp dir, results dir en graphs dir
tempDir <- file.path(outputDir, "temp")
resultDir <- file.path(outputDir, "results")
graphsDir <- file.path(outputDir, "graphs")

if(!dir.exists(tempDir)) dir.create(tempDir)
if(!dir.exists(resultDir)) dir.create(resultDir)
if(!dir.exists(graphsDir)) dir.create(graphsDir)

csvFiles <- list.files(failedRois, pattern = ".csv")
Failed_Rois <- data.frame()
  
# Read bedFile
roisData <- read.table(roisFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  
# Loop over algorithms
for (csvFile in csvFiles){
  algorithm <- sub("failedROIs_", "", csvFile)
  algorithm <- sub(".csv", "", algorithm)
    
  # Import and adapt algorithm data (failedRois file)
  algorithmData <- read.table(file.path(failedRois, csvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  algorithmData <- algorithmData[, c(2:5, 1)]
  algorithmData$SampleID <- sub("X", "", algorithmData$SampleID)
    
  # Loop over samples
  samples <- sort(unique(algorithmData$SampleID))
  for (sample in samples){
      
    # Prepare sample failedRois
    s <- which(algorithmData$SampleID == sample)
    sampleData <- algorithmData[s, ]
    sampleBed <- file.path(tempDir, "failed.bed")
    write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
      
    # Subset roisfile for sample
    roisSample <- subset(roisData, roisData$sampleID == sample)
    roisSampleFile <- file.path(tempDir, "roisample.bed")
    write.table(roisSample, roisSampleFile, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
      
    # intersect
    setwd(tempDir)
    system(paste("bedtools intersect -wa -wb -a",sampleBed , "-b", roisSampleFile, "> Failed_included.bed"))
      
    # read files
    if (file.size("Failed_included.bed") != 0) {FailedIncluded <- read.table("Failed_included.bed", sep = "\t", stringsAsFactors=FALSE)} else {FailedIncluded <- NULL}
      
    # Edit files
    if(!is.null(FailedIncluded)){
      n <- nrow(FailedIncluded)
      algorithmID <- rep(algorithm, n)
      FailedAlgorithm <- cbind(FailedIncluded, algorithmID)
      Failed_Rois <- rbind(Failed_Rois, FailedAlgorithm)
    }
  }
}

Failed_Rois <- Failed_Rois[, c(6:16)]
colnames(Failed_Rois) <- c(colnames(roisSample), "algorithmID")
  
write.table(Failed_Rois, file.path(resultDir, "Failed_ROIS.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  
  
# number of failedRois per algorithm
n_failed <- summary(as.factor(Failed_Rois$algorithmID))
write.table(n_failed, file.path(resultDir, "nfailedrois.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  

# Extract number of ROIs per gene
rois_gene <- data.frame(gene = c(roisData$gene))
rois_gene_count <- rois_gene %>% group_by(gene) %>% summarise(n=n())

# Obtain failed ROI at gene level
Failed_Rois$gene <- factor(Failed_Rois$gene, levels = unique(Failed_Rois$gene))
Failed_Rois$algorithmID <- factor(Failed_Rois$algorithmID, levels = unique(Failed_Rois$algorithmID))
  
# Summary by gene
failed_gene <- Failed_Rois %>% group_by(gene, algorithmID, .drop = FALSE) %>% summarize(n=n())
failed_gene <- cast(failed_gene, gene ~ algorithmID)
failed_gene <- merge(failed_gene, rois_gene_count, by = "gene")
  
failed_gene_perc <- failed_gene
failed_gene_perc[, 2:5] <- round(failed_gene_perc[, 2:5]/ failed_gene_perc[, 6] * 100, 4)
  
# Export summary by gene
write.table(failed_gene, file.path(resultDir, "failedgene.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(failed_gene_perc, file.path(resultDir, "failedgene_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  

# Summary by sample
Failed_Rois$sampleID <- factor(Failed_Rois$sampleID)
rois_sample <- data.frame(sampleID = c(roisData$sampleID))
rois_sample_count <- rois_sample %>% group_by(sampleID) %>% summarise(n=n())
  
failed_sample <- Failed_Rois %>% group_by(sampleID, algorithmID, .drop = FALSE) %>% summarize(n=n())
failed_sample <- cast(failed_sample, sampleID ~ algorithmID)
failed_sample <- merge(failed_sample, rois_sample_count, by = "sampleID")
  
failed_sample_perc <- failed_sample
failed_sample_perc[, 2:5] <- round(failed_sample_perc[, 2:5]/ failed_sample_perc[, 6] * 100, 4)
  
# Export summary by gene
write.table(failed_sample, file.path(resultDir, "failedsample.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(failed_sample_perc, file.path(resultDir,"failedsample_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  
# Failed ROIs included ExonCNV
exoncnv_failed <- subset(Failed_Rois, Failed_Rois$cnv == "ExonCNV")
exoncnv_summary <- exoncnv_failed %>%  group_by(algorithmID) %>% summarize(n=n())
  
# Export failed exoncnv
write.table(exoncnv_summary, file.path(resultDir, "exoncnvfailed.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  

