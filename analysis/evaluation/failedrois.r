# Failed ROIs
analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
bedDir <- file.path(analysisDir, "bedfiles")
evaluationDir <- file.path(analysisDir, "evaluation")
tempDir <- file.path(evaluationDir, "temp")
resultDir <- file.path(evaluationDir, "results")
graphsDir <- file.path(evaluationDir, "graphs")

for (dataset in c("all", "single")){
  csvDir <- file.path(analysisDir, "cnvfounds", dataset)
  csvFiles <- list.files(csvDir, pattern = ".csv")
  
  allFailed <- data.frame()
  
  # Read bedFile
  bedFile <- file.path(bedDir, paste0(dataset, "_rois.bed"))
  bedData <- read.table(bedFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  
  # Loop over algorithms
  for (csvFile in csvFiles){
    algorithm <- sub("failedROIs_", "", csvFile)
    algorithm <- sub(".csv", "", algorithm)
    
    # Import and adapt algorithm data (cnvFounds file)
    algorithmData <- read.table(file.path(csvDir, csvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
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
      
      # Subset bedfile for sample
      bedSample <- subset(bedData, bedData$sampleID == sample)
      bedSampleFile <- file.path(tempDir, "roisample.bed")
      write.table(bedSample, bedSampleFile, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
      
      # intersect
      setwd(tempDir)
      system(paste("bedtools intersect -wa -wb -a",sampleBed , "-b", bedSampleFile, "> Failed_included.bed"))
      
      # read files
      if (file.size("Failed_included.bed") != 0) {FailedIncluded <- read.table("Failed_included.bed", sep = "\t", stringsAsFactors=FALSE)} else {FailedIncluded <- NULL}
      
      # Edit files
      if(!is.null(FailedIncluded)){
        n <- nrow(FailedIncluded)
        algorithmID <- rep(algorithm, n)
        FailedAlgorithm <- cbind(FailedIncluded, algorithmID)
        allFailed <- rbind(allFailed, FailedAlgorithm)
      } 
    }
  }
  allFailed <- allFailed[, c(6:16)]
  colnames(allFailed) <- c(colnames(bedSample), "algorithmID")
  
  write.table(allFailed, file.path(resultDir, paste0(dataset, "_failedrois.txt")), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  
  assign(paste0(dataset, "_failedRois"), allFailed)
  
  # number of failedRois per algorithm
  n_failed <- summary(as.factor(allFailed$algorithmID))
  write.table(n_failed, file.path(resultDir, paste0(dataset, "_nfailedrois.txt")), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  
}

# Failed ROIs at gene level

comment(all_failedRois) <- "all"
comment(single_failedRois) <- "single"

for (dataset in list(all_failedRois, single_failedRois)){
  datasetName <- comment(dataset)
  
  # Read bedFile
  bedFile <- file.path(bedDir, paste0(datasetName, "_rois.bed"))
  bedData <- read.table(bedFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  
  # Extract number of ROIs per gene
  rois_gene <- data.frame(gene = c(bedData$gene))
  rois_gene_count <- rois_gene %>% group_by(gene) %>% summarise(n=n())

  # Obtain failed ROI at gene level
  dataset$gene <- factor(dataset$gene, levels = unique(dataset$gene))
  dataset$algorithmID <- factor(dataset$algorithmID, levels = unique(dataset$algorithmID))
  
  # Summary by gene
  failed_gene <- dataset %>% group_by(gene, algorithmID, .drop = FALSE) %>% summarize(n=n())
  failed_gene <- cast(failed_gene, gene ~ algorithmID)
  failed_gene <- merge(failed_gene, rois_gene_count, by = "gene")
  
  failed_gene_perc <- failed_gene
  failed_gene_perc[, 2:5] <- round(failed_gene_perc[, 2:5]/ failed_gene_perc[, 6] * 100, 4)
  
  # Export summary by gene
  write.table(failed_gene, file.path(resultDir, paste0(datasetName, "_failedgene.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  write.table(failed_gene_perc, file.path(resultDir, paste0(datasetName, "_failedgene_perc.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  

  # Summary by sample
  dataset$sampleID <- factor(dataset$sampleID)
  rois_sample <- data.frame(sampleID = c(bedData$sampleID))
  rois_sample_count <- rois_sample %>% group_by(sampleID) %>% summarise(n=n())
  
  failed_sample <- dataset %>% group_by(sampleID, algorithmID, .drop = FALSE) %>% summarize(n=n())
  failed_sample <- cast(failed_sample, sampleID ~ algorithmID)
  failed_sample <- merge(failed_sample, rois_sample_count, by = "sampleID")
  
  failed_sample_perc <- failed_sample
  failed_sample_perc[, 2:5] <- round(failed_sample_perc[, 2:5]/ failed_sample_perc[, 6] * 100, 4)
  
  # Export summary by gene
  write.table(failed_sample, file.path(resultDir, paste0(datasetName, "_failedsample.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  write.table(failed_sample_perc, file.path(resultDir, paste0(datasetName, "_failedsample_perc.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  
  # Failed ROIs included ExonCNV
  exoncnv_failed <- subset(dataset, dataset$cnv == "ExonCNV")
  exoncnv_summary <- exoncnv_failed %>%  group_by(algorithmID) %>% summarize(n=n())
  
  # Export failed exoncnv
  write.table(exoncnv_summary, file.path(resultDir, paste0(datasetName, "_exoncnvfailed.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
}
