

# Failed ROIs

## allDataset

csvFiles <- list.files(allDir, pattern = ".csv")

allFailed <- data.frame(matrix(ncol = 11, nrow = 0))

for (csvFile in csvFiles){
  algorithm <- sub("failedROIs_", "", csvFile)
  algorithm <- sub(".csv", "", algorithm)
  
  # Import and adapt algorithm data (cnvFounds file)
  algorithmData <- read.table(file.path(allDir, csvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  algorithmData <- algorithmData[, c(2:5, 1)]
  algorithmData$SampleID <- sub("X", "", algorithmData$SampleID)
  
  samples <- sort(unique(algorithmData$SampleID))
  
  for (sample in samples){
    s <- which(algorithmData$SampleID == sample)
    sampleData <- algorithmData[s, ]
    sampleBed <- file.path(tempDir, paste0(algorithm,"_failed_", sample, ".bed"))
    
    write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # Prepare sample bed file, filter by sample
    alldsSample <- subset(allDatasetBed, allDatasetBed$V4 == sample)
    alldsSampleBed <- file.path(tempDir, "allsample.bed")
    write.table(alldsSample, alldsSampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
    
    # intersect
    setwd(tempDir)
    system(paste("bedtools intersect -wa -wb -a",sampleBed , "-b", alldsSampleBed, "> Failed_included.bed"))
    
    # read files
    if (file.size("Failed_included.bed") != 0) {FailedIncluded <- read.table("Failed_included.bed", sep = "\t", stringsAsFactors=FALSE)} else {FailedIncluded <- NULL}
    
    # Edit files
    if(!is.null(FailedIncluded)){
      n <- nrow(FailedIncluded)
      algorithmID <- rep(algorithm, n)
      FailedAlgorithm <- cbind(algorithmID, FailedIncluded)
      allFailed <- rbind(allFailed, FailedAlgorithm)
    } 
  }
}






# Failed ROIs

rois_gene <- data.frame(genes = c(allData[, 4], allNormal[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FailedROIs_df <- allFailed[, c(1, 5, 6)]
FailedROIs_df$V4 <- factor(FailedROIs_df$V4, levels = unique(FailedROIs_df$V4))
FailedROIs_df$algorithmID <- factor(FailedROIs_df$algorithmID)

FailedROIs_summary <- FailedROIs_df %>% group_by(V4, algorithmID, .drop = FALSE) %>% summarize(n=n())
FailedROIs_summary <- cast(FailedROIs_summary, V4 ~ algorithmID)

names(FailedROIs_summary)[1] <- "genes"
FailedROIs_summary <- merge(FailedROIs_summary, rois_gene_count, by = "genes")

FailedROIs_summary_perc <- FailedROIs_summary
FailedROIs_summary_perc[, 2:5] <- round(FailedROIs_summary_perc[, 2:5]/ FailedROIs_summary_perc[, 6] * 100, 4)

write.table(FailedROIs_summary, file.path(resultDir, "FailedROIS.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FailedROIs_summary_perc, file.path(resultDir, "FailedROIS_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  



FailedROIs_df$V5 <- factor(FailedROIs_df$V5, levels = unique(FailedROIs_df$V5))

FailedROIs_samples <- FailedROIs_df %>% group_by(V5, algorithmID, .drop = FALSE) %>% summarize(n=n())
FailedROIs_samples <- cast(FailedROIs_samples, V5 ~ algorithmID)

write.table(FailedROIs_samples, file.path(resultDir, "FailedROIs_samples.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  










all_failedRois <- summary(as.factor(allFailed$algorithmID))
single_failedRois <- summary(as.factor(singleFailed$algorithmID))
clinic_failedRois <- summary(as.factor(clinicFailed$algorithmID))
