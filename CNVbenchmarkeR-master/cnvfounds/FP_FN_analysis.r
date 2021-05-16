
# Define dirs
analysisDir <- "/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/cnvfounds"
allDir <- file.path(analysisDir, "all")
singleDir <- file.path(analysisDir, "single")
clinicDir <- file.path(analysisDir, "clinic")
tempDir <- file.path(analysisDir, "temp")

# load bed files positive samples
allFile <- file.path(analysisDir, "all_positives.bed")
allData <- read.table(allFile, sep = "\t", stringsAsFactors=FALSE)
allNormalFile <- file.path(analysisDir, "all_normal.bed")
allNormal <- read.table(allNormalFile, sep = "\t", stringsAsFactors=FALSE)

clinicFile <- file.path(analysisDir, "clinic_positives.bed")
clinicData <- read.table(clinicFile, sep = "\t", stringsAsFactors=FALSE)
clinicNormalFile <- file.path(analysisDir, "clinic_normal.bed")
clinicNormal <- read.table(clinicNormalFile, sep = "\t", stringsAsFactors=FALSE)

singleFile <- file.path(analysisDir, "single_positives.bed")
singleData <- read.table(singleFile, sep = "\t", stringsAsFactors=FALSE)
singleNormalFile <- file.path(analysisDir, "single_normal.bed")
singleNormal <- read.table(singleNormalFile, sep = "\t", stringsAsFactors=FALSE)

# Transform bedfiles
allData <- allData[, c(10, 11, 12, 13, 4, 6, 8, 9)]
clinicData <- clinicData[, c(10, 11, 12, 13, 4, 6, 8, 9)]
singleData <-  singleData[, c(10, 11, 12, 13, 4, 6, 8, 9)]

allNormal <- allNormal[, c(8, 9, 10, 11, 4, 6)]
clinicNormal <- clinicNormal[, c(8, 9, 10, 11, 4, 6)]
singleNormal <- singleNormal[, c(8, 9, 10, 11, 4, 6)]

# Add identifier
allData$ID <- 1:nrow(allData)
allNormal$ID <- 1:nrow(allNormal)
allNormal$ID <- allNormal$ID + nrow(allData)

clinicData$ID <- 1:nrow(clinicData)
clinicNormal$ID <- 1:nrow(clinicNormal)
clinicNormal$ID <- clinicNormal$ID + nrow(clinicData)

singleData$ID <- 1:nrow(singleData)
singleNormal$ID <- 1:nrow(singleNormal)
singleNormal$ID <- singleNormal$ID + nrow(singleData)


# Make bedfile with all ROIS per dataset (normal + positives)
allDatasetBed <- rbind(allData[, 1:5], setNames(allNormal[, 1:5], names(allData)[1:5]))
clinicDatasetBed <- rbind(clinicData[, 1:5], setNames(clinicNormal[, 1:5], names(clinicData)[1:5]))
singleDatasetBed <- rbind(singleData[, 1:5], setNames(singleNormal[, 1:5], names(singleData)[1:5]))

write.table(allDatasetBed, file.path(tempDir, "allDataset.bed"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(clinicDatasetBed, file.path(tempDir, "clinicDataset.bed"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(singleDatasetBed, file.path(tempDir, "singleDataset.bed"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


TP_all <- data.frame(matrix(ncol = 15, nrow = 0))
FP_all <- data.frame(matrix(ncol = 13, nrow = 0))
FN_all <- data.frame(matrix(ncol = 11, nrow = 0))

# All
cnvFiles <- list.files(allDir, pattern = ".txt")
for (cnvFile in cnvFiles){
  algorithm <- sub("cnvFounds_", "", cnvFile)
  algorithm <- sub(".txt", "", algorithm)
  
  # Import and adapt algorithm data (cnvFounds file)
  algorithmData <- read.table(file.path(allDir, cnvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  names(algorithmData) <- tolower(names(algorithmData))
  names(algorithmData) <- sub("chromosome", "chr", names(algorithmData))
  names(algorithmData) <- sub("cnv.type", "type", names(algorithmData))
  colNames <- c("sample", "chr", "start", "end", "type")
  algorithmData <- algorithmData[, colNames] 
  algorithmData$sample <- as.character(algorithmData$sample)
  algorithmData$sample <- sub("X", "", algorithmData$sample)
  
  # For all
  samples <- unique(algorithmData$sample)
  
  for (sample in samples){
    s <- which(algorithmData$sample == sample)
    sampleData <- algorithmData[s, ]
    sampleData <- sampleData[, -1]
    sampleBed <- file.path(tempDir, paste0(algorithm,"_", sample, ".bed"))
    write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # Prepare sample bed file, filter by sample
    allSample <- subset(allData, allData$V4 == sample)
    allBed <- file.path(tempDir, "all.bed")
    write.table(allSample, allBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
    
    allNormalsample <- subset(allNormal, allNormal$V4== sample)
    allNormalbed <- file.path(tempDir, "allNormal.bed")
    write.table(allNormalsample, allNormalbed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # intersect
    setwd(tempDir)
    system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", allBed, "> TP.bed"))
    system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", allNormalbed, "> FP.bed"))
    system(paste("bedtools intersect -wa -wb -a", allBed, "-b", sampleBed, "-v > FN.bed"))

    # read files
    if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
    if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}
    if (file.size("FN.bed") != 0) {FNsample <- read.table("FN.bed", sep = "\t", stringsAsFactors=FALSE)} else {FNsample <- NULL}
    
    
    # Edit files
    if(!is.null(TPsample)){
      n <- nrow(TPsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      TPsample <- cbind(sampleID, algorithmID, TPsample)
      TP_all <- rbind(TP_all, TPsample)
    }
    
    if(!is.null(FPsample)){
      n <- nrow(FPsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      FPsample <- cbind(sampleID, algorithmID, FPsample)
      FP_all <- rbind(FP_all, FPsample)
    }
    
    if(!is.null(FNsample)){
      n <- nrow(FNsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      FNsample <- cbind(sampleID, algorithmID, FNsample)
      FN_all <- rbind(FN_all, FNsample)
    }
  }
}


TP_single <- data.frame(matrix(ncol = 15, nrow = 0))
FP_single <- data.frame(matrix(ncol = 13, nrow = 0))
FN_single <- data.frame(matrix(ncol = 11, nrow = 0))

# Single
cnvFiles <- list.files(singleDir, pattern = ".txt")
for (cnvFile in cnvFiles){
  algorithm <- sub("cnvFounds_", "", cnvFile)
  algorithm <- sub(".txt", "", algorithm)
  
  # Import and adapt algorithm data (cnvFounds file)
  algorithmData <- read.table(file.path(singleDir, cnvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  names(algorithmData) <- tolower(names(algorithmData))
  names(algorithmData) <- sub("chromosome", "chr", names(algorithmData))
  names(algorithmData) <- sub("cnv.type", "type", names(algorithmData))
  colNames <- c("sample", "chr", "start", "end", "type")
  algorithmData <- algorithmData[, colNames] 
  algorithmData$sample <- as.character(algorithmData$sample)
  algorithmData$sample <- sub("X", "", algorithmData$sample)
  
  # For all
  samples <- unique(algorithmData$sample)
  
  for (sample in samples){
    s <- which(algorithmData$sample == sample)
    sampleData <- algorithmData[s, ]
    sampleData <- sampleData[, -1]
    sampleBed <- file.path(tempDir, paste0(algorithm,"_", sample, ".bed"))
    write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # Prepare sample bed file, filter by sample
    singleSample <- subset(singleData, singleData$V4 == sample)
    singleBed <- file.path(tempDir, "single.bed")
    write.table(singleSample, singleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
    
    singleNormalsample <- subset(singleNormal, singleNormal$V4== sample)
    singleNormalbed <- file.path(tempDir, "singleNormal.bed")
    write.table(singleNormalsample, singleNormalbed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # intersect
    setwd(tempDir)
    system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", singleBed, "> TP.bed"))
    system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", singleNormalbed, "> FP.bed"))
    system(paste("bedtools intersect -wa -wb -a", singleBed, "-b", sampleBed, "-v > FN.bed"))

    # read files
    if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
    if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}
    if (file.size("FN.bed") != 0) {FNsample <- read.table("FN.bed", sep = "\t", stringsAsFactors=FALSE)} else {FNsample <- NULL}
    
    
    # Edit files
    if(!is.null(TPsample)){
      n <- nrow(TPsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      TPsample <- cbind(sampleID, algorithmID, TPsample)
      TP_single <- rbind(TP_single, TPsample)
    }
    
    if(!is.null(FPsample)){
      n <- nrow(FPsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      FPsample <- cbind(sampleID, algorithmID, FPsample)
      FP_single <- rbind(FP_single, FPsample)
    }
    
    if(!is.null(FNsample)){
      n <- nrow(FNsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      FNsample <- cbind(sampleID, algorithmID, FNsample)
      FN_single <- rbind(FN_single, FNsample)
    }
  }
}


TP_clinic <- data.frame(matrix(ncol = 15, nrow = 0))
FP_clinic <- data.frame(matrix(ncol = 13, nrow = 0))
FN_clinic <- data.frame(matrix(ncol = 11, nrow = 0))


cnvFiles <- list.files(clinicDir, pattern = ".txt")
for (cnvFile in cnvFiles){
  algorithm <- sub("cnvFounds_", "", cnvFile)
  algorithm <- sub(".txt", "", algorithm)
  
  # Import and adapt algorithm data (cnvFounds file)
  algorithmData <- read.table(file.path(clinicDir, cnvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  names(algorithmData) <- tolower(names(algorithmData))
  names(algorithmData) <- sub("chromosome", "chr", names(algorithmData))
  names(algorithmData) <- sub("cnv.type", "type", names(algorithmData))
  colNames <- c("sample", "chr", "start", "end", "type")
  algorithmData <- algorithmData[, colNames] 
  algorithmData$sample <- as.character(algorithmData$sample)
  algorithmData$sample <- sub("X", "", algorithmData$sample)
  
  # For all
  samples <- unique(algorithmData$sample)
  
  for (sample in samples){
    s <- which(algorithmData$sample == sample)
    sampleData <- algorithmData[s, ]
    sampleData <- sampleData[, -1]
    sampleBed <- file.path(tempDir, paste0(algorithm,"_", sample, ".bed"))
    write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # Prepare sample bed file, filter by sample
    clinicSample <- subset(clinicData, clinicData$V4 == sample)
    clinicBed <- file.path(tempDir, "clinic.bed")
    write.table(clinicSample, clinicBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
    
    clinicNormalsample <- subset(clinicNormal, clinicNormal$V4== sample)
    clinicNormalbed <- file.path(tempDir, "clinicNormal.bed")
    write.table(clinicNormalsample, clinicNormalbed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # intersect
    setwd(tempDir)
    system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", clinicBed, "> TP.bed"))
    system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", clinicNormalbed, "> FP.bed"))
    system(paste("bedtools intersect -wa -wb -a", clinicBed, "-b", sampleBed, "-v > FN.bed"))
    
    # read files
    if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
    if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}
    if (file.size("FN.bed") != 0) {FNsample <- read.table("FN.bed", sep = "\t", stringsAsFactors=FALSE)} else {FNsample <- NULL}
    
    # Edit files
    if(!is.null(TPsample)){
      n <- nrow(TPsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      TPsample <- cbind(sampleID, algorithmID, TPsample)
      TP_clinic <- rbind(TP_clinic, TPsample)
    }
    
    if(!is.null(FPsample)){
      n <- nrow(FPsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      FPsample <- cbind(sampleID, algorithmID, FPsample)
      FP_clinic <- rbind(FP_clinic, FPsample)
    }
    
    if(!is.null(FNsample)){
      n <- nrow(FNsample)
      sampleID <- rep(sample, n)
      algorithmID <- rep(algorithm, n)
      FNsample <- cbind(sampleID, algorithmID, FNsample)
      FN_clinic <- rbind(FN_clinic, FNsample)
    }
  }
}


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

## singleDataset

csvFiles <- list.files(singleDir, pattern = ".csv")

singleFailed <- data.frame(matrix(ncol = 11, nrow = 0))

for (csvFile in csvFiles){
  algorithm <- sub("failedROIs_", "", csvFile)
  algorithm <- sub(".csv", "", algorithm)
  
  # Import and adapt algorithm data (cnvFounds file)
  algorithmData <- read.table(file.path(singleDir, csvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  algorithmData <- algorithmData[, c(2:5, 1)]
  algorithmData$SampleID <- sub("X", "", algorithmData$SampleID)
  
  samples <- sort(unique(algorithmData$SampleID))
  
  for (sample in samples){
    s <- which(algorithmData$SampleID == sample)
    sampleData <- algorithmData[s, ]
    sampleBed <- file.path(tempDir, paste0(algorithm,"_failed_", sample, ".bed"))
    
    write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # Prepare sample bed file, filter by sample
    singledsSample <- subset(singleDatasetBed, singleDatasetBed$V4 == sample)
    singledsSampleBed <- file.path(tempDir, "singlesample.bed")
    write.table(singledsSample, singledsSampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
    
    # intersect
    setwd(tempDir)
    system(paste("bedtools intersect -wa -wb -a",sampleBed , "-b", singledsSampleBed, "> Failed_included.bed"))
    
    # read files
    if (file.size("Failed_included.bed") != 0) {FailedIncluded <- read.table("Failed_included.bed", sep = "\t", stringsAsFactors=FALSE)} else {FailedIncluded <- NULL}
    
    # Edit files
    if(!is.null(FailedIncluded)){
      n <- nrow(FailedIncluded)
      algorithmID <- rep(algorithm, n)
      FailedAlgorithm <- cbind(algorithmID, FailedIncluded)
      singleFailed <- rbind(singleFailed, FailedAlgorithm)
    } 
  }
}

## Clinic

csvFiles <- list.files(clinicDir, pattern = ".csv")

clinicFailed <- data.frame(matrix(ncol = 11, nrow = 0))

for (csvFile in csvFiles){
  algorithm <- sub("failedROIs_", "", csvFile)
  algorithm <- sub(".csv", "", algorithm)
  
  # Import and adapt algorithm data (cnvFounds file)
  algorithmData <- read.table(file.path(clinicDir, csvFile), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  algorithmData <- algorithmData[, c(2:5, 1)]
  algorithmData$SampleID <- sub("X", "", algorithmData$SampleID)
  
  samples <- sort(unique(algorithmData$SampleID))
  
  for (sample in samples){
    s <- which(algorithmData$SampleID == sample)
    sampleData <- algorithmData[s, ]
    sampleBed <- file.path(tempDir, paste0(algorithm,"_failed_", sample, ".bed"))
    
    write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
    
    # Prepare sample bed file, filter by sample
    clinicdsSample <- subset(clinicDatasetBed, clinicDatasetBed$V4 == sample)
    clinicdsSampleBed <- file.path(tempDir, "clinicsample.bed")
    write.table(clinicdsSample, clinicdsSampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
    
    # intersect
    setwd(tempDir)
    system(paste("bedtools intersect -wa -wb -a",sampleBed , "-b", clinicdsSampleBed, "> Failed_included.bed"))
    
    # read files
    if (file.size("Failed_included.bed") != 0) {FailedIncluded <- read.table("Failed_included.bed", sep = "\t", stringsAsFactors=FALSE)} else {FailedIncluded <- NULL}
    
    # Edit files
    if(!is.null(FailedIncluded)){
      n <- nrow(FailedIncluded)
      algorithmID <- rep(algorithm, n)
      FailedAlgorithm <- cbind(algorithmID, FailedIncluded)
      clinicFailed <- rbind(clinicFailed, FailedAlgorithm)
    } 
  }
}


# Compare del i dup
TP_all$V11 <- tolower(TP_all$V11)
TP_clinic$V11 <- tolower(TP_clinic$V11)
TP_single$V11 <- tolower(TP_single$V11)

a <- which(TP_all$V4 != TP_all$V11)
TP_all[a, ]

b <- which(TP_clinic$V4 != TP_clinic$V11)
TP_clinic[b, ]

c <- which(TP_single$V4 != TP_single$V11)
TP_single[c, ]

# FInal data frame
allData$cnvkit5 <- 0
allData$convading <- 0
allData$decon <- 0
allData$exomedepth <- 0
allData$manta <- 0
allData$panelcn <- 0

algorithms <- unique(TP_all$algorithmID)

for(i in 1:nrow(allData)){
  for(j in algorithms){
    a <- which(TP_all$algorithmID == j & TP_all$V13 == i)
    if(length(a) > 0){
      allData[i, j] <- 1
    }
  }
}

# FInal data frame Single
singleData$cnvkit5 <- 0
singleData$convading <- 0
singleData$decon <- 0
singleData$exomedepth <- 0
singleData$manta <- 0
singleData$panelcn <- 0

algorithms <- unique(TP_all$algorithmID)

for(i in 1:nrow(singleData)){
  for(j in algorithms){
    a <- which(TP_single$algorithmID == j & TP_single$V13 == i)
    if(length(a) > 0){
      singleData[i, j] <- 1
    }
  }
}

# FInal data frame clinic
clinicData$cnvkit5 <- 0
clinicData$convading <- 0
clinicData$decon <- 0
clinicData$exomedepth <- 0
clinicData$manta <- 0
clinicData$panelcn <- 0

algorithms <- unique(TP_all$algorithmID)

for(i in 1:nrow(clinicData)){
  for(j in algorithms){
    a <- which(TP_clinic$algorithmID == j & TP_clinic$V13 == i)
    if(length(a) > 0){
      clinicData[i, j] <- 1
    }
  }
}

resultDir <- file.path(analysisDir, "results")

write.table(allData, file.path(resultDir, "allData.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicData, file.path(resultDir, "clinicData.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleData, file.path(resultDir, "singleData.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  


colSums(allData[, 10:15])


# False positives

# FInal data frame
allNormal$cnvkit5 <- 0
allNormal$convading <- 0
allNormal$decon <- 0
allNormal$exomedepth <- 0
allNormal$manta <- 0
allNormal$panelcn <- 0

algorithms <- unique(FP_all$algorithmID)

for(i in 1:nrow(allNormal)){
  for(j in algorithms){
    a <- which(FP_all$algorithmID == j & FP_all$V11 == i)
    if(length(a) > 0){
      allNormal[i, j] <- 1
    }
  }
}

# FInal data frame Single
singleNormal$cnvkit5 <- 0
singleNormal$convading <- 0
singleNormal$decon <- 0
singleNormal$exomedepth <- 0
singleNormal$manta <- 0
singleNormal$panelcn <- 0

algorithms <- unique(FP_single$algorithmID)

for(i in 1:nrow(singleNormal)){
  for(j in algorithms){
    a <- which(FP_single$algorithmID == j & FP_single$V11 == i)
    if(length(a) > 0){
      singleNormal[i, j] <- 1
    }
  }
}

# FInal data frame clinic
clinicNormal$cnvkit5 <- 0
clinicNormal$convading <- 0
clinicNormal$decon <- 0
clinicNormal$exomedepth <- 0
clinicNormal$manta <- 0
clinicNormal$panelcn <- 0

algorithms <- unique(FP_clinic$algorithmID)

for(i in 1:nrow(clinicNormal)){
  for(j in algorithms){
    a <- which(FP_clinic$algorithmID == j & FP_clinic$V11 == i)
    if(length(a) > 0){
      clinicNormal[i, j] <- 1
    }
  }
}

resultDir <- file.path(analysisDir, "results")

write.table(allNormal, file.path(resultDir, "allNormal.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicNormal, file.path(resultDir, "clinicNormal.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleNormal, file.path(resultDir, "singleNormal.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  




write.table(FN_all, file.path(resultDir, "FN_all.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(FN_single, file.path(resultDir, "FN_single.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(FN_clinic, file.path(resultDir, "FN_clinic.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  






# Result table
TP <- as.vector(colSums(allData[, 10:15]))
FN <- as.vector(nrow(allData) - as.vector(colSums(allData[, 10:15])))
FP <- as.vector(colSums(allNormal[, 8:13]))
TN <- as.vector(nrow(allNormal) - as.vector(colSums(allNormal[, 8:13])))

Total <- TP + FN + FP + TN
Spec <- round(TN / (TN + FP), 4)
Sens <- round(TP / (TP + FN), 4)

allDataset <- cbind(TP, TN, FP, FN, Total, Sens, Spec, names(allData[, 10:15]))

TP <- as.vector(colSums(singleData[, 10:15]))
FN <- as.vector(nrow(singleData) - as.vector(colSums(singleData[, 10:15])))
FP <- as.vector(colSums(singleNormal[, 8:13]))
TN <- as.vector(nrow(singleNormal) - as.vector(colSums(singleNormal[, 8:13])))

Total <- TP + FN + FP + TN
Spec <- round(TN / (TN + FP), 4)
Sens <- round(TP / (TP + FN), 4)

singleDataset <- cbind(TP, TN, FP, FN, Total, Sens, Spec, names(allData[, 10:15]))

TP <- as.vector(colSums(clinicData[, 10:15]))
FN <- as.vector(nrow(clinicData) - as.vector(colSums(clinicData[, 10:15])))
FP <- as.vector(colSums(clinicNormal[, 8:13]))
TN <- as.vector(nrow(clinicNormal) - as.vector(colSums(clinicNormal[, 8:13])))

Total <- TP + FN + FP + TN
Spec <- round(TN / (TN + FP), 4)
Sens <- round(TP / (TP + FN), 4)

clinicDataset <- cbind(TP, TN, FP, FN, Total, Sens, Spec, names(allData[, 10:15]))

all_failedRois <- summary(as.factor(allFailed$algorithmID))
single_failedRois <- summary(as.factor(singleFailed$algorithmID))
clinic_failedRois <- summary(as.factor(clinicFailed$algorithmID))

write.table(allDataset, file.path(resultDir, "allDataset_results.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(clinicDataset, file.path(resultDir, "clinicDataset_results.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(singleDataset, file.path(resultDir, "singleDataset_results.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  



# VennDiagram
setwd(resultDir)
library("VennDiagram")

## Alldata

## All calls
all_venn <- rbind(allData[, 9:15], allNormal[, 7:13])
all_venn$ID <- as.character(all_venn$ID)
all_venn <- all_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(all_venn)[i]
  all_venn_algo <- all_venn[, c(1, i)]
  all_venn_algo <- subset(all_venn_algo, all_venn_algo[,2] != 0)
  venn_list[[algo]] <- all_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_allcalls.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TP
tp_venn <- rbind(allData[, 9:15])
tp_venn$ID <- as.character(tp_venn$ID)
tp_venn <- tp_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(tp_venn)[i]
  tp_venn_algo <- tp_venn[, c(1, i)]
  tp_venn_algo <- subset(tp_venn_algo, tp_venn_algo[,2] != 0)
  venn_list[[algo]] <- tp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tp.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TN
tn_venn <- rbind(allNormal[, 7:13])
tn_venn$ID <- as.character(tn_venn$ID)
tn_venn <- tn_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(tn_venn)[i]
  tn_venn_algo <- tn_venn[, c(1, i)]
  tn_venn_algo <- subset(tn_venn_algo, tn_venn_algo[,2] == 0)
  venn_list[[algo]] <- tn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FN
fn_venn <- rbind(allData[, 9:15])
fn_venn$ID <- as.character(fn_venn$ID)
fn_venn <- fn_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(fn_venn)[i]
  fn_venn_algo <- fn_venn[, c(1, i)]
  fn_venn_algo <- subset(fn_venn_algo, fn_venn_algo[,2] == 0)
  venn_list[[algo]] <- fn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FP
fp_venn <- rbind(allNormal[, 7:13])
fp_venn$ID <- as.character(fp_venn$ID)
fp_venn <- fp_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(fp_venn)[i]
  fp_venn_algo <- fp_venn[, c(1, i)]
  fp_venn_algo <- subset(fp_venn_algo, fp_venn_algo[,2] != 0)
  venn_list[[algo]] <- fp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fp.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## Single dataset

## All calls
all_venn <- rbind(singleData[, 9:15], singleNormal[, 7:13])
all_venn$ID <- as.character(all_venn$ID)
all_venn <- all_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(all_venn)[i]
  all_venn_algo <- all_venn[, c(1, i)]
  all_venn_algo <- subset(all_venn_algo, all_venn_algo[,2] != 0)
  venn_list[[algo]] <- all_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_allcalls_single.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TP
tp_venn <- rbind(singleData[, 9:15])
tp_venn$ID <- as.character(tp_venn$ID)
tp_venn <- tp_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(tp_venn)[i]
  tp_venn_algo <- tp_venn[, c(1, i)]
  tp_venn_algo <- subset(tp_venn_algo, tp_venn_algo[,2] != 0)
  venn_list[[algo]] <- tp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tp_single.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TN
tn_venn <- rbind(singleNormal[, 7:13])
tn_venn$ID <- as.character(tn_venn$ID)
tn_venn <- tn_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(tn_venn)[i]
  tn_venn_algo <- tn_venn[, c(1, i)]
  tn_venn_algo <- subset(tn_venn_algo, tn_venn_algo[,2] == 0)
  venn_list[[algo]] <- tn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tn_single.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FN
fn_venn <- rbind(singleData[, 9:15])
fn_venn$ID <- as.character(fn_venn$ID)
fn_venn <- fn_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(fn_venn)[i]
  fn_venn_algo <- fn_venn[, c(1, i)]
  fn_venn_algo <- subset(fn_venn_algo, fn_venn_algo[,2] == 0)
  venn_list[[algo]] <- fn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fn_single.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FP
fp_venn <- rbind(singleNormal[, 7:13])
fp_venn$ID <- as.character(fp_venn$ID)
fp_venn <- fp_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(fp_venn)[i]
  fp_venn_algo <- fp_venn[, c(1, i)]
  fp_venn_algo <- subset(fp_venn_algo, fp_venn_algo[,2] != 0)
  venn_list[[algo]] <- fp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fp_single.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)


## Clinic dataset

## All calls
all_venn <- rbind(clinicData[, 9:15], clinicNormal[, 7:13])
all_venn$ID <- as.character(all_venn$ID)
all_venn <- all_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(all_venn)[i]
  all_venn_algo <- all_venn[, c(1, i)]
  all_venn_algo <- subset(all_venn_algo, all_venn_algo[,2] != 0)
  venn_list[[algo]] <- all_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_allcalls_clinic.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TP
tp_venn <- rbind(clinicData[, 9:15])
tp_venn$ID <- as.character(tp_venn$ID)
tp_venn <- tp_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(tp_venn)[i]
  tp_venn_algo <- tp_venn[, c(1, i)]
  tp_venn_algo <- subset(tp_venn_algo, tp_venn_algo[,2] != 0)
  venn_list[[algo]] <- tp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tp_clinic.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TN
tn_venn <- rbind(clinicNormal[, 7:13])
tn_venn$ID <- as.character(tn_venn$ID)
tn_venn <- tn_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(tn_venn)[i]
  tn_venn_algo <- tn_venn[, c(1, i)]
  tn_venn_algo <- subset(tn_venn_algo, tn_venn_algo[,2] == 0)
  venn_list[[algo]] <- tn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tn_clinic.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FN
fn_venn <- rbind(clinicData[, 9:15])
fn_venn$ID <- as.character(fn_venn$ID)
fn_venn <- fn_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(fn_venn)[i]
  fn_venn_algo <- fn_venn[, c(1, i)]
  fn_venn_algo <- subset(fn_venn_algo, fn_venn_algo[,2] == 0)
  venn_list[[algo]] <- fn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fn_clinic.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FP
fp_venn <- rbind(clinicNormal[, 7:13])
fp_venn$ID <- as.character(fp_venn$ID)
fp_venn <- fp_venn[, -6]

venn_list <- list()

for(i in 2:6){
  algo <- names(fp_venn)[i]
  fp_venn_algo <- fp_venn[, c(1, i)]
  fp_venn_algo <- subset(fp_venn_algo, fp_venn_algo[,2] != 0)
  venn_list[[algo]] <- fp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fp_clinic.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)




# MDS ## 750x750
## All
library("MASS")

all_mds <- rbind(allData[, 9:15], allNormal[, 7:13])
all_mds <- all_mds[, c(-1,-6)]

d <- dist(t(all_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(13, 13, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

## Positives
TP_mds <- rbind(allData[, 9:15])
TP_mds <- TP_mds[, c(-1,-6)]

d <- dist(t(TP_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(2, 2.2, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

##Normal
TN_mds <- rbind(allNormal[, 7:13])
TN_mds <- TN_mds[, c(-1,-6)]

d <- dist(t(TN_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(13, 12, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

## Single
library("MASS")
single_mds <- rbind(singleData[, 9:15], singleNormal[, 7:13])
single_mds <- single_mds[, c(-1,-6)]

d <- dist(t(single_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(6, 12, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

## Positives
TP_mds <- rbind(singleData[, 9:15])
TP_mds <- TP_mds[, c(-1,-6)]

d <- dist(t(TP_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(1.5, 2, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

##Normal
TN_mds <- rbind(singleNormal[, 7:13])
TN_mds <- TN_mds[, c(-1,-6)]

d <- dist(t(TN_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(13, 12, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

## Clinic
library("MASS")
clinic_mds <- rbind(clinicData[, 9:15], clinicNormal[, 7:13])
clinic_mds <- clinic_mds[, c(-1,-6)]

d <- dist(t(clinic_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(13, 9, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

## Positives
TP_mds <- rbind(clinicData[, 9:15])
TP_mds <- TP_mds[, c(-1,-6)]

d <- dist(t(TP_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(1.5, 2, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

##Normal
TN_mds <- rbind(clinicNormal[, 7:13])
TN_mds <- TN_mds[, c(-1,-6)]

d <- dist(t(TN_mds), method = "euclidean")
MDS <- isoMDS(d)

plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5)
labels <- names(all_mds)
legend(13, 12, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)






# Corplot
library(corrplot)
corr.all <- rbind(allData[, 10:15], allNormal[, 8:13])
corr.all <- corr.all[, -5]

cor.data <- cor(corr.all)
corrplot(cor.data, method = "number", cl.lim = c(0,1))


corr.single <- rbind(singleData[, 10:15], singleNormal[, 8:13])
corr.single <- corr.single[, -5]

cor.data <- cor(corr.single)
corrplot(cor.data, method = "number", cl.lim = c(0,1))


corr.clinic <- rbind(clinicData[, 10:15], clinicNormal[, 8:13])
corr.clinic <- corr.clinic[, -5]

cor.data <- cor(corr.clinic)
corrplot(cor.data, method = "number", cl.lim = c(0,1))






# Cluster
clust <- rbind(allData[, 9:15], allNormal[, 7:13])
clust <- clust[, c(-1,-6)]

d <- dist(t(clust), method = "euclidean")
all.hc <- hclust(d, method = "ward.D2")
plot(all.hc)


clust <- rbind(singleData[, 9:15], singleNormal[, 7:13])
clust <- clust[, c(-1,-6)]

d <- dist(t(clust), method = "euclidean")
all.hc <- hclust(d, method = "ward.D2")
plot(all.hc)


clust <- rbind(clinicData[, 9:15], clinicNormal[, 7:13])
clust <- clust[, c(-1,-6)]

d <- dist(t(clust), method = "euclidean")
all.hc <- hclust(d, method = "ward.D2")
plot(all.hc)


# Model hibdrid

allData$hibrid <- 0
which(sum(allData[, 9:13] >= 3) & (allData$exomedepth == 1 | allData$decon == 1))
which(sum(allData[, 9:13] >= 3) & (allData$exomedepth == 1 | allData$decon == 1))

read.table(file.path())
ICR96bed <- read.table(file.path(analysisDir, "bedfiles", "ICR96_hg38_noSNP_"), sep = "\t", stringsAsFactors=FALSE)













# ROIS per gene
library(dplyr)
library(reshape)

## FP

rois_gene <- data.frame(genes = c(allNormal[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FP_all$V8 <- factor(FP_all$V8,  levels = rois_gene_count$genes)
FP_all$algorithmID <- factor(FP_all$algorithmID, levels = algorithms)

FP_analysis <- FP_all %>% group_by(V8, algorithmID, .drop = FALSE) %>% summarize(n=n())
FP_analysis <- cast(FP_analysis, V8 ~ algorithmID)
names(FP_analysis)[1] <- "genes"

FP_analysis <- merge(FP_analysis, rois_gene_count, by = "genes")
FP_analysis_perc <- FP_analysis
FP_analysis_perc[, 2:7] <- round(FP_analysis_perc[, 2:7]/ FP_analysis_perc[, 8] * 100, 4)

write.table(FP_analysis, file.path(resultDir, "FP_analysis_all.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FP_analysis_perc, file.path(resultDir, "FP_analysis_all_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  



rois_gene <- data.frame(genes = c(singleNormal[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FP_single$V8 <- factor(FP_single$V8,  levels = rois_gene_count$genes)
FP_single$algorithmID <- factor(FP_single$algorithmID, levels = algorithms)

FP_analysis <- FP_single %>% group_by(V8, algorithmID, .drop = FALSE) %>% summarize(n=n())
FP_analysis <- cast(FP_analysis, V8 ~ algorithmID)
names(FP_analysis)[1] <- "genes"

FP_analysis <- merge(FP_analysis, rois_gene_count, by = "genes")
FP_analysis_perc <- FP_analysis
FP_analysis_perc[, 2:7] <- round(FP_analysis_perc[, 2:7]/ FP_analysis_perc[, 8] * 100, 4)

write.table(FP_analysis, file.path(resultDir, "FP_analysis_single.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FP_analysis_perc, file.path(resultDir, "FP_analysis_single_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  




rois_gene <- data.frame(genes = c(clinicNormal[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FP_clinic$V8 <- factor(FP_clinic$V8,  levels = rois_gene_count$genes)
FP_clinic$algorithmID <- factor(FP_clinic$algorithmID, levels = algorithms)

FP_analysis <- FP_clinic %>% group_by(V8, algorithmID, .drop = FALSE) %>% summarize(n=n())
FP_analysis <- cast(FP_analysis, V8 ~ algorithmID)
names(FP_analysis)[1] <- "genes"

FP_analysis <- merge(FP_analysis, rois_gene_count, by = "genes")
FP_analysis_perc <- FP_analysis
FP_analysis_perc[, 2:7] <- round(FP_analysis_perc[, 2:7]/ FP_analysis_perc[, 8] * 100, 4)

write.table(FP_analysis, file.path(resultDir, "FP_analysis_clinic.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FP_analysis_perc, file.path(resultDir, "FP_analysis_clinic_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  



## FN
rois_gene <- data.frame(genes = c(allData[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FN_all$V4 <- factor(FN_all$V4,  levels = rois_gene_count$genes)
FN_all$algorithmID <- factor(FN_all$algorithmID, levels = algorithms)

FN_analysis <- FN_all %>% group_by(V4, algorithmID, .drop = FALSE) %>% summarize(n=n())
FN_analysis <- cast(FN_analysis, V4 ~ algorithmID)
names(FN_analysis)[1] <- "genes"

FN_analysis <- merge(FN_analysis, rois_gene_count, by = "genes")
FN_analysis_perc <- FN_analysis
FN_analysis_perc[, 2:7] <- round(FN_analysis_perc[, 2:7]/ FN_analysis_perc[, 8] * 100, 4)

write.table(FN_analysis, file.path(resultDir, "FN_analysis_all.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FN_analysis_perc, file.path(resultDir, "FN_analysis_all_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  




rois_gene <- data.frame(genes = c(singleData[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FN_single$V4 <- factor(FN_single$V4,  levels = rois_gene_count$genes)
FN_single$algorithmID <- factor(FN_single$algorithmID, levels = algorithms)

FN_analysis <- FN_single %>% group_by(V4, algorithmID, .drop = FALSE) %>% summarize(n=n())
FN_analysis <- cast(FN_analysis, V4 ~ algorithmID)
names(FN_analysis)[1] <- "genes"

FN_analysis <- merge(FN_analysis, rois_gene_count, by = "genes")
FN_analysis_perc <- FN_analysis
FN_analysis_perc[, 2:7] <- round(FN_analysis_perc[, 2:7]/ FN_analysis_perc[, 8] * 100, 4)

write.table(FN_analysis, file.path(resultDir, "FN_analysis_single.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FN_analysis_perc, file.path(resultDir, "FN_analysis_single_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  




rois_gene <- data.frame(genes = c(clinicData[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FN_clinic$V4 <- factor(FN_clinic$V4,  levels = rois_gene_count$genes)
FN_clinic$algorithmID <- factor(FN_clinic$algorithmID, levels = algorithms)

FN_analysis <- FN_clinic %>% group_by(V4, algorithmID, .drop = FALSE) %>% summarize(n=n())
FN_analysis <- cast(FN_analysis, V4 ~ algorithmID)
names(FN_analysis)[1] <- "genes"

FN_analysis <- merge(FN_analysis, rois_gene_count, by = "genes")
FN_analysis_perc <- FN_analysis
FN_analysis_perc[, 2:7] <- round(FN_analysis_perc[, 2:7]/ FN_analysis_perc[, 8] * 100, 4)

write.table(FN_analysis, file.path(resultDir, "FN_analysis_clinic.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FN_analysis_perc, file.path(resultDir, "FN_analysis_clinic_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  



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

