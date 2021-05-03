
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

write.table(allData, file.path(tempDir, "all.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicData, file.path(tempDir, "clinic.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleData, file.path(tempDir, "single.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  

write.table(allNormal, file.path(tempDir, "allNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicNormal, file.path(tempDir, "clinicNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleNormal, file.path(tempDir, "singleNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  

TP_all <- data.frame(matrix(ncol = 15, nrow = 0))
FP_all <- data.frame(matrix(ncol = 13, nrow = 0))

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
  
  # For all
  samples <- unique(algorithmData$sample)
  
  for (sample in samples){
    sample <- sub("X", "", sample)
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
    
    # read files
    if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
    if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}

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
  }
}


TP_single <- data.frame(matrix(ncol = 15, nrow = 0))
FP_single <- data.frame(matrix(ncol = 13, nrow = 0))

# All
cnvFiles <- list.files(singleDir, pattern = ".txt")
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
  
  # For all
  samples <- unique(algorithmData$sample)
  
  for (sample in samples){
    sample <- sub("X", "", sample)
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
    
    # read files
    if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
    if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}
    
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
  }
}


TP_clinic <- data.frame(matrix(ncol = 15, nrow = 0))
FP_clinic <- data.frame(matrix(ncol = 13, nrow = 0))

cnvFiles <- list.files(singleDir, pattern = ".txt")
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
  
  # For all
  samples <- unique(algorithmData$sample)
  
  for (sample in samples){
    sample <- sub("X", "", sample)
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
    
    # read files
    if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
    if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}
    
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
  }
}


# Compare del i dup

TP_all$V12 <- tolower(TP_all$V11)
TP_clinic$V12 <- tolower(TP_clinic$V11)
TP_single$V12 <- tolower(TP_single$V11)

a <- which(TP_all$V4 != TP_all$V12)
TP_all[a, ]

b <- which(TP_clinic$V4 != TP_clinic$V12)
TP_clinic[b, ]

c <- which(TP_single$V4 != TP_single$V12)
TP_single[c, ]

# FInal data frame
allData$cnvkit <- 0
allData$cnvkit2 <- 0
allData$cnvkit3 <- 0
allData$cnvkit4 <- 0
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
singleData$cnvkit <- 0
singleData$cnvkit2 <- 0
singleData$cnvkit3 <- 0
singleData$cnvkit4 <- 0
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
clinicData$cnvkit <- 0
clinicData$cnvkit2 <- 0
clinicData$cnvkit3 <- 0
clinicData$cnvkit4 <- 0
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

write.table(allData, file.path(resultDir, "all.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicData, file.path(resultDir, "clinic.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleData, file.path(resultDir, "single.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  

# False positives

# FInal data frame
allNormal$cnvkit <- 0
allNormal$cnvkit2 <- 0
allNormal$cnvkit3 <- 0
allNormal$cnvkit4 <- 0
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
singleNormal$cnvkit <- 0
singleNormal$cnvkit2 <- 0
singleNormal$cnvkit3 <- 0
singleNormal$cnvkit4 <- 0
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
clinicNormal$cnvkit <- 0
clinicNormal$cnvkit2 <- 0
clinicNormal$cnvkit3 <- 0
clinicNormal$cnvkit4 <- 0
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

write.table(allNormal, file.path(resultDir, "allNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicNormal, file.path(resultDir, "clinicNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleNormal, file.path(resultDir, "singleNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
