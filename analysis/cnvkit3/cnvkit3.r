# CNVKIT3 analysis

setwd("/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/cnvfounds/cnvkit3")
dir <- getwd()
analysisDir <- "/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/cnvfounds"


cnvFounds <- read.table(file.path(dir, "cnvFounds_cnvkit5.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)

# load bed files positive samples
allFile <- file.path(analysisDir, "all_positives.bed")
allData <- read.table(allFile, sep = "\t", stringsAsFactors=FALSE)
allNormalFile <- file.path(analysisDir, "all_normal.bed")
allNormal <- read.table(allNormalFile, sep = "\t", stringsAsFactors=FALSE)

clinicFile <- file.path(analysisDir, "clinic_positives.bed")
clinicData <- read.table(clinicFile,  sep = "\t", stringsAsFactors=FALSE)
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


# Prepare cnvFounds
cnvFounds <- cnvFounds[, c(3:9, 1:2)]

# Obtain FP and TP cnvkit3
TP_all <- data.frame(matrix(ncol = 19, nrow = 0))
FP_all <- data.frame(matrix(ncol = 17, nrow = 0))

samples <- as.character(sort(unique(c(allNormal$V4, allData$V4))))

for (sample in samples){
  sample <- sub("X", "", sample)
  s <- which(cnvFounds$Sample == sample)
  sampleData <- cnvFounds[s, ]
  sampleBed <- file.path(dir, paste0("cnvkit3_", sample, ".bed"))
  write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
  
  # Prepare sample bed file, filter by sample
  allSample <- subset(allData, allData$V4 == sample)
  allBed <- file.path(dir, "all.bed")
  write.table(allSample, allBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
  
  allNormalsample <- subset(allNormal, allNormal$V4== sample)
  allNormalbed <- file.path(dir, "allNormal.bed")
  write.table(allNormalsample, allNormalbed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
  
  # intersect
  setwd(dir)
  system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", allBed, "> TP.bed"))
  system(paste("bedtools intersect -wa -wb -a", sampleBed, "-b", allNormalbed, "> FP.bed"))
  system(paste("bedtools intersect -v -wa -a", allBed, "-b", sampleBed, "> FN.bed"))
  
  # read files
  if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
  if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}
  
  # Edit files
  if(!is.null(TPsample)){
    n <- nrow(TPsample)
    sampleID <- rep(sample, n)
    TPsample <- cbind(sampleID, TPsample)
    TP_all <- rbind(TP_all, TPsample)
  }
  
  if(!is.null(FPsample)){
    n <- nrow(FPsample)
    sampleID <- rep(sample, n)
    FPsample <- cbind(sampleID, FPsample)
    FP_all <- rbind(FP_all, FPsample)
  }
}

weight_FP <- FP_all$V8
weight_TP <- TP_all$V8
weights <- data.frame(Value = c(weight_FP, weight_TP), Grup = c(rep("FP", length(weight_FP)), rep("TP", length(weight_TP))))


log_FP <- FP_all$V4
log_TP <- TP_all$V4
logs <- data.frame(Value = c(log_FP, log_TP), Grup = c(rep("FP", length(log_FP)), rep("TP", length(log_TP))))


library(ggplot2)
library(dplyr)


ggplot(weights, aes(Value, color = Grup, fill = Grup)) + 
  geom_density(alpha = 0.2) + 
  theme_minimal()

ggplot(logs, aes(Value, color = Grup, fill = Grup)) + 
  geom_density(alpha = 0.2) + 
  theme_minimal()


# Logs summary
Del <- subset(logs, logs$Value < 0)
Dup <- subset(logs, logs$Value > 0)

DelDF <- Del %>% group_by(Grup) %>% summarise(mean = mean(Value), median = median(Value), min = min(Value), max = max(Value), IQR25 = quantile(Value, 0.25), IQR75 = quantile(Value, 0.75))
DupDF <- Dup %>% group_by(Grup) %>% summarise(mean = mean(Value), median = median(Value), min = min(Value), max = max(Value), IQR25 = quantile(Value, 0.25), IQR75 = quantile(Value, 0.75))

