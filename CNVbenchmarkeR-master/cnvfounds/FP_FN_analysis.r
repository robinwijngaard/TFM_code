
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
    system(paste("bedtools intersect -wa -v -a", allBed, "-b", sampleBed, "> FN.bed"))
    
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
    system(paste("bedtools intersect -wa -v -a", allBed, "-b", sampleBed, "> FN.bed"))
    
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
  }
}


TP_clinic <- data.frame(matrix(ncol = 15, nrow = 0))
FP_clinic <- data.frame(matrix(ncol = 13, nrow = 0))
FN_clinic <- data.frame(matrix(ncol = 11, nrow = 0))

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
    system(paste("bedtools intersect -wa -v -a", allBed, "-b", sampleBed, "> FN.bed"))
    
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

write.table(allData, file.path(resultDir, "all.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicData, file.path(resultDir, "clinic.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleData, file.path(resultDir, "single.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  

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

write.table(allNormal, file.path(resultDir, "allNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(clinicNormal, file.path(resultDir, "clinicNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
write.table(singleNormal, file.path(resultDir, "singleNormal.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  

# VennDiagram
setwd(resultDir)

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

library("VennDiagram")

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

# MDS
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
legend(1.5, 2, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)

##Normal
TN_mds <- rbind(allNormal[, 7:13])
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

# Cluster
clust <- rbind(allData[, 9:15], allNormal[, 7:13])
clust <- clust[, c(-1,-6)]

d <- dist(t(clust), method = "euclidean")
all.hc <- hclust(d, method = "ward.D2")
plot(all.hc)


# ROIS per gene

rois_gene <- data.frame(genes = c(allData[, 4], allNormal[, 4]))
rois_gene_count <- rois_gene %>% group_by(genes) %>% summarise(n=n())

FP_analysis <- FP_all %>% group_by(V8, algorithmID) %>% summarize(n=n())
FN_analysis <- FN_all %>% group_by(V4, algorithmID) %>% summarize(n=n())

FP_analysis$roicount <- NA

for(i in 1:nrow(FP_analysis)){
  for(j in 1:nrow(rois_gene_count)){
    if(FP_analysis$V8[i] == rois_gene_count$genes[j]){
      FP_analysis$roicount[i] <- rois_gene_count$n[j]
    }
  }
}

FN_analysis$roicount <- NA

for(i in 1:nrow(FN_analysis)){
  for(j in 1:nrow(rois_gene_count)){
    if(FN_analysis$V4[i] == rois_gene_count$genes[j]){
      FN_analysis$roicount[i] <- rois_gene_count$n[j]
    }
  }
}

FP_analysis$perc <- FP_analysis$n/FP_analysis$roicount * 100
FN_analysis$perc <- FN_analysis$n/FN_analysis$roicount * 100

write.table(FP_analysis, file.path(resultDir, "FP_analysis.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FN_analysis, file.path(resultDir, "FN_analysis.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)

# Search for specific ROIs frequently FP / FN

FP_analysis <- FP_all %>% group_by(V11, V8) %>% summarize(n=n())
FN_analysis <- FN_all %>% group_by(V9, V4) %>% summarize(n=n())

write.table(FP_analysis, file.path(resultDir, "FP_analysis_ROI.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FN_analysis, file.path(resultDir, "FN_analysis_ROI.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)


# Failed ROIS vs FP

deconFailedFile <- file.path(analysisDir, "all", "failedROIs_decon.csv")
deconFailed <- read.table(deconFailedFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
deconFailed <- deconFailed[, c(2:4,1)]

samples <- sort(unique(deconFailed$SampleID)) 

Failed_all <- data.frame(matrix(ncol = 10, nrow = 0))

for(sample in samples){
  s <- which(deconFailed$SampleID == sample)
  sampleData <- deconFailed[s, ]
  sampleBed <- file.path(tempDir, paste0("deconfailed_", sample, ".bed"))
  write.table(sampleData, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
  
  # Prepare sample bed file, filter by sample
  allSample <- subset(allData[, c(1:5,9)], allData$V4 == sample)
  allNormalsample <- subset(allNormal[, c(1:5, 7)], allNormal$V4== sample)
  colnames(allSample) <- c("chr", "start", "end", "gene", "sample", "ID")
  colnames(allNormalsample) <- c("chr", "start", "end", "gene", "sample", "ID")
  
  allcompletesample <- rbind(allSample, allNormalsample)
  
  allBed <- file.path(tempDir, "all.bed")
  write.table(allcompletesample, allBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
  
  # intersect
  setwd(tempDir)
  system(paste("bedtools intersect -wb -a", sampleBed, "-b", allBed, "> failed.bed"))
  
  if (file.size("failed.bed") != 0) {FailedSample <- read.table("failed.bed", sep = "\t", stringsAsFactors=FALSE)} else {FailedSample <- NULL}
  
  # Edit files
  if(!is.null(FailedSample)){
    Failed_all <- rbind(Failed_all, FailedSample)
  }
}

Failed_all <- Failed_all[, 5:10]
length(Failed_all$V10[Failed_all$V10 < 292])

FP_cnvkit <- subset(FP_all, FP_all$algorithmID == "cnvkit5")
FP_cnvkit <- FP_cnvkit[, c(3:6,1, 13)]

write.table(Failed_all, file.path(resultDir, "failed_decon.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
write.table(FP_cnvkit, file.path(resultDir, "fp_cnvkit.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 

setwd(tempDir)

Failed_FP.df <- data.frame(matrix(ncol = 12, nrow = 0))

for(sample in samples){
  s <- which(Failed_all$V9 == sample)
  sampleFailed <- Failed_all[s, ]
  sampleBed <- file.path(tempDir, paste0("deconfailed_", sample, ".bed"))
  write.table(sampleFailed, sampleBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)  
  
  # Prepare sample FP cnvkit
  fpSample <- subset(FP_cnvkit, FP_cnvkit$sampleID == sample)
  
  fpBed <- file.path(tempDir, paste0("cnvkitfp_", sample, ".bed"))
  write.table(fpSample, fpBed, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
  
  # intersect
  setwd(tempDir)
  system(paste("bedtools intersect -wb -wa -a", sampleBed, "-b", fpBed, "> failed_fp.bed"))
  
  if (file.size("failed_fp.bed") != 0) {FailedFPSample <- read.table("failed_fp.bed", sep = "\t", stringsAsFactors=FALSE)} else {FailedFPSample <- NULL}
  
  # Edit files
  if(!is.null(FailedFPSample)){
    Failed_FP.df <- rbind(Failed_FP.df, FailedFPSample)
  }
}

## Venn
venn_list <- list(Failed_decon = as.character(sort(unique(Failed_all$V10))),
                  FP_cnvkit = as.character(sort(unique(FP_cnvkit$V11))))

venn.diagram(venn_list, filename = "venn_failed_fp.tiff", col = "black",fill = c("dodgerblue", "goldenrod1"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

write.table(Failed_all %>% group_by(V9) %>% summarize(n=n()), file.path(resultDir, "Failed_samples.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
write.table(FP_cnvkit %>% group_by(sampleID) %>% summarize(n=n()), file.path(resultDir, "FP_samples.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 
