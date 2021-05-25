# Analysis of cnvfound results
library(dplyr)
library(reshape)
library(ggplot2)
library(gridExtra)

# Define dirs
analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
bedDir <- file.path(analysisDir, "bedfiles")
evaluationDir <- file.path(analysisDir, "evaluation")
tempDir <- file.path(evaluationDir, "temp")
resultDir <- file.path(evaluationDir, "results")
graphsDir <- file.path(evaluationDir, "graphs")

# Run over dataset and obtain list with FP, FN, TP and TN
for (dataset in c("all", "single")){
  bedFile <- file.path(bedDir, paste0(dataset, "_rois.bed"))
  bedData <- read.table(bedFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  
  # Generate separate bedfiles for positives and negatives
  positiveData <- subset(bedData, bedData$cnv == "ExonCNV")
  negativeData <- subset(bedData, bedData$cnv == "Normal")
  
  # Add ROI ID
  positiveData$ID <- 1:nrow(positiveData)
  negativeData$ID <- (1+nrow(positiveData)):nrow(bedData)
  
  # Import cnvfounds file
  cnvDir <- file.path(analysisDir, "cnvfounds", dataset)
  cnvFiles <- list.files(cnvDir, pattern = ".txt")
  
  # Create dataframes for result saving
  FPdata <- data.frame()
  TPdata <- data.frame()
  
  # Run over results algorithms
  for (cnvFile in cnvFiles){
    algorithm <- sub("cnvFounds_", "", cnvFile)
    algorithm <- sub(".txt", "", algorithm)
    
    # Import cnvfound data from algorithm
    algorithmData <- read.table(file.path(cnvDir, cnvFile), sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    
    # Change colnames
    names(algorithmData) <- tolower(names(algorithmData))
    names(algorithmData) <- sub("chromosome", "chr", names(algorithmData))
    names(algorithmData) <- sub("cnv.type", "type", names(algorithmData))
    
    # Select columns
    colNames <- c("chr", "start", "end", "sample","type")
    algorithmData <- algorithmData[, colNames] 
    
    # Eliminate X in sample name
    algorithmData$sample <- as.character(algorithmData$sample)
    algorithmData$sample <- sub("X", "", algorithmData$sample)
    
    # Run of sample and check for FP, FN, TP and TN
    samples <- unique(algorithmData$sample)
    for (sample in samples){
      
      # filter cnvfounds by sample
      s <- which(algorithmData$sample == sample)
      sampleData <- algorithmData[s, ]
      sampleFile <- file.path(tempDir, "sample.bed")
      write.table(sampleData, sampleFile, sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)  
      
      # Make positive roi and negative roi bed file for sample
      samplepositiveData <- subset(positiveData, positiveData$sampleID == sample)
      samplepositiveFile <- file.path(tempDir, "positive.bed")
      write.table(samplepositiveData, samplepositiveFile, sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)  
      
      samplenegativeData <- subset(negativeData, negativeData$sampleID == sample)
      samplenegativeFile <- file.path(tempDir, "negative.bed")
      write.table(samplenegativeData, samplenegativeFile, sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)  
      
      # Compare sampleCNVs with validated results to obtain TP and FP
      setwd(tempDir)
      system(paste("bedtools intersect -wa -wb -a", sampleFile, "-b", samplepositiveFile, "> TP.bed"))
      system(paste("bedtools intersect -wa -wb -a", sampleFile, "-b", samplenegativeFile, "> FP.bed"))
      
      # read files
      if (file.size("TP.bed") != 0) {TPsample <- read.table("TP.bed", sep = "\t", stringsAsFactors=FALSE)} else {TPsample <- NULL}
      if (file.size("FP.bed") != 0) {FPsample <- read.table("FP.bed", sep = "\t", stringsAsFactors=FALSE)} else {FPsample <- NULL}
  
      # Add found TP and FP to dataframe
      if(!is.null(TPsample)){
        n <- nrow(TPsample)
        algorithmID <- rep(algorithm, n)
        TPsample <- cbind(TPsample, algorithmID)
        TPdata <- rbind(TPdata, TPsample)
      }
      
      if(!is.null(FPsample)){
        n <- nrow(FPsample)
        algorithmID <- rep(algorithm, n)
        FPsample <- cbind(FPsample, algorithmID)
        FPdata <- rbind(FPdata, FPsample)
      }
    }
  }
  
  # Delete redundant columns
  FPdata <- FPdata[, c(6:17, 5)]
  TPdata <- TPdata[, c(6:17, 5)]
  colnames(FPdata) <- colnames(TPdata) <- c(colnames(bedData), "roiID", "algorithmID", "detected_cnv")
  PosData <- rbind(TPdata, FPdata)
  
  # Make final result dataframe
  ResultsDataframe <- rbind(positiveData, negativeData)
  ResultsDataframe$cnvkit5 <- 0
  ResultsDataframe$convading <- 0
  ResultsDataframe$decon <- 0
  ResultsDataframe$exomedepth <- 0
  ResultsDataframe$manta <- 0
  ResultsDataframe$panelcn <- 0
  
  algorithms <- unique(c(TPdata$algorithmID, FPdata$algorithmID))
  
  for(i in 1:nrow(ResultsDataframe)){
    for(j in algorithms){
      a <- which(PosData$algorithmID == j & PosData$roiID == i)
      if(length(a) > 0){
        ResultsDataframe[i, j] <- 1
      }
    }
  }
  
  # For hybdrid model
  HybridDataframe <- rbind(positiveData, negativeData)
  HybridDataframe$cnvkit5 <- 0
  HybridDataframe$convading <- 0
  HybridDataframe$decon <- 0
  HybridDataframe$exomedepth <- 0
  HybridDataframe$manta <- 0
  HybridDataframe$panelcn <- 0
  
  algorithms <- unique(c(TPdata$algorithmID, FPdata$algorithmID))
  
  # There is one ROI for DeCON (FP) where both DEL and DUP is predicted. This ROI is considered as normal for the hybrid model
  for(i in 1:nrow(HybridDataframe)){
    for(j in algorithms){
      a <- which(PosData$algorithmID == j & PosData$roiID == i)
      if(length(a) == 1){
        if(PosData$detected_cnv[a] == "deletion") {
          HybridDataframe[i, j] <- -1
      }
        else if(PosData$detected_cnv[a] == "duplication") {
          HybridDataframe[i, j] <- 1
        }
      }
    }
  }
  
  # Save and assign results to dataframe for dataset
  write.table(ResultsDataframe, file.path(resultDir, paste0(dataset, "results.txt")), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  
  write.table(HybridDataframe, file.path(resultDir, paste0(dataset, "hybridresults.txt")), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  
  assign(paste0(dataset, "Results"), ResultsDataframe)
}


PosDataFail <- PosData[, c("algorithmID", "roiID")]
PosDataFail_algo <- subset(PosDataFail, PosDataFail$algorithmID == "decon")
nrow(unique(PosDataFail_algo))

# Result table

## Add identifier to result dataframes
comment(allResults) <- "all"
comment(singleResults) <- "single"

for(resultframe in list(allResults, singleResults)) {
  
  # Obtain name dataframe
  resultName <- comment(resultframe) 
  
  # Calculate TP, FN, FP, TN
  Pos <- length(resultframe$cnv[resultframe$cnv == "ExonCNV"])
  TP <- as.vector(colSums(resultframe[1:Pos, 12:17]))
  FN <- as.vector(nrow(resultframe[1:Pos, ]) - as.vector(colSums(resultframe[1:Pos, 12:17])))
  FP <- as.vector(colSums(resultframe[(Pos+1):nrow(resultframe), 12:17]))
  TN <- as.vector(nrow(resultframe[(1+Pos):nrow(resultframe), ]) - as.vector(colSums(resultframe[(Pos+1):nrow(resultframe), 12:17])))
  
  # Calculate Specificity and sensitivity
  Total <- TP + FN + FP + TN
  Spec <- round(TN / (TN + FP), 4)
  Sens <- round(TP / (TP + FN), 4)
  
  # Export result table
  resultTable <- cbind(TP, TN, FP, FN, Total, Sens, Spec, names(resultframe[, 12:17]))
  write.table(resultTable, file.path(resultDir, paste0("summary_", resultName, ".txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
}

# Obtain the number of FP, FN ROIs at gene level
for(resultframe in list(allResults, singleResults)) {
  
  # Obtain name resultframe
  resultName <- comment(resultframe) 
  
  # Separate posive and negative ROIs
  positiveData <- subset(resultframe, resultframe$cnv == "ExonCNV")
  negativeData <- subset(resultframe, resultframe$cnv == "Normal")
  
  # Number of positive and negative ROIs per gene
  rois_gene_pos <- data.frame(gene = c(positiveData$gene))
  rois_gene_pos_count <- rois_gene_pos %>% group_by(gene) %>% summarise(n=n())
  rois_gene_neg <- data.frame(gene = c(negativeData$gene))
  rois_gene_neg_count <- rois_gene_neg %>% group_by(gene) %>% summarise(n=n())
  
  # Convert to factor
  positiveData$gene <- factor(positiveData$gene,  levels = rois_gene_pos_count$gene)
  negativeData$gene <- factor(negativeData$gene,  levels = rois_gene_neg_count$gene)
  
  # Obtain FN dataframe
  FNdataframe <- data.frame()
  for(algorithm in algorithms){
    positiveAlgorithm <- cbind(positiveData[, c(1:11)], positiveData[, algorithm])
    positiveAlgorithm <- subset(positiveAlgorithm, positiveAlgorithm[, 12] == 0)
    positiveAlgorithm[, 12] <- rep(algorithm, nrow(positiveAlgorithm))
    colnames(positiveAlgorithm)[12] <- "algorithmID"
    FNdataframe <- rbind(FNdataframe, positiveAlgorithm)
  } 
  
  # Obtain FP dataframe
  FPdataframe <- data.frame()
  for(algorithm in algorithms){
    negativeAlgorithm <- cbind(negativeData[, c(1:11)], negativeData[, algorithm])
    negativeAlgorithm <- subset(negativeAlgorithm, negativeAlgorithm[, 12] == 1)
    negativeAlgorithm[, 12] <- rep(algorithm, nrow(negativeAlgorithm))
    colnames(negativeAlgorithm)[12] <- "algorithmID"
    FPdataframe <- rbind(FPdataframe, negativeAlgorithm)
  } 
  
  # Convert algorithmID to factor
  FNdataframe$algorithmID <- factor(FNdataframe$algorithmID, levels = algorithms)
  FPdataframe$algorithmID <- factor(FPdataframe$algorithmID, levels = algorithms)
  
  # Analysis FN
  FN_analysis <- FNdataframe %>% group_by(gene, algorithmID, .drop = FALSE) %>% summarize(n=n())
  FN_analysis <- cast(FN_analysis, gene ~ algorithmID)
  
  # Add total ROI and calculate %
  FN_analysis <- merge(FN_analysis, rois_gene_neg_count, by = "gene")
  FN_analysis_perc <- FN_analysis
  FN_analysis_perc[, 2:7] <- round(FN_analysis_perc[, 2:7] / FN_analysis_perc[, 8] * 100, 4)
  
  write.table(FN_analysis, file.path(resultDir, paste0(resultName, "_FN.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  write.table(FN_analysis_perc, file.path(resultDir, paste0(resultName, "_FN_perc.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  
  # Analysis FP
  FP_analysis <- FPdataframe %>% group_by(gene, algorithmID, .drop = FALSE) %>% summarize(n=n())
  FP_analysis <- cast(FP_analysis, gene ~ algorithmID)
  
  # Add total ROI and calculate %
  FP_analysis <- merge(FP_analysis, rois_gene_pos_count, by = "gene")
  FP_analysis_perc <- FP_analysis
  FP_analysis_perc[, 2:7] <- round(FP_analysis_perc[, 2:7] / FP_analysis_perc[, 8] * 100, 4)
  
  write.table(FP_analysis, file.path(resultDir, paste0(resultName, "_FP.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  write.table(FP_analysis_perc, file.path(resultDir, paste0(resultName, "_FP_perc.txt")), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  
  # Assign results to object
  assign(paste0(resultName, "_FN"), FN_analysis)
  assign(paste0(resultName, "_FN_perc"), FN_analysis_perc)
  assign(paste0(resultName, "_FP"), FP_analysis)
  assign(paste0(resultName, "_FP_perc"), FP_analysis_perc)
  
  # Assign FP and FN lists to object
  assign(paste0(resultName, "_FPlist"), FNdataframe)
  assign(paste0(resultName, "_FNlist"), FPdataframe)
}


# Compare differences in length, single / multi and del / dup between FP and TP, and FN and TP

p_list_FN <- list()
p_list_FP <- list()

FN_comp_length <- data.frame()
FN_comp_exontype <- data.frame()
FN_comp_cnvtype <- data.frame()

FP_comp_length <- data.frame()

# Loop over algorithms
for(algorithm in algorithms){
  
  # TP and TN for all Dataset
  allResults$cnv_type <- factor(allResults$cnv_type, levels = c("Deletion", "Duplication"))
  allResults$exon_type <- factor(allResults$exon_type, levels = c("Single", "Multi"))
  positiveData <- subset(allResults, allResults$cnv == "ExonCNV")
  negativeData <- subset(allResults, allResults$cnv == "Normal")
  
  # Select positive cases for algorithm
  positiveAlgorithm <- cbind(positiveData[, c(1:11)], detected = positiveData[, algorithm])
  positiveAlgorithm$detected[positiveAlgorithm$detected == 1] <- "TP"
  positiveAlgorithm$detected[positiveAlgorithm$detected == 0] <- "FN"
  positiveAlgorithm$detected <- factor(positiveAlgorithm$detected)
  positiveAlgorithm[, 13] <- rep(algorithm, nrow(positiveAlgorithm))
  colnames(positiveAlgorithm)[13] <- "algorithmID"
  
  # Select negative cases for algorithm
  negativeAlgorithm <- cbind(negativeData[, c(1:11)], detected = negativeData[, algorithm])
  negativeAlgorithm$detected[negativeAlgorithm$detected == 1] <- "FP"
  negativeAlgorithm$detected[negativeAlgorithm$detected == 0] <- "TN"
  negativeAlgorithm$detected <- factor(negativeAlgorithm$detected)
  negativeAlgorithm[, 13] <- rep(algorithm, nrow(negativeAlgorithm))
  colnames(negativeAlgorithm)[13] <- "algorithmID"
  
  # Compare TP and FN
  
  ## Compare PB_length
  comp_length <- positiveAlgorithm %>% group_by(detected) %>% summarize(mean = mean(pb_length), 
                                                         min = min(pb_length), max = max(pb_length),
                                                         median = median(pb_length),
                                                         IQR25 = quantile(pb_length, 0.25), 
                                                         IQR75 = quantile(pb_length, 0.75))
  
  comp_length$algorithm <- rep(algorithm, nrow(comp_length))
  
  ## Make distribution profile PB_length
  p <- ggplot(positiveAlgorithm, aes(log(pb_length), color = detected, fill = detected)) + 
    geom_density(alpha = 0.2) + 
    ggtitle(algorithm) +
    xlab("log(Base pairs)") +
    theme_minimal()
  
  p_list_FN[[algorithm]] <- p
  
  ## Compare Exontype
  comp_exontype <- positiveAlgorithm %>% group_by(detected, exon_type, .drop = FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
  comp_exontype$algorithm <- rep(algorithm, nrow(comp_exontype))
  
  ## Compare cnvtype
  comp_cnvtype <- positiveAlgorithm %>% group_by(detected, cnv_type, .drop = FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
  comp_cnvtype$algorithm <- rep(algorithm, nrow(comp_cnvtype))
  
  # Save results in common dataframe
  FN_comp_length <- rbind(FN_comp_length, comp_length)
  FN_comp_exontype <- rbind(FN_comp_exontype, comp_exontype)
  FN_comp_cnvtype <- rbind(FN_comp_cnvtype, comp_cnvtype)
  
  # Compare TP and FP
  
  ## Compare PB_length
  TP_FP_Algorithm <- rbind(subset(positiveAlgorithm, positiveAlgorithm$detected == "TP"), subset(negativeAlgorithm, negativeAlgorithm$detected == "FP"))
   
  comp_length <- TP_FP_Algorithm %>% group_by(detected) %>% summarize(mean = mean(pb_length), 
                                                                        min = min(pb_length), max = max(pb_length),
                                                                        median = median(pb_length),
                                                                        IQR25 = quantile(pb_length, 0.25), 
                                                                        IQR75 = quantile(pb_length, 0.75))
  
  comp_length$algorithm <- rep(algorithm, nrow(comp_length))
  
  ## Make distribution profile PB_length
  p <- ggplot(TP_FP_Algorithm, aes(log(pb_length), color = detected, fill = detected)) + 
    geom_density(alpha = 0.2) + 
    ggtitle(algorithm) +
    xlab("log(Base pairs)") +
    theme_minimal()
  
  p_list_FP[[algorithm]] <- p
  
  # Save results in common dataframe
  FP_comp_length <- rbind(FP_comp_length, comp_length)
}


# PB length plots
setwd(graphsDir)
tiff("BP_FP.tiff", units="in", width=7, height=5, res=150)

n <- length(p_list_FP)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(p_list_FP, ncol=nCol))

dev.off()

tiff("BP_FN.tiff", units="in", width=7, height=5, res=150)

n <- length(p_list_FN)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(p_list_FN, ncol=nCol))

dev.off()

# Plots for FN per exontype and cnvtype
FN_cnv_plot <- subset(FN_comp_cnvtype, FN_comp_cnvtype$detected == "FN" & FN_comp_cnvtype$algorithm != "manta")
FN_cnv_plot$detected <- factor(FN_cnv_plot$detected)
FN_cnv_plot$algorithm <- factor(FN_cnv_plot$algorithm)

tiff("FN_cnv.tiff", units="in", width=7, height=5, res=150)

ggplot(FN_cnv_plot, aes(x=algorithm, y=freq*100, fill=cnv_type, color=cnv_type))+
  geom_bar(stat="identity",position=position_dodge(), alpha = 0.2)+
  xlab("Algorithm")+
  scale_y_continuous(limits = c(0,100))+
  ylab("Frequency per type")+ 
  geom_text(aes(label=n, fontface="bold"), position=position_dodge(width = 0.9), vjust=-0.25, size=4)+
  theme_minimal()

dev.off()

FN_exon_plot <- subset(FN_comp_exontype, FN_comp_exontype$detected == "FN" & FN_comp_exontype$algorithm != "manta")
FN_exon_plot$detected <- factor(FN_exon_plot$detected)
FN_exon_plot$algorithm <- factor(FN_exon_plot$algorithm)

tiff("FN_exon.tiff", units="in", width=7, height=5, res=150)

ggplot(FN_exon_plot, aes(x=algorithm, y=freq*100, fill=exon_type, color=exon_type))+
  geom_bar(stat="identity",position=position_dodge(), alpha = 0.2) +
  xlab("Algorithm")+
  scale_y_continuous(limits = c(0,100))+
  ylab("Frequency per type")+ 
  geom_text(aes(label=n, fontface="bold"), position=position_dodge(width = 0.9), vjust=-0.25, size=4)+
  theme_minimal()

dev.off()


