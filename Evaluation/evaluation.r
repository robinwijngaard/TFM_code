# Analysis of cnvfound results
library(dplyr)
library(reshape)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(R.utils)

args=commandArgs(asValues = TRUE)

# Define dir
outputDir <- args$outputDir
cnvFounds <- args$cnvFounds
bedFile <- args$bedFile
roisFile <- args$roisFile

# make temp dir, results dir en graphs dir
tempDir <- file.path(outputDir, "temp")
resultDir <- file.path(outputDir, "results")
graphsDir <- file.path(outputDir, "graphs")

if(!dir.exists(tempDir)) dir.create(tempDir)
if(!dir.exists(resultDir)) dir.create(resultDir)
if(!dir.exists(graphsDir)) dir.create(graphsDir)

# Run over dataset and obtain list with FP, FN, TP and TN
roisData <- read.table(roisFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  
# Generate separate bedfiles for positives and negatives
positiveData <- subset(roisData, roisData$cnv == "ExonCNV")
negativeData <- subset(roisData, roisData$cnv == "Normal")
  
# Add ROI ID
positiveData$ID <- 1:nrow(positiveData)
negativeData$ID <- (1+nrow(positiveData)):nrow(roisData)
  
# Import cnvfounds file
cnvFiles <- list.files(cnvFounds, pattern = ".txt")
  
# Create dataframes for result saving
FPdata <- data.frame()
TPdata <- data.frame()
  
# Run over results algorithms
for (cnvFile in cnvFiles){
  algorithm <- sub("cnvFounds_", "", cnvFile)
  algorithm <- sub(".txt", "", algorithm)
    
  # Import cnvfound data from algorithm
  algorithmData <- read.table(file.path(cnvFounds, cnvFile), sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    
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
    
  # Run over sample and check for FP, FN, TP and TN
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
FPdata <- FPdata[, c(6:ncol(FPdata), 5)]
TPdata <- TPdata[, c(6:ncol(TPdata), 5)]
colnames(FPdata) <- colnames(TPdata) <- c(colnames(roisData), "roiID", "algorithmID", "detected_cnv")
PosData <- rbind(TPdata, FPdata)
  
# Make final result dataframe
resultsData <- rbind(positiveData, negativeData)
algorithms <- unique(c(TPdata$algorithmID, FPdata$algorithmID))

for (i in 1:length(algorithms)){
  algorithmColumn <- rep(0, nrow(resultsData))
  resultsData <- cbind(resultsData, algorithmColumn)
  names(resultsData)[ncol(resultsData)] <- algorithms[i]
}

# Add 1 if del or dup detected
for(i in 1:nrow(resultsData)){
  for(j in algorithms){
    a <- which(PosData$algorithmID == j & PosData$roiID == i)
    if(length(a) > 0){
      resultsData[i, j] <- 1
    }
  }
}
  
# Idem but with del = -1 and dup = 1
combinedData <- rbind(positiveData, negativeData)

for (i in 1:length(algorithms)){
  algorithmColumn <- rep(0, nrow(combinedData))
  combinedData <- cbind(combinedData, algorithmColumn)
  names(combinedData)[ncol(combinedData)] <- algorithms[i]
}

# There is one ROI for DeCON (FP) where both DEL and DUP is predicted. This ROI is considered as normal for the hybrid model
for(i in 1:nrow(combinedData)){
  for(j in algorithms){
    a <- which(PosData$algorithmID == j & PosData$roiID == i)
    if(length(a) == 1){
      if(PosData$detected_cnv[a] == "deletion") {
        combinedData[i, j] <- -1
    }
      else if(PosData$detected_cnv[a] == "duplication") {
        combinedData[i, j] <- 1
      }
    }
  }
}
  
# Save and assign results to dataframe for dataset
write.table(resultsData, file.path(resultDir, "resultData.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  
write.table(combinedData, file.path(resultDir, "combinedData.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  

# Result table

# Calculate TP, FN, FP, TN
Pos <- length(resultsData$cnv[resultsData$cnv == "ExonCNV"])
start <- 12
end <- 12 + length(algorithms) -1

TP <- as.vector(colSums(resultsData[1:Pos, start:end]))
FN <- as.vector(nrow(resultsData[1:Pos, ]) - as.vector(colSums(resultsData[1:Pos, start:end])))
FP <- as.vector(colSums(resultsData[(Pos+1):nrow(resultsData), start:end]))
TN <- as.vector(nrow(resultsData[(1+Pos):nrow(resultsData), ]) - as.vector(colSums(resultsData[(Pos+1):nrow(resultsData), start:end])))
  
# Calculate Specificity and sensitivity
Total <- TP + FN + FP + TN
Spec <- round(TN / (TN + FP), 4)
Sens <- round(TP / (TP + FN), 4)
Accuracy <- round((TP + TN)/(Total), 4)
  
# Export result table
resultTable <- as.data.frame(cbind(TP, TN, FP, FN, Total, Sens, Spec, Accuracy, names(resultsData[, start:end])))
write.table(resultTable, file.path(resultDir, "result_table.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  

# Obtain the number of FP, FN ROIs at gene level

# Separate posive and negative ROIs
positiveData <- subset(resultsData, resultsData$cnv == "ExonCNV")
negativeData <- subset(resultsData, resultsData$cnv == "Normal")
  
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
FN_analysis <- merge(FN_analysis, rois_gene_pos_count, by = "gene")
nlength <-  length(algorithms)+1
FN_analysis_perc <- FN_analysis
FN_analysis_perc[, 2:nlength] <- round(FN_analysis_perc[, 2:nlength] / FN_analysis_perc[, (nlength+1)] * 100, 4)

write.table(FN_analysis, file.path(resultDir, "FN_gene.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FN_analysis_perc, file.path(resultDir,"FN_by_gene.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  

# Analysis FP
FP_analysis <- FPdataframe %>% group_by(gene, algorithmID, .drop = FALSE) %>% summarize(n=n())
FP_analysis <- cast(FP_analysis, gene ~ algorithmID)
  
# Add total ROI and calculate %
FP_analysis <- merge(FP_analysis, rois_gene_neg_count, by = "gene")
FP_analysis_perc <- FP_analysis
FP_analysis_perc[, 2:nlength] <- round(FP_analysis_perc[, 2:nlength] / FP_analysis_perc[, (nlength+1)] * 100, 4)

write.table(FP_analysis, file.path(resultDir, "FP_by_gene.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FP_analysis_perc, file.path(resultDir, "FP_by_gene_perc.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
  
# Compare differences in length, single / multi and del / dup between FP and TP, and FN and TP
p_list_FN <- list()
p_list_FP <- list()

FN_comp_length <- data.frame()
FN_comp_type <- data.frame()
FP_comp_length <- data.frame()

# Delete manta
algorithms <- algorithms[-5]

# Loop over algorithms
for(algorithm in algorithms){
  
  # TP and TN for all Dataset
  resultsData$cnv_type <- factor(resultsData$cnv_type, levels = c("Deletion", "Duplication"))
  resultsData$exon_type <- factor(resultsData$exon_type, levels = c("Single", "Multi"))
  positiveData <- subset(resultsData, resultsData$cnv == "ExonCNV")
  negativeData <- subset(resultsData, resultsData$cnv == "Normal")
  
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
  
  ## Compare type
  comp_type <- positiveAlgorithm %>% group_by(detected, cnv_type, exon_type, .drop = FALSE) %>% summarise(n = n())
  comp_type$algorithm <- rep(algorithm, nrow(comp_type))
  
  # Save results in common dataframe
  FN_comp_length <- rbind(FN_comp_length, comp_length)
  FN_comp_type <- rbind(FN_comp_type, comp_type)

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

do.call("grid.arrange", c(p_list_FP, ncol=2))

dev.off()

tiff("BP_FN.tiff", units="in", width=7, height=5, res=150)

do.call("grid.arrange", c(p_list_FN, ncol=2))

dev.off()

# Plots for FN per exontype and cnvtype
FN_cnv_plot <- subset(FN_comp_type, FN_comp_type$detected == "FN" & FN_comp_type$algorithm != "manta" & FN_comp_type$n != 0)
FN_cnv_plot$detected <- factor(FN_cnv_plot$detected)
FN_cnv_plot$algorithm <- factor(FN_cnv_plot$algorithm)
FN_cnv_plot$type <- paste(FN_cnv_plot$cnv_type, FN_cnv_plot$exon_type)

tiff("FN_cnv.tiff", units="in", width=9, height=5, res=150)

colors <- brewer.pal(n = 8, name = "Spectral")[c(3,8)]

algo.labs <- c("CNVkit", "CoNVaDING", "DECoN", "ExomeDepth", "panelcn.MOPS")
names(algo.labs) <- levels(FN_cnv_plot$algorithm)

ggplot(FN_cnv_plot, aes(x=cnv_type, y=n, fill=exon_type, color=exon_type))+
  geom_bar(stat="identity",position="stack", alpha = 0.2)+
  facet_grid( ~algorithm, labeller = labeller(algorithm = algo.labs)) +
  scale_y_continuous(limits = c(0,14))+
  ylab("Number of false negatives")+ 
  scale_x_discrete(name ="CNV type", labels = c("Del", "Dup")) +
  theme_classic() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  theme(legend.title = element_blank(), text = element_text(size=16),
          axis.text.x = element_text(size=12))+
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors)

dev.off()
