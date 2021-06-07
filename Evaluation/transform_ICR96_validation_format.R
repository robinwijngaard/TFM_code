# Obtain ROI files for datasets
library(readxl)
library(xlsx)
library(readr)
library(R.utils)

args=commandArgs(asValues = TRUE)

# Define dir
outputDir <- args$outputDir
bedFile <- args$bedFile
mlpaFile <- args$mlpaFile

# make temp dir, results dir en graphs dir
tempDir <- file.path(outputDir, "temp")
resultDir <- file.path(outputDir, "results")
graphsDir <- file.path(outputDir, "graphs")

if(!dir.exists(tempDir)) dir.create(tempDir)
if(!dir.exists(resultDir)) dir.create(resultDir)
if(!dir.exists(graphsDir)) dir.create(graphsDir)

# Obtain ROI for negative validated genes
MLPA <- read_excel(mlpaFile, sheet = "negatives")

peti <- MLPA$`NÚMERO DE PETICIÓ`
mlpas <- MLPA[, -c(2)]

negatives <- data.frame()

# Check for negatively validated genes
for (i in 1:nrow(mlpas)){
  for (j in 2:ncol(mlpas)){
    if(!is.na(as.character(mlpas[i,j]))){
      negatives <- rbind(negatives, c(as.character(mlpas[i, 1]), as.character(names(mlpas)[j])))
    }
  }
}

colnames(negatives) <- c("peti", "gen")

# Merge by gene
bedData <- read.table(bedFile, sep = "\t", stringsAsFactors=FALSE)
validatedNegative <- merge(negatives, bedData, by.x = "gen", by.y = "V4")

# Fix negatives bed file
validatedNegative <- validatedNegative[, c(3,4,5,1,2)]
names(validatedNegative) <- c("chr", "start", "end", "gene", "sampleID")

# Add missing columns
validatedNegative$cnv_description <- rep("Normal", nrow(validatedNegative))
validatedNegative$cnv <- rep("Normal", nrow(validatedNegative))
validatedNegative$cnv_type <- rep("", nrow(validatedNegative))
validatedNegative$exon_type <- rep("", nrow(validatedNegative))

# Load and fix positives
validatedPositive <- read_excel(mlpaFile, sheet = "positives", col_names = FALSE)

# Prepare positive data
validatedPositive$cnv <- rep("ExonCNV", nrow(validatedPositive))
validatedPositive <- validatedPositive[, c(1,2,3,5,4,6,10, 7,8)]
names(validatedPositive) <- c("chr", "start", "end", "gene", "sampleID", "cnv_description", "cnv", "cnv_type", "exon_type")

# Combine
validatedData <- rbind(validatedPositive, validatedNegative)
write.table(validatedData, file.path(outputDir, "validated_ICR96format.bed"), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE) 