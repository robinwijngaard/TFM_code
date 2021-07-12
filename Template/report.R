library(gggenes)
library(ggplot2)
library(rmarkdown)
library(DT)
library(ExomeDepth)
library(dplyr)
library(tiff)
library(patchwork)
library(RColorBrewer)
library(R.utils)
library(grid)
library(gridExtra)
library(stringr)

args=commandArgs(asValues = TRUE)

# Define dir
outputDir <- args$outputDir
annotatedFile <- args$annotatedFile
panellsFile <- args$panellsFile
samplesFile <- args$samplesFile
RDSfolder <- args$RDSfolder
markdownfile <- args$markdownfile

# make temp dir, results dir en graphs dir
tempDir <- file.path(outputDir, "temp")
resultDir <- file.path(outputDir, "results")
graphsDir <- file.path(outputDir, "graphs")
reportsDir <- file.path(outputDir, "reports")

if(!dir.exists(tempDir)) dir.create(tempDir)
if(!dir.exists(resultDir)) dir.create(resultDir)
if(!dir.exists(graphsDir)) dir.create(graphsDir)
if(!dir.exists(reportsDir)) dir.create(reportsDir)

# Read annotatedFile
annotatedFile <- read.delim(annotatedFile, sep="\t", header = TRUE)  

# Read samples and linked panel
samplesData <- read.delim(samplesFile, header=FALSE)

# Run over samples and create samplereport
for (i in 1:nrow(samplesData)){
  render(markdownfile, output_format = "html_document", params = list(sample = samplesData$V1[i], panell = samplesData$V2[i]), output_file = file.path(reportsDir, paste0(samplesData$V1[i],"report.html")))
}



outputDir="/Users/robinwijngaard/Dropbox/Master_UOC/TFM/PostTFM/clinic_template"#outputDir
annotatedFile="/Users/robinwijngaard/Dropbox/Master_UOC/TFM/PostTFM/clinic_template/annotatedFile.bed" #annotatedfile (output annotate.R)
panellsFile="/Users/robinwijngaard/Dropbox/Master_UOC/TFM/PostTFM/gens_nm_panell.xlsx" #file idnicating genes per panell and NM
samplesFile="/Users/robinwijngaard/Dropbox/Master_UOC/TFM/PostTFM/clinic_template/samples.txt" #file indicating sample and panell association
RDSfolder="/Users/robinwijngaard/Dropbox/Master_UOC/TFM/PostTFM/clinic_template/rds" #folder containing output RDS files from the ExomeDepth algorithm#
markdownfile="/Users/robinwijngaard/Dropbox/Master_UOC/TFM/TFM_code/Template/report.Rmd" #rmarkdown file


