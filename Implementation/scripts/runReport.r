library(R.utils)
library(gggenes)
library(ggplot2)
library(rmarkdown)
library(DT)
library(ExomeDepth)
library(dplyr)
library(tiff)
library(patchwork)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(stringr)

# Read args
args <- commandArgs(asValues = TRUE)

## sample report files
annotatedFile <- args$annotatedFile
panellsFile <- args$panellsFile
samplesFile <- args$samplesFile
markdownfile <- args$markdownfile

resultDir <- args$resultDir
run <- args$run

# build output folder and file
runDir <- file.path(resultDir, run)
rdsDir <- file.path(runDir, "RDSfiles")
tempDir <- file.path(runDir, "temp")
graphsDir <- file.path(tempDir, "graphs")
reportsDir <- file.path(runDir, "reports")

# make temp dir, results dir en graphs dir
if(!dir.exists(resultDir)) dir.create(resultDir)
if(!dir.exists(runDir)) {dir.create(runDir)}
if(!dir.exists(rdsDir)) {dir.create(rdsDir)}
if(!dir.exists(tempDir)) dir.create(tempDir)
if(!dir.exists(graphsDir)) dir.create(graphsDir)
if(!dir.exists(reportsDir)) dir.create(reportsDir)

############# Sample template report

# Read annotatedFile
annotatedFile <- read.delim(annotatedFile, sep="\t", header = TRUE)  

# Read samples and linked panel
samplesData <- read.delim(samplesFile, header=FALSE)

# Run over samples and create samplereport
for (i in 1:nrow(samplesData)){
  render(markdownfile, output_format = "html_document", params = list(sample = samplesData$V1[i], panell = samplesData$V2[i]), output_file = file.path(reportsDir, paste0(samplesData$V1[i],"report.html")))
}





