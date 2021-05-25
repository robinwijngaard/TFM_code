library(gggenes)
library(ggplot2)
library(rmarkdown)
library(DT)
library(ExomeDepth)
library(dplyr)
library(tiff)
library(patchwork)
library(RColorBrewer)

# mama, colon, melanoma o endocri
analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
reportDir <- file.path(analysisDir, "report")
bedDir <- file.path(analysisDir, "bedfiles")
tempDir <- file.path(reportDir, "temp")
exomedepthDir <- file.path(reportDir, "exomedepth")
outputDir <- file.path(reportDir, "output")
graphsDir <- file.path(reportDir, "graphs")

# Read annotatedFile
annotatedFile <- read.delim(file.path(reportDir, "annotatedFile.bed"), sep="\t", header = TRUE)  

# Read samples and linked panel
samplesData <- read.delim(file.path(reportDir, "samples.txt"), header=FALSE)

# Run over samples and create samplereport
for (i in 1:nrow(samplesData)){
  render(file.path(reportDir, "reportmrkdown.Rmd"), output_format = "html_document", params = list(sample = samplesData$V1[i], panell = samplesData$V2[i]), output_file = file.path(outputDir, paste0(samplesData$V1[i],"report.html")))
}



