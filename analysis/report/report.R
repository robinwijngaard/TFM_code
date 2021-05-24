library(gggenes)
library(ggplot2)
library(rmarkdown)
library(DT)

# mama, colon, melanoma o endocri
analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
reportDir <- file.path(analysisDir, "report")
bedDir <- file.path(analysisDir, "bedfiles")
tempDir <- file.path(reportDir, "temp")
exomedepthDir <- file.path(reportDir, "exomedepth")
outputDir <- file.path(reportDir, "output")

# Read cnvFounds file and list .RDS files
cnvFounds_exomedepth <- read.delim(file.path(exomedepthDir, "cnvFounds_exomedepth.txt"))[, c(7,5,6,3,13,14,9,12)]
cnvNames <- colnames(cnvFounds_exomedepth)

# Correct "X" in sample Name
cnvFounds_exomedepth$sample <- as.character(sub("X", "", cnvFounds_exomedepth$sample))

# Read annotatedFile
annotatedFile <- read.delim(file.path(reportDir, "annotatedFile.bed"), sep="\t", header = TRUE)  

# Read samples and linked panel
samplesData <- read.delim(file.path(reportDir, "samples.txt"), header=FALSE)

# Run over samples and create samplereport
for (i in 1:nrow(samplesData)){
  render(file.path(reportDir, "reportmrkdown.Rmd"), output_format = "html_document", params = list(sample = samplesData$V1[i], panell = samplesData$V2[i]), output_file = file.path(outputDir, paste0(samplesData$V1[i],"report.html")))
}



