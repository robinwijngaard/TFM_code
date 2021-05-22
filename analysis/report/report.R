library(gggenes)
library(ggplot2)
library(rmarkdown)
library(DT)

reportDir <- "~/Dropbox/Master UOC/TFM/Report"
setwd(reportDir)

cnvFounds_exomedepth <- read.delim("cnvFounds_exomedepth.txt")[, c(7,5,6,3,13,14,9,12)]
cnvNames <- colnames(cnvFounds_exomedepth)

cnvFounds_exomedepth$sample <- sub("X", "", cnvFounds_exomedepth$sample)
cnvFounds_exomedepth$sample <- as.character(cnvFounds_exomedepth$sample)

annotatedFile <- read.delim(file.path(reportDir, "annotatedFile.bed"), sep="\t", header = TRUE)  

samplesData <- read.delim("~/Dropbox/Master UOC/TFM/Report/samples.txt", header=FALSE)

for (i in 1:nrow(samplesData)){
  render("reportmrkdown.Rmd", output_format = "html_document", params = list(sample = samplesData$V1[i], panell = samplesData$V2[i]), output_file = file.path(reportDir, paste0(samplesData$V1[i],"report.html")))
}



