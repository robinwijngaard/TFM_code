setwd("~/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/decon-dataset1")
grPositives <- readRDS("grPositives.rds")


setwd("~/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/panelcn-dataset1")
grPositives2 <- readRDS("grPositives.rds")



outputFolder <- "/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/decon-datasetsingle"
setwd("/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master")

source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

saveResultsFileToGR(outputFolder, "calls_all.txt", chrColumn = "Chromosome", startColumn = "Start", endColumn = "End")




resultDir <- "/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/output/manta-datasetsinglefilter/results/variants"
svTools <- "/home/robin/Documents/Project/TFM_code/Files/svtools-master/"

system(paste("gzip -d", file.path(resultDir, "diploidSV.vcf.gz")))
system(paste("python", file.path(svTools, "vcfToBedpe"), "-i diploidSV.vcf -o bedpeSV"))

system(paste("vcftools --vcf", file.path(resultDir, "diploidSV.vcf"), "--keep-only-indels > indels.vcf")) 
