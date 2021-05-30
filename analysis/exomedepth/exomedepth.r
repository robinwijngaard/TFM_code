library(ggplot2)
library(dplyr)

analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
exomeDir <- file.path(analysisDir, "exomedepth")
bedDir <- file.path(analysisDir, "bedfiles")
evalDir <- file.path(analysisDir, "evaluation")
resDir <- file.path(evalDir, "results")
tempDir <- file.path(exomeDir, "temp")
graphsDir <- file.path(exomeDir, "graphs")
resultsDir <- file.path(exomeDir, "results")

# Load result file
resFile <- file.path(resDir, "allresults.txt")
resData <- read.table(resFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
resData <- resData[, c(1:11,15)]

# Extract TP and FP
TP_exome <- subset(resData, resData$exomedepth == 1 & resData$cnv == "ExonCNV")
FP_exome <- subset(resData, resData$exomedepth == 1 & resData$cnv != "ExonCNV")

# Write tables
write.table(TP_exome, file.path(tempDir, "TP_exome.bed"), sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)  
write.table(FP_exome, file.path(tempDir, "FP_exome.bed"), sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)  

# Load cnvfound file
cnvFile <- file.path(exomeDir, "cnvFounds_exomedepth.txt")
cnvData <- read.table(cnvFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)
cnvData$sample <- sub("X", "", cnvData$sample)

# Order columns to adapt to bed format
cnvData <- cnvData[, c(7,5,6,3,4,8:14)]

# Write table
write.table(cnvData, file.path(tempDir, "cnvData.bed"), sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)  

# Intersect
samples <- unique(resData$sampleID)

system(paste("bedtools intersect -wao -a", file.path(tempDir, "TP_exome.bed"), "-b", file.path(tempDir, "cnvData.bed"), ">", file.path(tempDir, "TP_exome_annotated.bed")))
system(paste("bedtools intersect -wao -a", file.path(tempDir, "FP_exome.bed"), "-b", file.path(tempDir, "cnvData.bed"), ">", file.path(tempDir, "FP_exome_annotated.bed")))

# Read table
TP_exome_annotated <- read.delim(file.path(tempDir, "TP_exome_annotated.bed"), header = FALSE)
FP_exome_annotated <- read.delim(file.path(tempDir, "FP_exome_annotated.bed"), header = FALSE)

# Add colnames
colnames(TP_exome_annotated) <- colnames(FP_exome_annotated) <- c(colnames(resData), colnames(cnvData), "overlap")

# Subset only sampleID cnvData == sampleID resData
TP_exome_annotated <- TP_exome_annotated[which(TP_exome_annotated$sampleID == TP_exome_annotated$sample), ]
FP_exome_annotated <- FP_exome_annotated[which(FP_exome_annotated$sampleID == FP_exome_annotated$sample), ]

# Compare read.ratio and BF between sets
plot.table <- data.frame(BF = c(TP_exome_annotated$BF, FP_exome_annotated$BF), reads.ratio = c(TP_exome_annotated$reads.ratio, FP_exome_annotated$reads.ratio), Pos = c(rep("TP", nrow(TP_exome_annotated)), rep("FP", nrow(FP_exome_annotated))), Exon = c(TP_exome_annotated$exon_type, FP_exome_annotated$exon_type))
plot.table$Group[plot.table$Pos == "TP" & plot.table$Exon == "Single"] <- "TP single"
plot.table$Group[plot.table$Pos == "TP" & plot.table$Exon == "Multi"] <- "TP multi"
plot.table$Group[plot.table$Pos == "FP"] <- "FP"

plot.table$Group <- factor(plot.table$Group, levels = c("TP single", "TP multi", "FP"))

BF_comp <- plot.table %>% group_by(Group) %>% summarize(mean = mean(BF), 
                                                                    min = min(BF), max = max(BF),
                                                                    median = median(BF),
                                                                    IQR25 = quantile(BF, 0.25), 
                                                                    IQR75 = quantile(BF, 0.75))


reads.ratio_comp <- plot.table %>% group_by(Group) %>% summarize(mean = mean(reads.ratio), 
                                                                    min = min(reads.ratio), max = max(reads.ratio),
                                                                    median = median(reads.ratio),
                                                                    IQR25 = quantile(reads.ratio, 0.25), 
                                                                    IQR75 = quantile(reads.ratio, 0.75))

write.table(BF_comp, file.path(resultsDir, "BF_comp.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(reads.ratio_comp, file.path(resultsDir, "reads.ratio_comp.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  

tiff(file.path(graphsDir, "BF_comp.tiff"), units="in", width=7, height=4, res=150)

ggplot(plot.table, aes(x=BF, color=Group, fill=Group)) +
  geom_histogram(aes(y=..density..), bins=50, position="identity", alpha=0.2)+
  geom_density(alpha=0.2)+
  facet_grid(Group ~ .) +
  labs(x="Bayes Factor")+
  theme_classic()+
  theme(legend.position = "none")

dev.off()

tiff(file.path(graphsDir, "reads.ratio_comp.tiff"), units="in", width=7, height=5, res=150)

ggplot(plot.table, aes(x=reads.ratio, color=Group, fill=Group)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2)+
  geom_density(alpha=0.2)+
  facet_grid(Group ~ .) +
  labs(title="Read-ratio between groups",x="read.ratio", y = "Density")+
  theme_classic()

dev.off()

# TP BF by single and multi
TP_exon.type <- TP_exome_annotated[, -c(14, 15)] %>% group_by(exon_type) %>% summarize(mean = mean(BF), 
                                                                   min = min(BF), max = max(BF),
                                                                   median = median(BF),
                                                                   IQR25 = quantile(BF, 0.25), 
                                                                   IQR75 = quantile(BF, 0.75))

FP_exon.type <- FP_exome_annotated[, -c(14, 15)] %>% group_by(exon_type) %>% summarize(mean = mean(BF), 
                                                                                       min = min(BF), max = max(BF),
                                                                                       median = median(BF),
                                                                                       IQR25 = quantile(BF, 0.25), 
                                                                                       IQR75 = quantile(BF, 0.75))

write.table(TP_exon.type, file.path(resultsDir, "TP_exontype.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  
write.table(FP_exon.type, file.path(resultsDir, "FP_exontype.txt"), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)  







