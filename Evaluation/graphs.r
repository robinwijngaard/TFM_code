# Analysis of cnvfound results
library(VennDiagram)
library(MASS)
library(corrplot)
library(R.utils)

args=commandArgs(asValues = TRUE)

# Define dir
outputDir <- args$outputDir
resultFile <- args$resultFile

# make temp dir, results dir en graphs dir
tempDir <- file.path(outputDir, "temp")
resultDir <- file.path(outputDir, "results")
graphsDir <- file.path(outputDir, "graphs")

if(!dir.exists(tempDir)) dir.create(tempDir)
if(!dir.exists(resultDir)) dir.create(resultDir)
if(!dir.exists(graphsDir)) dir.create(graphsDir)

# Load Data
resultData <- read.table(resultFile, sep = "\t",  header=TRUE, stringsAsFactors = FALSE)

# VennDiagram
setwd(graphsDir)
venn_graphs <- resultData[, c(7, 11:15, 17)]
venn_graphs$ID <- as.character(venn_graphs$ID)
names(venn_graphs)[3:7] <- c("CNVkit", "CoNVaDING", "DECoN", "ExomeDepth", "panelcn.MOPS")

make_venn <- function(data, value){
  venn_list <- list()
  for(algorithm in names(data)[3:7]){
    data_algorithm <- data[, c("ID", algorithm)]
    data_algorithm <- subset(data_algorithm, data_algorithm[, algorithm] == value)
    venn_list[[algorithm]] <- data_algorithm$ID
  }
  return(venn_list)
}

## TP
tp_venn <- subset(venn_graphs, venn_graphs$cnv == "ExonCNV")
tp_venn_list <- make_venn(tp_venn, 1)
venn.diagram(tp_venn_list, filename = "tp_venn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TN
tn_venn <- subset(venn_graphs, venn_graphs$cnv == "Normal")
tn_venn_list <- make_venn(tn_venn, 0)
venn.diagram(tn_venn_list, filename = "tn_venn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FN
fn_venn <- subset(venn_graphs, venn_graphs$cnv == "ExonCNV")
fn_venn_list <- make_venn(fn_venn, 0)
venn.diagram(fn_venn_list, filename = "fn_venn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FP
fp_venn <- subset(venn_graphs, venn_graphs$cnv == "Normal")
fp_venn_list <- make_venn(fp_venn, 1)
venn.diagram(fp_venn_list, filename = "fp_venn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

# Multi-scale
mds_data <- resultData[, c(12:15, 17)]
names(mds_data) <- c("CNVkit", "CoNVaDING", "DECoN", "ExomeDepth", "panelcn.MOPS")
d <- dist(t(mds_data), method = "euclidean")
MDS <- isoMDS(d)

tiff("MDS.tiff", units="in", width=6, height=5, res=150)
plot(MDS$points, pch = 16, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cex = 1.5, xlab = "Dimension 1", ylab = "Dimension 2")
labels <- names(mds_data)
legend("topright", inset=.05, legend=labels, col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), pch = 16)
dev.off()

# Corplot
corData <- resultData[, c(12:15, 17)]
names(corData) <- c("CNVkit", "CoNVaDING", "DECoN", "ExomeDepth", "panelcn.MOPS")
corData <- cor(corData)

tiff("corrplot.tiff", units="in", width=6, height=5, res=150)
corrplot(corData, method = "number", cl.lim = c(0,1))
dev.off()

# Cluster
clustData <- resultData[, c(12:15, 17)]
names(clustData) <- c("CNVkit", "CoNVaDING", "DECoN", "ExomeDepth", "panelcn.MOPS")
d <- dist(t(clustData), method = "euclidean")
all.hc <- hclust(d, method = "ward.D2")

tiff("cluster.tiff", units="in", width=6, height=5, res=150)
plot(all.hc)
dev.off()

