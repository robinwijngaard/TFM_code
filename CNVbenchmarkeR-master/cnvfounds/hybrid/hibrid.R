analysisDir <- "/home/robin/Documents/Project/TFM_code/CNVbenchmarkeR-master/cnvfounds"
allDir <- file.path(analysisDir, "all")
singleDir <- file.path(analysisDir, "single")
clinicDir <- file.path(analysisDir, "clinic")
tempDir <- file.path(analysisDir, "temp")
resultDir <- file.path(analysisDir, "results")
hybridDir <- file.path(analysisDir, "hybrid")


allData <- read.table(file.path(resultDir, "allData.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
clinicData <- read.table(file.path(resultDir, "clinicData.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
singleData <-  read.table(file.path(resultDir, "singleData.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)

allNormal <- read.table(file.path(resultDir, "allNormal.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
clinicNormal <- read.table(file.path(resultDir, "clinicNormal.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
singleNormal <-  read.table(file.path(resultDir, "singleNormal.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)

All4 <- rbind(allData[, -c(6,8, 12, 14)], setNames(allNormal[, -c(10, 12)], names(allData[-c(6, 8,12, 14)])))



# VennDiagram
setwd(hybridDir)
library("VennDiagram")

## TP
tp_venn <- All4[1:296, 7:11]
tp_venn$ID <- as.character(tp_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(tp_venn)[i]
  tp_venn_algo <- tp_venn[, c(1, i)]
  tp_venn_algo <- subset(tp_venn_algo, tp_venn_algo[,2] != 0)
  venn_list[[algo]] <- tp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tp.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TN
tn_venn <- All4[297:27347, 7:11]
tn_venn$ID <- as.character(tn_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(tn_venn)[i]
  tn_venn_algo <- tn_venn[, c(1, i)]
  tn_venn_algo <- subset(tn_venn_algo, tn_venn_algo[,2] == 0)
  venn_list[[algo]] <- tn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_tn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FN
fn_venn <- All4[1:296, 7:11]
fn_venn$ID <- as.character(fn_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(fn_venn)[i]
  fn_venn_algo <- fn_venn[, c(1, i)]
  fn_venn_algo <- subset(fn_venn_algo, fn_venn_algo[,2] == 0)
  venn_list[[algo]] <- fn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FP
fp_venn <- All4[297:27347, 7:11]
fp_venn$ID <- as.character(fp_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(fp_venn)[i]
  fp_venn_algo <- fp_venn[, c(1, i)]
  fp_venn_algo <- subset(fp_venn_algo, fp_venn_algo[,2] != 0)
  venn_list[[algo]] <- fp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_fp.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)




Single4 <- rbind(singleData[, -c(6,8, 12, 14)], setNames(singleNormal[, -c(10, 12)], names(singleData[-c(6, 8,12, 14)])))



# VennDiagram
setwd(hybridDir)
library("VennDiagram")

## TP
tp_venn <- Single4[1:24, 7:11]
tp_venn$ID <- as.character(tp_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(tp_venn)[i]
  tp_venn_algo <- tp_venn[, c(1, i)]
  tp_venn_algo <- subset(tp_venn_algo, tp_venn_algo[,2] != 0)
  venn_list[[algo]] <- tp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_single_tp.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## TN
tn_venn <- Single4[25:17119, 7:11]
tn_venn$ID <- as.character(tn_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(tn_venn)[i]
  tn_venn_algo <- tn_venn[, c(1, i)]
  tn_venn_algo <- subset(tn_venn_algo, tn_venn_algo[,2] == 0)
  venn_list[[algo]] <- tn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_single_tn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FN
fn_venn <- Single4[1:24, 7:11]
fn_venn$ID <- as.character(fn_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(fn_venn)[i]
  fn_venn_algo <- fn_venn[, c(1, i)]
  fn_venn_algo <- subset(fn_venn_algo, fn_venn_algo[,2] == 0)
  venn_list[[algo]] <- fn_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_single_fn.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

## FP
fp_venn <- Single4[25:17119, 7:11]
fp_venn$ID <- as.character(fp_venn$ID)

venn_list <- list()

for(i in 2:5){
  algo <- names(fp_venn)[i]
  fp_venn_algo <- fp_venn[, c(1, i)]
  fp_venn_algo <- subset(fp_venn_algo, fp_venn_algo[,2] != 0)
  venn_list[[algo]] <- fp_venn_algo[,1]
}

venn.diagram(venn_list, filename = "venn_single_fp.tiff", col = "black",fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),alpha = 0.50,cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),cat.cex = 1.5,cat.fontface = "bold",margin = 0.2)

