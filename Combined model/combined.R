library(caret)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)

## There is one ROI for DeCON (FP) where both DEL and DUP is predicted. This ROI is considered as normal for the hybrid model

# Load dirs
analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
hybridDir <- file.path(analysisDir, "hybrid")
bedDir <- file.path(analysisDir, "bedfiles")
tempDir <- file.path(hybridDir, "temp")
resultDir <- file.path(hybridDir, "results")
evalDir <- file.path(analysisDir, "evaluation")
resDir <- file.path(evalDir, "results")
clinicDir <-  "~/Dropbox/Master_UOC/TFM/Clinic/evaluation/results"

# Load result file
hybridFile <- file.path(resDir, "allhybridresults.txt")
hybridData <- read.table(hybridFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)[, -c(14, 16)]

# Obtain GC content
ranges <- cbind(hybridData[, 1:3], rep(".", nrow(hybridData)))
colnames(ranges) <- c("chr", "start", "end", "strand")   
ranges$chr <- paste0("chr", ranges$chr)

#Build a GRanges from your matrix
ranges <- toGRanges(data.frame(ranges))

#Get the sequences and compute the GC content
freqs <- alphabetFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, ranges))
gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)
hybridData$gc <- gc

# Transform data to numeric levels (0: normal, 1:dup, 2: del)
modelData <- hybridData[, c(7,10,12:16)]
modelData$cnvkit5[modelData$cnvkit5 == -1] <- 2
modelData$panelcn[modelData$panelcn == -1] <- 2
modelData$exomedepth[modelData$exomedepth == -1] <- 2
modelData$convading[modelData$convading == -1] <- 2
modelData$cnv <- as.factor(modelData$cnv)

# Separate in train and test set
set.seed(123)
s <- createDataPartition(y = modelData$cnv, p = 0.67, list = FALSE)
train.data <- modelData[s, ]
test.data <- modelData[-s, ]

# Set ctrl: 10-fold CV
ctrl <- trainControl(classProbs = TRUE, method = "cv", number =  10, summaryFunction = twoClassSummary)

# Model 1: artificial neural network
ann.mdl <- train(cnv ~ ., data = train.data, method = "nnet", trControl = ctrl, tuneLength = 4,
                 preProc = c("range"),trace = FALSE, metric = "Sens", maximize = TRUE)


ann.mdl
ann.pred <- predict(ann.mdl, newdata = test.data[, -1])
ann.cm <- confusionMatrix(ann.pred, test.data$cnv)
ann.cm

# Model 2: random forest
rf.mdl <- train(cnv ~ ., data = train.data, method = "rf",  trControl = ctrl, tuneLength = 4,
                trace = FALSE, preProcess = c("center", "scale"), metric = "Sens", maximize = TRUE)

rf.mdl
rf.pred <- predict(rf.mdl, newdata = test.data[, -1])
(rf.cm <- confusionMatrix(rf.pred, test.data$cnv))

# Model 3: naive-bayes
grid <- expand.grid(.adjust = 1,
                    .usekernel = c(FALSE, TRUE),
                    .laplace = c(0, 1))
nb.mdl <- train(cnv ~ ., data = train.data, method = "naive_bayes", trControl = ctrl, tuneGrid = grid,
                metric = "Sens", maximize = TRUE)

nb.mdl
nb.pred <- predict(nb.mdl, newdata = test.data[, -1])
(nb.cm <- confusionMatrix(nb.pred, test.data$cnv))

# Model 4: SVM radial
svmR.mdl <- train(cnv ~ ., data = train.data, method = "svmRadial", trControl = ctrl, tuneLength = 4,
                  preProc = c("center", "scale"), trace = FALSE, metric = "Sens", maximize = TRUE)
svmR.mdl

svm.pred <- predict(svmR.mdl, newdata = test.data[, -1])
(svm.cm <- confusionMatrix(svm.pred, test.data$cnv))




#### On clinic dataste
clinicFile <- file.path(clinicDir, "clinicMLresults.txt")
clinicData <- read.table(clinicFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)

# Obtain GC content
ranges <- cbind(clinicData[, 1:3], rep(".", nrow(clinicData)))
colnames(ranges) <- c("chr", "start", "end", "strand")   
ranges$chr <- paste0("chr", ranges$chr)

#Build a GRanges from your matrix
ranges <- toGRanges(data.frame(ranges))

#Get the sequences and compute the GC content
freqs <- alphabetFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, ranges))
gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)
clinicData$gc <- gc




ann.pred <- predict(ann.mdl, newdata = test.data[, -1])
ann.cm <- confusionMatrix(ann.pred, test.data$cnv)
ann.cm

rf.pred <- predict(rf.mdl, newdata = test.data[, -1])
(rf.cm <- confusionMatrix(rf.pred, test.data$cnv))

nb.pred <- predict(nb.mdl, newdata = test.data[, -1])
(nb.cm <- confusionMatrix(nb.pred, test.data$cnv))

svm.pred <- predict(svmR.mdl, newdata = test.data[, -1])
(svm.cm <- confusionMatrix(svm.pred, test.data$cnv))



