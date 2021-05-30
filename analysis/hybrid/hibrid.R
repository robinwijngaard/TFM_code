library(caret)

# Load dirs
analysisDir <- "~/Dropbox/Master_UOC/TFM/TFM_code/analysis"
hybridDir <- file.path(analysisDir, "hybrid")
bedDir <- file.path(analysisDir, "bedfiles")
tempDir <- file.path(hybridDir, "temp")
resultDir <- file.path(hybridDir, "results")
evalDir <- file.path(analysisDir, "evaluation")
resDir <- file.path(evalDir, "results")

# Load result file
hybridFile <- file.path(resDir, "allhybridresults.txt")
hybridData <- read.table(hybridFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE)[, -c(14, 16)]

# Adapt input format
## There is one ROI for DeCON (FP) where both DEL and DUP is predicted. This ROI is considered as normal for the hybrid model
NNData <- hybridData[, c(4, 7, 8, 10, 12:15)]
NNData$cnvkit5[NNData$cnvkit5 == -1] <- "Deletion"
NNData$cnvkit5[NNData$cnvkit5 == 1] <- "Duplication"
NNData$cnvkit5[NNData$cnvkit5 == 0] <- "Normal"
NNData$cnvkit5 <- as.factor(NNData$cnvkit5)

NNData$panelcn[NNData$panelcn == -1] <- "Deletion"
NNData$panelcn[NNData$panelcn == 1] <- "Duplication"
NNData$panelcn[NNData$panelcn == 0] <- "Normal"
NNData$panelcn <- as.factor(NNData$panelcn)

NNData$exomedepth[NNData$exomedepth == -1] <- "Deletion"
NNData$exomedepth[NNData$exomedepth == 1] <- "Duplication"
NNData$exomedepth[NNData$exomedepth == 0] <- "Normal"
NNData$exomedepth <- as.factor(NNData$exomedepth)

NNData$convading[NNData$convading == -1] <- "Deletion"
NNData$convading[NNData$convading == 1] <- "Duplication"
NNData$convading[NNData$convading == 0] <- "Normal"
NNData$convading <- as.factor(NNData$convading)

NNData$cnv <- as.factor(NNData$cnv)
NNData$cnv_type[NNData$cnv == "Normal"] <- "Normal"
NNData$cnv_type <- as.factor(NNData$cnv_type)

## Artificial Neural Network

### Entrenar el modelo

NNData_1 <- NNData[, c(2,4:8)]

set.seed(123)

s <- createDataPartition(y = NNData_1$cnv, p = 0.67, list = FALSE)
train.data <- NNData_1[s, ]
test.data <- NNData_1[-s, ]

ctrl <- trainControl(classProbs = TRUE, method = "cv", number =  10, summaryFunction = twoClassSummary)
grid <- expand.grid(size = c(6, 8, 10, 12), decay = 0)

ann.mdl <- train(cnv ~ ., data = train.data, method = "nnet", tuneGrid = grid, trControl = ctrl, 
                 preProc = c("center", "scale"),trace = FALSE, metric = "Sens", maximize = TRUE)

ann.mdl
plot(ann.mdl)
ann <- ann.mdl$finalModel$tuneValue
ann.pred <- predict(ann.mdl, newdata = test.data[, -1])
ann.cm <- confusionMatrix(ann.pred, test.data$cnv)
ann.cm

ann.pred <- predict(ann.mdl, newdata = NNData_1[, -1])
ann.cm <- confusionMatrix(ann.pred, NNData_1$cnv)
ann.cm

## Random Forest

set.seed(123)

ctrl <- trainControl(classProbs = TRUE, method = "cv", number =  10, summaryFunction = twoClassSummary)
grid <- expand.grid(mtry = c(50, 100))

rf.mdl <- train(cnv ~ ., data = train.data, method = "rf", tuneGrid = grid, trControl = ctrl, trace = FALSE, preProcess = c("center", "scale"),
                metric = "Sens", maximize = TRUE)

rf.mdl
