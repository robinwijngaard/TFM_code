# Runs Manta over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runManta.r [Manta_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ExomeDepth))
library(methods)
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

# Convert to deletion, duplication
auxConvert <- function(x) {
  if (x == "DEL") return("deletion") 
  else if (x == "DUP") return("duplication")
  else return("")
}

# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  mantaParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  mantaParamsFile <- "algorithms/manta/mantaParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
}

#Load the parameters file  
params <- yaml.load_file(mantaParamsFile)
datasets <- yaml.load_file(datasetsParamsFile)

# extract manta params
mantaFolder <- file.path(params$mantaFolder)
svTools <- file.path(params$svTools)
print(paste("Params for this execution:", list(params)))

# create configManta.py.ini file

configManta.py.ini <- c("[manta]",
                        paste("referenceFasta =", params$referenceFasta),
                        paste("minCandidateVariantSize =", params$minCandidateVariantSize),
                        paste("rnaMinCandidateVariantSize =", params$rnaMinCandidateVariantSize),
                        paste("minEdgeObservations =", params$minEdgeObservations),
                        paste("graphNodeMaxEdgeCount =", params$graphNodeMaxEdgeCount),
                        paste("minCandidateSpanningCount =", params$minCandidateSpanningCount),
                        paste("minScoredVariantSize =", params$minScoredVariantSize),
                        paste("minDiploidVariantScore =", params$minDiploidVariantScore),
                        paste("minPassDiploidVariantScore =", params$minPassDiploidVariantScore),
                        paste("minPassDiploidGTScore =", params$minPassDiploidGTScore),
                        paste("minSomaticScore =", params$minSomaticScore),
                        paste("minPassSomaticScore =", params$minPassSomaticScore),
                        paste("enableRemoteReadRetrievalForInsertionsInGermlineCallingModes =", params$enableRemoteReadRetrievalForInsertionsInGermlineCallingModes),
                        paste("enableRemoteReadRetrievalForInsertionsInCancerCallingModes =", params$enableRemoteReadRetrievalForInsertionsInCancerCallingModes),
                        paste("useOverlapPairEvidence =", params$useOverlapPairEvidence),
                        paste("enableEvidenceSignalFilter =", params$enableEvidenceSignalFilter))

configManta.py.ini <- as.data.frame(configManta.py.ini)

# Execute manta 
processMantaBody <- function(launchFile, bamList, fastaFile, subsetFolder, bedFile,configFile, params){
  system(paste0("python ", launchFile," ", bamList, " --referenceFasta=", fastaFile, " --config=", configFile, " --exome", " --runDir=", subsetFolder))
  system(paste("ulimit -n 75000;", paste("python", file.path(subsetFolder, "runWorkflow.py"))))
}

#  " --callRegions=", bedFile,

# go over datasets and run manta
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if(dataset$include){
    print(paste("Starting manta for", name, "dataset", sep=" "))
    
    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- paste0(file.path(dataset$bed_file), ".gz")
    fastaFile <- file.path(dataset$fasta_file)
    
    # create output folder
    outputFolder <- file.path(getwd(), "output", paste0("manta-", name))  
    dir.create(outputFolder)
    
    # save configuration file
    configFile <- file.path(outputFolder, "configManta.py.ini")
    write.table(configManta.py.ini, configFile, sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
    
    # build launch file path
    launchFile <- file.path(mantaFolder, "configManta.py")
    
    # List bam files
    bamFiles <- list.files(bamsDir)
    bamFiles <- bamFiles[!grepl(".bai", bamFiles)]
    
    # Run Manta in group of 6 samples
    sets <- length(bamFiles) / 6

    #cnvfound file
    cnvFounds <- data.frame(matrix(ncol = 7, nrow = 0)) 
    colnames(cnvFounds) <- c("Sample", "Chr", "Start", "End", "CNV.type", "Filter", "SampleFilter")
    
    # Loops over subsets 
    for(i in 1:sets){
      
      #Create subsetFolder
      subsetFolder <- file.path(outputFolder, paste0("set", i))
      dir.create(subsetFolder)
      
      # Obtain Bams of set and run manta
      inici <- (i - 1)*6 + 1
      final <- i*6
      
      bamList <- c()
      for (bam in bamFiles[inici:final]){
        bamList <- c(bamList, paste0("--bam=", file.path(bamsDir, bam)))
      }
      bamList <- paste(bamList, collapse = " ")
      
      # Process Manta
      processMantaBody(launchFile, bamList, fastaFile, subsetFolder, bedFile,configFile, params)
      
      # Decompress result VCF file
      resultFolder <- file.path(subsetFolder, "results", "variants")
      system(paste("gzip -d", file.path(resultFolder, "diploidSV.vcf.gz")))
      
      # Transform VCF to bed
      system(paste("python", file.path(svTools, "vcfToBedpe"), "-i",file.path(resultFolder, "diploidSV.vcf"), ">", file.path(resultFolder, "diploidSV.bedpe")))
      
      # Read bed Result File
      resultFile <- file.path(resultFolder, "diploidSV.bedpe")
      resultData <- read.table(resultFile, sep="\t", stringsAsFactors=FALSE, comment.char = "", header = TRUE)
      
      # Only DUP and DEL and filter useless columns
      resultData <- subset(resultData, resultData$TYPE == "DUP" | resultData$TYPE == "DEL")
      ncols <- ncol(resultData)
      resultData <- resultData[, c(1, 2, 5, 11, 12, 15:ncols)]
      
      # Add 1 to start.A and start.B 
      resultData$START_A <- resultData$START_A + 1
      resultData$START_B <- resultData$START_B + 1
      
      # Edit col names
      sampleNames <- sub(".bam", "", bamFiles[inici:final])
      colnames(resultData) <- c("Chr", "Start", "End", "CNV.type", "Filter", sampleNames)
      
      # Add CNV to cnvFounds per sample
      for (sample in sampleNames){
        sampleData <- resultData[, c(1:5)]
        sampleData <- cbind(sampleData, resultData[, sample])
        colnames(sampleData)[6] <- sample
        
        # Select positive cases
        positivePattern <- c("0/1", "1/1")
        positiveCases <- grep(paste(positivePattern, collapse = "|"), sampleData[, 6])
        sampleData <- sampleData[positiveCases, ]
        
        # Add lines to CNVfound
        sampleID <- rep(sample, nrow(sampleData))
        cnvSample <- cbind(sampleID, sampleData)
        
        # Rename colnames
        colnames(cnvSample) <- c("Sample", "Chr", "Start", "End", "CNV.type", "Filter", "SampleFilter")
        
        #Add CNVs to CNVfound files
        cnvFounds <- rbind(cnvFounds, cnvSample)
      }
      
    } 
    
    # Convert CNVtype
    cnvFounds$CNV.type <- lapply(cnvFounds$CNV.type, function(x) sapply(x, auxConvert))

    # Write result table
    cnvFounds$CNV.type <- as.character(cnvFounds$CNV.type)
    write.table(cnvFounds, file.path(outputFolder, "cnvFounds.txt"), sep="\t", row.names=FALSE, quote = FALSE)  
    
    # Save results in GRanges format
    message("Saving GenomicRanges results")
    saveResultsFileToGR(outputFolder, "cnvFounds.txt")
    
    print(paste("manta for", name, "dataset finished", sep=" "))
    print(paste("Finishing at", Sys.time()))
    cat("\n\n\n")
  }
}


print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)
