# Runs Manta over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runManta.r [Manta_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ExomeDepth))
library(methods)
library(reticulate)
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

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
print(paste("Params for this execution:", list(params)))

# create configManta.py.ini file

configManta.py.ini <- c("[manta]",
                        paste("referenceFasta = ", params$referenceFasta),
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
processMantaBody <- function(launchFile, bamList, fastaFile, outputFolder, bedFile,configFile, params){
  system(paste0("python ", launchFile," ", bamList, " --referenceFasta=", fastaFile, " --config=", configFile, " --exome", " --runDir=", outputFolder))
  system(paste("ulimit -n 75000;", paste("python", file.path(outputFolder, "runWorkflow.py"))))
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
    bamList <- c()
    for (bam in bamFiles){
      bamList <- c(bamList, paste0("--bam=", file.path(bamsDir, bam)))
    }
    bamList <- paste(bamList, collapse = " ")
    
    #bedFile <- file.path(dataset$bed_file)
    #bedData <- read.table(bedFile, sep="\t", stringsAsFactors=FALSE)
    
    #regions <- c()
    #for (i in 1:nrow(bedData)){
    #  regions <- c(regions, paste0("--region=", bedData[i, 1], ":", bedData[i, 2], "-", bedData[i, 3]))
    #}
    
    #regionsList <- paste(regions, collapse = " ")

    
    # Execute Manta
    processMantaBody(launchFile, bamList, fastaFile, outputFolder, bedFile,configFile, params)
    
    # Decompress result VCF file
    
    #system(paste("gzip -d", file.path(outputFolder, "results", "variants", "diploidSV.vcf.gz")))
    
    # read two times the vcf file, first for the columns names, second for the data
    #vcfFile <- readLines(file.path(outputFolder, "results", "variants", "diploidSV.vcf"))
    #vcfData <- read.table(file.path(outputFolder, "results", "variants", "diploidSV.vcf"), stringsAsFactors = FALSE)
    
    # filter for the columns names
    #vcfFile <- vcfFile[-(grep("#CHROM", vcfFile) + 1):-(length(vcfFile))]
    #vcf_names <- unlist(strsplit(vcfFile[length(vcfFile)], "\t"))
    #names(vcfData) <- vcf_names
    
    # Rename sample
    #sampleNames <- sub(".bam", "", bamFiles)
    #nSamples <- length(sampleNames)
    
    # create CNV.Type
    #vcfData$CNV.Type <- ifelse(grepl("DEL", vcfData$ID),"deletion", ifelse(grepl("DUP", vcfData$ID), "duplication", NA))
    #vcfData <- vcfData[complete.cases(vcfData), ]
    
    # separate samples
    #sampleCols <- vcfData[, c(10:(10+nSamples-1))]
    #colnames(sampleCols) <- sampleNames
    #vcfData <- vcfData[, -c(10:(10+nSamples-1))]
  }
}




print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)
