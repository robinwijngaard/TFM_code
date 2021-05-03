# Runs CNVkit over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runCNVkit.r [cnvkit_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(yaml))
library(methods)
library(dplyr)
source(if (basename(getwd()) == "optimizers") "../utils/utils.r" else "utils/utils.r") # Load utils functions

auxCNname <- function(x) {
  if (x %in% c(0, 1)) return("deletion") 
  else if (x %in% c(3, 4)) return("duplication")
}

# execute cnvkit

processCNVkitBody <- function(launchFile, bedFile, outputFolder,bamsDir, params, fastaFile, refbamsDir){
  # load bams
  bamFiles <- list.files(bamsDir)
  bamFiles <- bamFiles[!grepl(".bai", bamFiles)]
  
  antitarget <- data.frame(matrix(ncol=5, nrow=0))
  colnames(antitarget) <- c("chromosome", "start", "end", "gene", "log2")
  
  # create ref file
  #system(paste("python", launchFile, "batch -m amplicon --normal", file.path(bamsDir, "*.bam"), "--targets", bedFile, "--fasta", fastaFile, "--output-reference", file.path(outputFolder, "reference.cnn"), "--output-dir", outputFolder))
  
  for(bam in bamFiles){
    sampleName <- sub(".bam", "", bam)
    write.table(antitarget, file.path(outputFolder, paste0(sampleName, ".antitargetcoverage.cnn")), sep="\t", row.names=FALSE, quote = FALSE)
    system(paste("python", launchFile, "coverage", file.path(bamsDir, bam), bedFile, "-o", file.path(outputFolder, paste0(sampleName, ".targetcoverage.cnn"))))
  }
  
  coverageFiles <- list.files(outputFolder, pattern = ".targetcoverage.cnn")
  coverageFilesPath <- file.path(outputFolder, coverageFiles)
  coverageFilesList <- paste(coverageFilesPath, collapse = " ")
  
  system(paste("python", launchFile, "reference", coverageFilesList, "-f", fastaFile, "--no-edge -o", file.path(outputFolder, "reference.cnn")))
  refFile <- file.path(outputFolder, "reference.cnn")
  
  # batch
  system(paste("python", launchFile, "batch -m amplicon", file.path(bamsDir, "*.bam"), "-r", file.path(outputFolder, "reference.cnn"), "--segment-method hmm --output-dir", outputFolder))
  
  # call on bintest
  cnsFiles <- list.files(outputFolder, pattern = ".bintest.cns")
  for (cnsFile in cnsFiles){
    sampleName <- sub(".bintest.cns", "", cnsFile)
    system(paste("python", launchFile, "call", file.path(outputFolder, cnsFile), "-m threshold -o", file.path(outputFolder, paste0(sampleName, ".bintest.call.cns"))))
  } 
}
  
# Read args
args <- commandArgs(TRUE)
print(args)
if(length(args)>0) {
  cnvkitParamsFile <- args[1]
  datasetsParamsFile <- args[2]
} else {
  cnvkitParamsFile <- "algorithms/cnvkit/cnvkitParams.yaml"
  datasetsParamsFile <- "datasets.yaml"
}

#Load the parameters file  
params <- yaml.load_file(cnvkitParamsFile) 
datasets <- yaml.load_file(datasetsParamsFile)

# extract cnvkit params
cnvkitFolder <- file.path(params$cnvkitFolder)
print(paste("Params for this execution:", list(params)))

# go over datasets and run cnvkit
for (name in names(datasets)) {
  dataset <- datasets[[name]]
  if(dataset$include){
    print(paste("Starting cnvkit for", name, "dataset", sep=" "))
    
    # extract fields
    bamsDir <- file.path(dataset$bams_dir)
    bedFile <- file.path(dataset$bed_file)
    bedData <- read.table(bedFile, sep = "\t", stringsAsFactors=FALSE)
    colnames(bedData) <- c("Chr", "Start", "End", "Gene")
    fastaFile <- file.path(dataset$fasta_file)
    refbamsDir <- file.path(params$refbamsDir)
    
    # create output folder
    outputFolder <- file.path(getwd(), "output", paste0("cnvkit2-", name))  
    dir.create(outputFolder)

    # build launch file path
    launchFile <- file.path(cnvkitFolder, "cnvkit.py")
    launchFile2 <- file.path(cnvkitFolder, "cnvlib", "bintest.py")
    
    # Execute CNVkit
    processCNVkitBody(launchFile, bedFile, outputFolder,bamsDir, params, fastaFile, refbamsDir)
    
    # Export results
    calls <- list.files(outputFolder, pattern = ".bintest.call.cns")
    
    # Create main files
    failedROIs <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(failedROIs) <- c("SampleID", "Chr", "Start", "End", "Gene")
    
    cnvFounds <- data.frame(matrix(ncol = 8, nrow = 0)) 
    colnames(cnvFounds) <- c("Sample","Gene", "Chr", "Start", "End", "log2", "CN", "CNV.type")

    # Obtain failedcalls and cnv
    for (call in calls){
      sample <- sub(".bintest.call.cns", "", call)
      callPath <- file.path(outputFolder, call)
      callData <-  read.table(callPath, sep="\t", stringsAsFactors=FALSE, header = TRUE)
      
      # delete normal (cn = 2)
      cnvCall <- subset(callData, callData$cn != 2)
      cnvCall$CNV.type <- lapply(cnvCall$cn, function(x) sapply(x, auxCNname))
      
      #Add sampleID column
      Sample <- rep(sample, nrow(cnvCall))
      cnvCall <- cbind(Sample, cnvCall)
      
      #Reorder and delete non-necessary colums
      cnvCall <- cnvCall[, c(1, 5, 2, 3, 4, 6, 7, 11)]
      colnames(cnvCall) <- c("Sample","Gene", "Chr", "Start", "End", "log2", "CN", "CNV.type") 
      
      cnvFounds <- rbind(cnvFounds, cnvCall)
    
      #failed calls
      failedFile <- paste0(sample, ".cnr")
      failedPath <- file.path(outputFolder, failedFile)
      failedData <-  read.table(failedPath, sep="\t", stringsAsFactors=FALSE, header = TRUE)
      
      failedcall <- failedData[, 1:4]
      colnames(failedcall) <- c("Chr", "Start", "End", "Gene")
      
      failedcall <- anti_join(bedData, failedcall, by=c("Chr", "Start", "End", "Gene"))
      
      # Add sample ID column
      SampleID <- rep(sample, nrow(failedcall))
      failedcall <- cbind(SampleID, failedcall)
      
      # Add to main file
      failedROIs <- rbind(failedROIs, failedcall)
      

    }
    
    #Write failedROIs
    write.table(failedROIs, file.path(outputFolder, "failedROIs.csv"), sep="\t", row.names=FALSE, quote = FALSE)
    
    # Write result table
    cnvFounds$CNV.type <- as.character(cnvFounds$CNV.type)
    write.table(cnvFounds, file.path(outputFolder, "cnvFounds.txt"), sep="\t", row.names=FALSE, quote = FALSE)  
    
    # Save results in GRanges format
    message("Saving GenomicRanges results")
    saveResultsFileToGR(outputFolder, "cnvFounds.txt")
    
    print(paste("cnvkit for", name, "dataset finished", sep=" "))
    print(paste("Finishing at", Sys.time()))
    cat("\n\n\n")
  }
}



print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)
