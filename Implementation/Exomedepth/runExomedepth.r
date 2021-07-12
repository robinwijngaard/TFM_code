# Runs ExomeDepth over the datasets cofigured at [datasets_params_file]
#USAGE: Rscript runExomeDepth.r [ExomeDepth_params_file] [datasets_params_file]
print(paste("Starting at", startTime <- Sys.time()))
suppressPackageStartupMessages(library(ExomeDepth))
library(methods)
library(R.utils)

# Read args
args <- commandArgs(asValues = TRUE)

## Input files
bamsDir <- args$bamsDir
bedFile <- args$bedFile
fastaFile <- args$fastaFile 

## output dirs
resultDir <- args$resultDir

## Parameters
run <- args$run
readLength <- args$readLength
phi.bins <- as.numeric(args$phi_bins)
transition.probability <- as.numeric(args$transition_probability)
expected.CNV.length <- as.numeric(args$expected_CNV_length)

# returns gene for position matching ROI in bed data
auxGetGene <- function(bedData, chr, pos){
  gene <- bedData[bedData$chr == chr & (pos >= bedData$start & pos <= bedData$end), "gene"]
  if (length(gene) == 0)
    return(NULL)
  else
    return(gene)
}
  
# Process ExomeDepth alg. body
processExomedepthBody <- function(testCountsDF, controlCountsDF, countsDef, phi.bins, transition.probability, expected.CNV.length){
  all <- data.frame()
  
  for (i in 1:ncol(testCountsDF)) {
    # obtain sampleName: remove first char "X"
    sampleName = colnames(testCountsDF)[i]
    sampleName <-  sub("X", "", sampleName)
    
    # build necessary matrix
    controlCounts <- as.matrix(controlCountsDF)[, -i]
    testCounts <- as.matrix(testCountsDF)[, i]
    
    # define final controls set
    references = select.reference.set(test.counts = testCounts,
                                      reference.count = controlCounts,
                                      bin.length = (countsDF$end - countsDF$start) / 1000, 
                                      n.bins.reduced = 10000,
                                      phi.bins = phi.bins)
    controls <- apply(X = as.matrix(controlCounts[, references$reference.choice]), MAR = 1, FUN = sum)

    
    # call cnvs
    all_exons = new('ExomeDepth',
                    test = testCounts,
                    reference = controls,
                    formula = 'cbind(test, reference) ~ 1')
    
    all_exons = CallCNVs(x = all_exons,
                         transition.probability = transition.probability,
                         expected.CNV.length = expected.CNV.length,
                         chromosome = countsDef$chromosome,
                         start = countsDef$start,
                         end = countsDef$end,
                         name = countsDef$exon)
    
    # Save RDS for sample report generation
    saveRDS(all_exons, file = file.path(rdsDir, paste0("all_exons_", sampleName, ".RDS")))
    
    if (nrow(all_exons@CNV.calls) > 0){
      # add sample column
      all_exons@CNV.calls$sample <- sampleName
      
      # add gene column
      for(i in 1:nrow(all_exons@CNV.calls)) {
        row <- all_exons@CNV.calls[i,]
        row$gene <- auxGetGene(bedData, row$chromosome, row$end)
        if (is.null(row$gene))
          stop(paste("Error: gene not found for chromosome", row$chromosome, ", start", row$start, ", end", row$end))
        all_exons@CNV.calls[i, "Gene"] <- row$gene
      }
    }
    
    all <- rbind(all, all_exons@CNV.calls)
  }
  
  return(all)
}

# go over datasets and run ExomeDepth for those which are active
print(paste("Starting ExomeDepth for", run, "dataset", sep=" "))
    
# read bed file
bedData <- read.table(bedFile, sep="\t", stringsAsFactors=FALSE, col.names = (c("chr", "start", "end", "gene")))
    
# build output folder and file
runDir <- file.path(resultDir, run)
rdsDir <- file.path(runDir, "RDSfiles")

if(!dir.exists(runDir)) {dir.create(runDir)}
if(!dir.exists(rdsDir)) {dir.create(rdsDir)}

# output Files
outputFile <- file.path(runDir, "all_cnv_calls.txt")
    
# Do pre-calc part of the algorithm
# read bam counts   
bamFiles <- list.files(bamsDir, "*.bam$", full.names=T)
counts <- getBamCounts(bed.file = bedFile, 
                       bam.files = bamFiles,
                       read.width = readLength,
                       referenceFasta = fastaFile)
countsDF  <- as.data.frame(counts)
names(countsDF) <- gsub(".bam", "", names(countsDF)) # remove .bam from sample name
      
# Fix rare bug: sometimes X is added at the beginning of input samples
parts <- strsplit(bamFiles[1], "/")
original <- gsub(".bam", "", parts[[1]][length(parts[[1]])])
if (substr(original, 1, 1) != "X" && substr(names(countsDF)[6], 1, 1) == "X")
names(countsDF)[6:ncol(countsDF)] <- substring(names(countsDF)[6:ncol(countsDF)] , 2)
      
# Process each sample
all <- data.frame()
all <- processExomedepthBody(testCountsDF = countsDF[6:ncol(countsDF)],
                              controlCountsDF = countsDF[6:ncol(countsDF)],
                              countsDef = countsDF[1:5],
                              phi.bins = phi.bins,
                              transition.probability = transition.probability,
                              expected.CNV.length = expected.CNV.length)

# write results
write.table(all, file = outputFile, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
    
df <- read.table(outputFile,sep="\t", stringsAsFactors=FALSE, header = TRUE) # REMOVE

print(paste("ExomeDepth for", run, "dataset finished", sep=" "))

print(paste("Finishing at", endTime <- Sys.time()))
cat("\nElapsed time:")
print(endTime - startTime)


# Generate sample template report
print(paste("Start generating sample template reports", startTime <- Sys.time()))




