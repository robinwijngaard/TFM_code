# Make a single CNVfounds file from all runs for every algorithm
library(R.utils)
args=commandArgs(asValues = TRUE)

outputbenchmarker <- args$outputbenchmarker
cnvFounds <- args$cnvFounds

convadingDir <- file.path(outputbenchmarker, "convading")
exomeDir <- file.path(outputbenchmarker, "exomedepth")
panelcnDir <- file.path(outputbenchmarker, "panelcn")
cnvkitDir <- file.path(outputbenchmarker, "cnvkit")

# Run over convading
convading_all <- data.frame()
subdirs <- list.files(convadingDir)
for(i in subdirs){
  cnvfounds <- read.table(file.path(convadingDir, i, "cnvFounds.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  convading_all <- rbind(convading_all, cnvfounds)
}

## Change start and end if end > start
a <- which(convading_all$Start > convading_all$End)
convading_all[a, c("Start", "End")] <- convading_all[a, c("End", "Start")] 

write.table(convading_all, file.path(cnvFounds, "cnvFounds_convading.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  

# Run over exomedepth
exome_all <- data.frame()

subdirs <- list.files(exomeDir)
for(i in subdirs){
  cnvfounds <- read.table(file.path(exomeDir, i, "all_cnv_calls.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  exome_all <- rbind(exome_all, cnvfounds)
}

write.table(exome_all, file.path(cnvFounds, "cnvFounds_exomedepth.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  

# run over panelcnmops
panelcn_all <- data.frame()

subdirs <- list.files(panelcnDir)
for(i in subdirs){
  cnvfounds <- read.table(file.path(panelcnDir, i, "cnvFounds.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  panelcn_all <- rbind(panelcn_all, cnvfounds)
}

## add chr to chromosome column
panelcn_all$Chr <- paste0("chr", panelcn_all$Chr)
write.table(panelcn_all, file.path(cnvFounds, "cnvFounds_panelcn.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  

# run over cnvkit5
cnvkit_all <- data.frame()

subdirs <- list.files(cnvkitDir)
for(i in subdirs){
  cnvfounds <- read.table(file.path(cnvkitDir, i, "cnvFounds.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  cnvkit_all <- rbind(cnvkit_all, cnvfounds)
}

write.table(cnvkit_all, file.path(cnvFounds, "cnvFounds_cnvkit.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  



