# Make a single CNVfounds file from all runs for every algorithm
library(R.utils)
args=commandArgs(asValues = TRUE)

cnvFounds <- args$cnvFounds

convadingDir <- file.path(cnvFounds, "convading")
exomeDir <- file.path(cnvFounds, "exomedepth")
panelcnDir <- file.path(cnvFounds, "panelcn")
#cnvkitDir <- file.path(cnvFounds, "cnvkit5")

# Run over convading
convading_all <- data.frame()
subdirs <- list.files(convadingDir)
for(i in subdirs){
  cnvfounds <- read.table(file.path(convadingDir, i, "cnvFounds.txt"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)
  convading_all <- rbind(convading_all, cnvfounds)
}

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

write.table(panelcn_all, file.path(cnvFounds, "cnvFounds_panelcn.txt"), sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)  

# run over cnvkit5