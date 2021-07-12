#!/bin/bash

outputDir=/home/robin/Documents/Project/Files/Template #outputDir
annotatedFile=/home/robin/Documents/Project/Files/Template/annotatedFile.bed #annotatedfile (output annotate.R)
panellsFile=/home/robin/Documents/Project/Files/Template/gens_nm_panell.xlsx #file idnicating genes per panell and NM
samplesFile=/home/robin/Documents/Project/Files/Template/samples.txt #file indicating sample and panell association
RDSfolder=/home/robin/Documents/Project/Files/Template/RDSfiles #folder containing output RDS files from the ExomeDepth algorithm
markdownfile=/home/robin/Documents/Project/TFM_code/Template/report.Rmd #rmarkdown file

# for file in *; do mv "${file}" "${file/X/}"; done   eliminate X for .RDS file names

Rscript report.R --outputDir $outputDir --annotatedFile $annotatedFile --panellsFile $panellsFile --samplesFile $samplesFile --RDSfolder $RDSfolder --markdownfile $markdownfile